#include "GenArch.h"

std::vector<double> GenArch::makeChromosomes(const Param &p) const
{

    std::vector<double> chromends(p.nchrom);

    // Chromosomes all have the same size
    for (size_t chrom = 0u; chrom < p.nchrom; ++chrom)
        chromends[chrom] = (chrom + 1.0) / p.nchrom;

    return chromends;

}

std::vector<size_t> GenArch::makeEncodedTraits(const Param &p) const
{

    std::vector<size_t> encoded(p.nloci);

    // Make an ordered vector of trait indices
    size_t i = 0u;
    for (size_t trait = 0u; trait < 3u; ++trait) {
        for (size_t locus = 0u; locus < p.nvertices[trait]; ++locus) {
            encoded[i] = trait;
            ++i;
        }
    }

    assert(encoded.size() == p.nloci);

    // Shuffle encoded traits randomly
    std::shuffle(encoded.begin(), encoded.end(), rnd::rng);

    assert(encoded.size() == p.nloci);

    std::vector<size_t> nvertices = std::vector<size_t>(3u, 0u);
    for (size_t locus = 0u; locus < p.nloci; ++locus)
        ++nvertices[encoded[locus]];

    for (size_t trait = 0u; trait < 3u; ++trait)
        assert(nvertices[trait] == p.nvertices[trait]);

    return encoded;

}



std::vector<double> GenArch::makeLocations(const Param &p) const
{
    std::vector<double> positions(p.nloci);

    // Locations are sampled from a uniform distribution between 0 and 1
    auto getlocation = rnd::uniform(0.0, 1.0);

    for (size_t locus = 0u; locus < p.nloci; ++locus)
        positions[locus] = getlocation(rnd::rng);

    std::sort(positions.begin(), positions.end());

    for (size_t locus = 1u; locus < p.nloci; ++locus)
        assert(positions[locus] > positions[locus - 1u]);

    return positions;
}


std::vector<double> GenArch::makeEffects(const Param &p) const
{

    if (p.effectshape == 0.0) return std::vector<double>(p.nloci, 0.0);
    if (p.effectscale == 0.0) return std::vector<double>(p.nloci, 1.0);

    std::vector<double> effectsizes(p.nloci);
    // square rooted sum of squares
    std::vector<double> sss = std::vector<double>(3u, 0.0);

    // Effect sizes are sampled from a two-sided Gamma distribution
    auto geteffect = rnd::gamma(p.effectshape, p.effectscale);
    auto isflipped = rnd::bernoulli(0.5);

    for (size_t locus = 0u; locus < p.nloci; ++locus) {

        double effect = geteffect(rnd::rng);
        if (isflipped(rnd::rng)) effect *= -1.0;
        effectsizes[locus] = effect;
        sss[traits[locus]] += utl::sqr(effect);
    }

    for (size_t trait = 0u; trait < 3u; ++trait) {
        sss[trait] = sss[trait] > 0.0 ? sqrt(sss[trait]) : 1.0;
        assert(sss[trait] > 0.0);
    }

    for (size_t locus = 0u; locus < p.nloci; ++locus)
        effectsizes[locus] /= sss[traits[locus]];

    return effectsizes;
}


std::vector<double> GenArch::makeDominances(const Param &p) const
{

    if (p.dominancevar == 0.0) return std::vector<double>(p.nloci, 1.0);

    std::vector<double> coefficients(p.nloci);
    std::vector<double> sss = std::vector<double>(3u, 0.0);
    // square rooted sum of squares

    // Dominance coefficients are sampled from a half-normal distribution
    auto getdominance = rnd::normal(0.0, p.dominancevar);

    for (size_t locus = 0u; locus < p.nloci; ++locus) {
        double dom = getdominance(rnd::rng);
        if (dom < 0.0) dom *= -1.0;
        assert(dom >= 0.0);
        coefficients[locus] = dom;
        sss[traits[locus]] += utl::sqr(dom);
    }

    for (size_t trait = 0u; trait < 3u; ++trait) {
        sss[trait] = sss[trait] > 0.0 ? sqrt(sss[trait]) : 1.0;
        assert(sss[trait] > 0.0);
    }

    for (size_t locus = 0u; locus < p.nloci; ++locus)
        coefficients[locus] /= sss[traits[locus]];

    return coefficients;
}

MultiNet GenArch::makeNetworks(const Param &p) const
{
    MultiNet multinet;

    for (size_t trait = 0u; trait < 3u; ++trait)
        multinet.push_back(Network(trait, p, traits));

    // The indices in these network maps are indices among the loci underlying
    // a given trait,
    // not absolute loci indices across the genome

    assert(multinet.size() == 3u);

    return multinet;
}

// Functions to write the content of the genetic architecture

// Write vector as a row in text file, with end of line
void GenArch::write(const std::vector<double> &v, std::ofstream &file, const char &sep) const
{
    for (auto x : v)
        file << x << sep;
    file << '\n';
}

// Same for vector of integers
void GenArch::write(const std::vector<size_t> &v, std::ofstream &file, const char &sep) const
{
    for (auto x : v)
        file << x << sep;
    file << '\n';
}

// Same for vector of pairs (i determines first or second member)
void GenArch::write(const std::vector<Edge> &v, std::ofstream &file, const bool &i, const char &sep) const
{
    for (size_t p = 0u; p < v.size(); ++p)
        file << (i ? v[p].second : v[p].first) << sep;
    file << '\n';
}

void GenArch::save(Param &pars) const
{
    std::ofstream archfile(pars.archfile); // should be arg

    if (!archfile.is_open())
        throw std::runtime_error("Unable to open file " + pars.archfile + '\n');

    // Write parameters first
    archfile << "Parameters used to generate the architecture:\n";
    pars.write(archfile);

    archfile << "\nArchitecture:\n";

    archfile << "chromosomes\n";
    write(chromosomes, archfile);

    archfile << "traits\n";
    write(traits, archfile);

    archfile << "locations\n";
    write(locations, archfile);

    archfile << "effects\n";
    write(effects, archfile);

    archfile << "dominances\n";
    write(dominances, archfile);

    for (size_t trait = 0u; trait < 3u; ++trait) {

        archfile << "\nnetwork " << trait << ' ';
        archfile << networks[trait].edges.size() << '\n';

        archfile << "edges0\n";
        write(networks[trait].edges, archfile, false);

        archfile << "edges1\n";
        write(networks[trait].edges, archfile, true);

        archfile << "weights\n";
        write(networks[trait].weights, archfile);

    }

    archfile.close();
}


void GenArch::read(std::vector<double> &v, const size_t &n, std::ifstream &file)
{
    for (size_t i = 0u; i < n; ++i)
        file >> v[i];
}

void GenArch::read(std::vector<size_t> &v, const size_t &n, std::ifstream &file)
{
    for (size_t i = 0u; i < n; ++i)
        file >> v[i];
}

void GenArch::read(std::vector<Edge> &v, const size_t &n, const bool &id, std::ifstream &file)
{
    double x;
    for (size_t p = 0u; p < n; ++p) {
        file >> x;
        if (id)
            v[p].second = static_cast<size_t>(x);
        else
            v[p].first = static_cast<size_t>(x);
    }

}


void GenArch::load(const Param &pars)
{

    // This function will overwrite the genetic architecture
    // with that found in the arhictecture file provided
    // and update the parameters accordingly

    const std::string filename = pars.archfile;

    // Open the architecture file
    std::ifstream file(filename.c_str());
    if (!file.is_open())
        throw std::runtime_error("Unable to open file " + filename + '\n');

    // Prepare to read parameters
    std::string field;
    size_t nchrom = 1u;
    size_t nloci = 0u;

    // Read in parameters of interest first
    do {

        file >> field;

        if (field == "nchrom") file >> nchrom;
        if (field == "nvertices") {

            for (size_t trait = 0u; trait < 3u; ++trait) {
                size_t nvertices;
                file >> nvertices;
                nloci += nvertices;
            }
        }
    }
    while (field != "Architecture:");

    // Reset the architecture
    chromosomes.resize(nchrom);
    traits.resize(nloci);
    locations.resize(nloci);
    effects.resize(nloci);
    dominances.resize(nloci);

    assert(chromosomes.size() == nchrom);
    assert(traits.size() == nloci);
    assert(locations.size() == nloci);
    assert(effects.size() == nloci);
    assert(dominances.size() == nloci);

    // Prepare to read architecture
    size_t trait = 0u;
    size_t nedges = 0u;

    // Read in architecture
    while (file >> field) {

        if (field == "chromosomes") read(chromosomes, nchrom, file);
        else if (field == "traits") read(traits, nloci, file);
        else if (field == "locations") read(locations, nloci, file);
        else if (field == "effects") read(effects, nloci, file);
        else if (field == "dominances") read(dominances, nloci, file);

        else if (field == "network") {
            file >> trait;
            file >> nedges;
            networks[trait].edges.resize(nedges);
            networks[trait].weights.resize(nedges);
            assert(networks[trait].edges.size() == nedges);
            assert(networks[trait].weights.size() == nedges);
        }

        else if (field == "edges0") read(networks[trait].edges, nedges, false, file);
        else if (field == "edges1") read(networks[trait].edges, nedges, true, file);
        else if (field == "weights") read(networks[trait].weights, nedges, file);

    }

    file.close();

    // Update relevant parameters
    pars.nchrom = chromosomes.size();
    pars.nloci = locations.size();
    for (size_t i = 0u; i < 3u; ++i) {
        pars.nvertices[i] = networks[i].loci.size();
        pars.nedges[i] = networks[i].map.size();
    }
}
