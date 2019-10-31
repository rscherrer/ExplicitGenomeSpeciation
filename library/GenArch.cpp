#include "GenArch.h"

bool GenArch::resetseed(const size_t &seed) const
{
    rnd::rng.seed(seed);
    return true;
}

vecDbl GenArch::makeChromosomes(const Param &p) const
{

    vecDbl chromends(p.nchrom);

    // Chromosomes all have the same size
    for (size_t chrom = 0u; chrom < p.nchrom; ++chrom)
        chromends[chrom] = (chrom + 1.0) / p.nchrom;

    return chromends;

}

vecUns GenArch::makeEncodedTraits(const Param &p) const
{

    vecUns encoded(p.nloci);

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

    vecUns nvertices = utl::uzeros(3u);
    for (size_t locus = 0u; locus < p.nloci; ++locus)
        ++nvertices[encoded[locus]];

    for (size_t trait = 0u; trait < 3u; ++trait)
        assert(nvertices[trait] == p.nvertices[trait]);

    return encoded;

}



vecDbl GenArch::makeLocations(const Param &p) const
{
    vecDbl positions(p.nloci);

    // Locations are sampled from a uniform distribution between 0 and 1
    auto getlocation = rnd::uniform(0.0, 1.0);

    for (size_t locus = 0u; locus < p.nloci; ++locus)
        positions[locus] = getlocation(rnd::rng);

    std::sort(positions.begin(), positions.end());

    for (size_t locus = 1u; locus < p.nloci; ++locus)
        assert(positions[locus] > positions[locus - 1u]);

    return positions;
}


vecDbl GenArch::makeEffects(const Param &p) const
{

    if (p.effectshape == 0.0 || p.effectscale == 0.0)
        return utl::zeros(p.nloci);

    vecDbl effectsizes(p.nloci);
    vecDbl sss = utl::zeros(3u); // square rooted sum of squares

    // Effect sizes are sampled from a two-sided Gamma distribution
    auto getffect = rnd::gamma(p.effectshape, p.effectscale);
    auto isflipped = rnd::bernoulli(0.5);

    for (size_t locus = 0u; locus < p.nloci; ++locus) {

        double effect = getffect(rnd::rng);
        if (isflipped(rnd::rng)) effect *= -1.0;
        effectsizes[locus] = effect;
        sss[traits[locus]] += utl::sqr(effect);
    }

    for (size_t trait = 0u; trait < 3u; ++trait) {
        sss[trait] = sqrt(sss[trait]);
        assert(sss[trait] > 0.0);
    }

    for (size_t locus = 0u; locus < p.nloci; ++locus)
        effectsizes[locus] /= sss[traits[locus]];

    return effectsizes;
}


vecDbl GenArch::makeDominances(const Param &p) const
{

    if (p.dominancevar == 0.0) return utl::zeros(p.nloci);

    vecDbl coefficients(p.nloci);
    vecDbl sss = utl::zeros(3u); // square rooted sum of squares

    // Dominance coefficients are sampled from a half-normal distribution
    auto getdominance = rnd::normal(0.0, p.dominancevar);

    for (size_t locus = 0u; locus < p.nloci; ++locus) {
        double dom = getdominance(rnd::rng);
        if (dom < 0.0) dom *= -1.0;
        assert(dom >= 0.0);
        coefficients[locus] = dom;
        sss[traits[locus]] += utl::sqr(dom);
    }

    for (size_t trait = 0u; trait < 3u; ++trait)
        sss[trait] = sqrt(sss[trait]);

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

vecStrings GenArch::whattosave() const
{
    return {

        "architecture_chromosomes",
        "architecture_traits",
        "architecture_locations",
        "architecture_effects",
        "architecture_dominances",
        "architecture_edges1",
        "architecture_edges2",
        "architecture_weights"

    };
}
