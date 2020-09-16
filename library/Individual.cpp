#include "Individual.h"

// Constructors
//--------------

// To generate an initial population
Individual::Individual(const Param &pars, const GenArch &arch) :
    genome(genomize(pars)),
    transcriptome(std::vector<double>(pars.nloci, 0.0)),
    locivalues(std::vector<double>(pars.nloci, 0.0)),
    genvalues(std::vector<double>(3u, 0.0)),
    traitvalues(std::vector<double>(3u, 0.0)),
    midparents(std::vector<double>(3u, 0.0)),
    fitness(1.0),
    feeding(std::vector<double>(3u, 0.0)),
    ecotype(0u),
    habitat(0u),
    gender(determinesex()),
    alive(true),
    adult(true)
{
    develop(pars, arch);

    assert(transcriptome.size() == pars.nloci);
    assert(traitvalues.size() == 3u);
    assert(fitness >= 0.0);
    assert(feeding[0u] >= 0.0);
    assert(feeding[1u] >= 0.0);
    assert(feeding[0u] <= 1.0);
    assert(feeding[1u] <= 1.0);
}

// Newborn
Individual::Individual(const Param &pars, const GenArch &arch,
 const Individual &mom, const Individual &dad) :
    genome(fecundate(mom, dad, pars, arch)),
    transcriptome(std::vector<double>(pars.nloci, 0.0)),
    locivalues(std::vector<double>(pars.nloci, 0.0)),
    genvalues(std::vector<double>(3u, 0.0)),
    traitvalues(std::vector<double>(3u, 0.0)),
    midparents(calcmidparent(mom, dad)),
    fitness(1.0),
    feeding(std::vector<double>(2u, 0.0)),
    ecotype(0u),
    habitat(mom.getHabitat()),
    gender(determinesex()),
    alive(true),
    adult(false)
{
    develop(pars, arch);

    assert(transcriptome.size() == pars.nloci);
    assert(traitvalues.size() == 3u);
    assert(fitness >= 0.0);
    assert(feeding[0u] >= 0.0);
    assert(feeding[1u] >= 0.0);
    assert(feeding[0u] <= 1.0);
    assert(feeding[1u] <= 1.0);
}

// Member functions
//-----------------

Genome Individual::genomize(const Param &p) const
{

    // Generate a bitset of size 2N (diploid genome)
    // For each position, sample based on allele frequency

    Genome sequence; // diploid genome

    assert(sequence.count() == 0u);

    if (p.allfreq == 0.0) return sequence;

    if (p.allfreq < 0.1) {

        mutate(sequence, p);

    }
    else {

        // Use a bernoulli if common
        auto ismutation = rnd::bernoulli(p.allfreq);
        for (size_t i = 0u; i < p.nloci * 2u; ++i)
            if (ismutation(rnd::rng)) sequence.set(i);

    }

    return sequence;
}


void Individual::recombine(Genome &zygote, const Param &p, const GenArch &arch)
 const
{

    size_t locus = 0u;
    size_t chrom = 0u;        

    // Haplotypes have equal chances to be transmitted
    auto gethaplotype = rnd::bernoulli(0.5);

    // Crossovers are sampled from an exponential distribution
    const double recombrate = p.recombination > 0.0 ? p.recombination : 100.0;
    auto nextcrossover = rnd::exponential(recombrate);

    double crossover = 1.1; // beyond the end of the genome
    if (p.recombination > 0.0) crossover = nextcrossover(rnd::rng);

    double position = arch.locations[0u];
    double chromend = arch.chromosomes[0u];

    size_t hap = gethaplotype(rnd::rng);

    while (locus < p.nloci) {

        // What is the thing coming up next?
        size_t next = static_cast<size_t>(crossover);
        if (crossover < chromend && crossover < position)
            next = 0u;
        else if (chromend < crossover && chromend < position)
            next = 1u;
        else
            next = 2u;

        switch (next) {

        // Upon crossover point, switch haplotype
        case 0u:
            hap = hap ? 0u : 1u;
            crossover += nextcrossover(rnd::rng);
            break;

        // Upon free recombination point, switch to random chromosome
        case 1u:
            hap = gethaplotype(rnd::rng);
            ++chrom;
            if (chrom < p.nchrom) chromend = arch.chromosomes[chrom];
            assert(chrom < p.nchrom);
            break;

        // Upon gene, transmit haplotype to the zygote
        default:
            assert(locus + hap * p.nloci < 2u * p.nloci);
            assert(locus + gender * p.nloci < 2u * p.nloci);
            if (genome.test(locus + hap * p.nloci))
                zygote.set(locus + gender * p.nloci);
            ++locus;
            if (locus < p.nloci) position = arch.locations[locus];
            break;
        }
    }

    assert(locus == p.nloci);
    assert(chrom == p.nchrom - 1u);

}

void Individual::mutate(Genome &seq, const Param &p) const
{

    if (p.mutation == 0.0) return;
    if (p.mutation == 1.0)
        for (size_t i = 0u; i < 2u * p.nloci; ++i)
            seq.set(i);

    // Mutations are sampled from a geometric distribution
    assert(p.mutation > 0.0);
    auto getnextmutant = rnd::iotagap(p.mutation);
    getnextmutant.reset(0u);
    for (;;) {
        const size_t mut = getnextmutant(rnd::rng);
        if (mut >= 2.0 * p.nloci) break;
        seq.flip(mut);
    }
}

Genome Individual::fecundate(const Individual &mom, const Individual &dad,
 const Param &p, const GenArch &arch) const
{
    Genome zygote; // diploid genome
    assert(zygote.count() == 0u);

    mom.recombine(zygote, p, arch);
    dad.recombine(zygote, p, arch);
    mutate(zygote, p);

    return zygote;
}

void Individual::develop(const Param &p, const GenArch &arch)
{

    // Accumulate independent effects for each locus
    for (size_t locus = 0u; locus < p.nloci; ++locus) {

        // Determine the encoded trait
        const size_t trait = arch.traits[locus];

        // Determine the genotype
        size_t genotype = genome.test(locus) + genome.test(locus + p.nloci);

        // Determine gene expression
        double expression;
        const double dominance = arch.dominances[locus];
        switch (genotype) {
            case 1u : expression = p.scaleD[trait] * dominance; break; // Aa
            case 2u : expression = 1.0; break; // AA
            default : expression = -1.0; break; // aa
        }

        assert(expression >= -1.0);

        transcriptome[locus] = expression; // record gene expression

        // Contribute to trait
        double locuseffect = arch.effects[locus] * expression;
        locuseffect *= p.scaleA[trait];
        locivalues[locus] = locuseffect;
        genvalues[trait] += locivalues[locus];

    }

    // For each network...
    for (size_t trait = 0u; trait < 3u; ++trait) {

        // Accumulate interaction effects for each edge
        for (size_t e = 0u; e < p.nedges[trait]; ++e) {

            assert(arch.networks[trait].edges.size() > 0u);

            // Level of expression of an interaction
            const Edge edge = arch.networks[trait].edges[e];
            double intexp = transcriptome[edge.first];
            intexp *= transcriptome[edge.second];

            assert(intexp >= -1.0);
            assert(intexp <= 1.0);

            double interaction = intexp * arch.networks[trait].weights[e];
            interaction *= p.scaleI[trait];
            locivalues[edge.first] += 0.5 * interaction;
            locivalues[edge.second] += 0.5 * interaction;
            genvalues[trait] += interaction;

        }
    }

    // Add normally distributed environmental effect for each trait
    for (size_t trait = 0u; trait < 3u; ++trait) {
        double envnoise = 0.0;
        if (p.scaleE[trait] > 0.0) {
            auto getenvnoise = rnd::normal(0.0, p.scaleE[trait]);
            envnoise = getenvnoise(rnd::rng);
        }
        traitvalues[trait] = genvalues[trait] + envnoise;
    }

    // Feeding rates
    const double selection = p.ecosel / utl::sqr(p.ecoscale);
    feeding[0u] = exp(-selection * utl::sqr(traitvalues[0u] + p.ecoscale));
    feeding[1u] = exp(-selection * utl::sqr(traitvalues[0u] - p.ecoscale));

    assert(feeding[0u] >= 0.0);
    assert(feeding[1u] >= 0.0);
    assert(feeding[0u] <= 1.0);
    assert(feeding[1u] <= 1.0);

}

std::vector<double> Individual::calcmidparent(const Individual &mom,
 const Individual &dad) const
{
    std::vector<double> midtraits = std::vector<double>(3u, 0.0);
    for (size_t trait = 0u; trait < 3u; ++trait) {
        midtraits[trait] = mom.getTraitValue(trait) + dad.getTraitValue(trait);
        midtraits[trait] /= 2.0;
    }
    return midtraits;

}

bool Individual::isalive() const
{
    return alive;
}

void Individual::disperse()
{
    habitat = habitat == 0u ? 1u : 0u;
}

void Individual::feed(const std::vector<double> &food)
{
    fitness = feeding[0u] * food[0u] + feeding[1u] * food[1u];
    assert(fitness >= 0.0);
}

void Individual::classify(const double &meanx)
{
    ecotype = traitvalues[0u] > meanx;
    assert(ecotype == 0u || ecotype == 1u);
}

double calcAssortProb(const double &y, const double &xi,
 const double &xj, const double &sexsel)
{
    const double d = xi - xj;
    const double prob = exp(- 0.5 * sexsel * utl::sqr(y * d));
    assert(prob >= 0.0);
    assert(prob <= 1.0);
    return prob;
}

double calcDisassortProb(const double &y, const double &xi,
 const double &xj, const double &sexsel)
{
    double prob = 1.0;
    prob -= utl::sqr(utl::sqr(y)) * calcAssortProb(y, xi, xj, sexsel);

    if (prob <= 0.0) prob = 0.0; // can happen if y goes below -1

    assert(prob >= 0.0);
    assert(prob <= 1.0);
    return prob;
}

double Individual::mate(const double &x, const Param &p) const
{
    double prob;
    if (traitvalues[1u] >= 0.0)
       prob = calcAssortProb(traitvalues[1u], traitvalues[0u], x, p.sexsel);
    else
       prob = calcDisassortProb(traitvalues[1u], traitvalues[0u], x, p.sexsel);
    assert(prob >= 0.0);
    assert(prob <= 1.0);
    return prob;
}

void Individual::survive(const bool &x)
{
    if (!adult) {
        adult = true;
        return;
    }
    alive = x;
}

bool Individual::determinesex() const
{
    auto getsex = rnd::bernoulli(0.5);
    return getsex(rnd::rng);
}

// Various getters
//------------------

// Getters called from outside
bool Individual::getGender() const
{
    return gender;
}
size_t Individual::getEcotype() const
{
    return ecotype;
}
size_t Individual::getHabitat() const
{
    return habitat;
}
double Individual::getFitness() const
{
    return fitness;
}
double Individual::getTraitValue(const size_t &trait) const
{
    return traitvalues[trait];
}
double Individual::getMidparent(const size_t &trait) const
{
    return midparents[trait];
}
double Individual::getGenValue(const size_t &trait) const
{
    return genvalues[trait];
}
double Individual::getFeeding(const size_t &r) const
{
    return feeding[r];
}
double Individual::getLocusValue(const size_t &locus) const
{
    return locivalues[locus];
}
size_t Individual::getZygosity(const size_t &locus, const size_t &nloci) const
{
    const size_t zyg = genome.test(locus) + genome.test(locus + nloci);
    assert(zyg == 0u || zyg == 1u || zyg == 2u);
    return zyg;
}
size_t Individual::getAlleleSum() const
{
    return genome.count();
}
double Individual::getExpression() const
{
    double sum = 0.0;
    for (size_t locus = 0u; locus < transcriptome.size(); ++locus) {
        sum += transcriptome[locus];
    }
    return sum;
}

// Return an individual genome in bitset format
Genome Individual::getFullGenome() const {

    return genome;

    // Caution: the length of the genome is 10,000 because the bitset
    // is of a constant, globally defined maximum size. Only the 2 * nloci
    // first values are relevant
}

// Get the ith 64bit-chunk of the genome
std::bitset<64u> Individual::getGenomeChunk(const size_t &i) const {

    std::bitset<64u> chunk;
    for (size_t l = 0u, k = l + i * 64u; l < 64u && k < 10000u; ++l, ++k)
        if (genome.test(k)) chunk.set(l);

    return chunk;

}


// Force resetters
//----------------

// Change the trait value of an individual
void Individual::resetTrait(const size_t &trait, const double &newvalue,
 const Param &p)
{
    traitvalues[trait] = newvalue;
    if (trait == 0u) {
        feeding[0u] = exp(-p.ecosel * utl::sqr(traitvalues[trait] + 1.0));
        feeding[1u] = exp(-p.ecosel * utl::sqr(traitvalues[trait] - 1.0));
        assert(feeding[0u] >= 0.0);
        assert(feeding[1u] >= 0.0);
        assert(feeding[0u] <= 1.0);
        assert(feeding[1u] <= 1.0);
    }
}

// Change the ecotype of an individual
void Individual::resetEcotype(const size_t &e)
{
    ecotype = e;
}

// Change the gender of an individual
void Individual::resetGender(const bool &sex)
{
    gender = sex;
}
