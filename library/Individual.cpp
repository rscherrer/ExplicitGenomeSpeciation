#include "Individual.h"

Genome Individual::genomize(const Param &p) const
{

    // Generate a bitset of size 2N (diploid genome)
    // For each position, sample based on allele frequency

    Genome sequence(2u * p.nloci); // diploid genome

    assert(sequence.size() == 2u * p.nloci);
    assert(sequence.count() == 0u);

    auto ismutation = rnd::bernoulli(p.allfreq);

    for (size_t i = 0u; i < sequence.size(); ++i)
        if (ismutation(rnd::rng)) sequence.set(i);

    assert(sequence.size() == 2u * p.nloci);

    return sequence;
}


void Individual::recombine(Genome &zygote, const Param &p, const GenArch &arch)
 const
{

    size_t locus = 0u;
    size_t chrom = 0u;

    // Crossover points are sampled from an exponential distribution
    auto nextcrossover = rnd::exponential(p.recombination);

    // Haplotypes have equal chances to be transmitted
    auto gethaplotype = rnd::bernoulli(0.5);

    double crossover = nextcrossover(rnd::rng);
    double position = arch.locations[0u];
    double chromend = arch.chromosomes[0u];

    size_t hap = gethaplotype(rnd::rng);

    while (locus < p.nloci) {

        // What is the thing coming up next?
        vecDbl closest = { crossover, chromend, position };
        size_t next = utl::argmin(closest);

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
    assert(genome.size() == 2u * p.nloci);

}

void Individual::mutate(Genome &zygote, const Param &p) const
{
    // The number of mutations is sampled from a Poisson distribution
    auto getnmutations = rnd::poisson(p.mutation * zygote.size());

    size_t nmut = getnmutations(rnd::rng);

    // Sample mutation targets across the genome
    auto gettarget = rnd::random(0u, zygote.size() - 1u);

    while (nmut) {
        zygote.flip(gettarget(rnd::rng));
        --nmut;
    }

    // NB: Maybe mutations should be sampled without replacement.
    // But typically there should be so few that the chances of a locus
    // being hit twice are negligible.
}

Genome Individual::fecundate(const Individual &mom, const Individual &dad,
 const Param &p, const GenArch &arch) const
{
    Genome zygote(2u * p.nloci); // diploid genome
    assert(zygote.size() == 2u * p.nloci);
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
        auto getenvnoise = rnd::normal(0.0, p.scaleE[trait]);
        const double envnoise = getenvnoise(rnd::rng);
        traitvalues[trait] = genvalues[trait] + envnoise;
    }

    // Phenotype
    ecotrait = traitvalues[0u];
    matepref = traitvalues[1u];
    neutrait = traitvalues[2u];

    // Feeding rates
    const double max = p.rdynamics ? 1.0 : p.maxfeed;
    feeding[0u] = max * exp(-p.ecosel * utl::sqr(ecotrait + 1.0));
    feeding[1u] = max * exp(-p.ecosel * utl::sqr(ecotrait - 1.0));

    assert(feeding[0u] >= 0.0);
    assert(feeding[1u] >= 0.0);
    assert(feeding[0u] <= 1.0);
    assert(feeding[1u] <= 1.0);

}

bool Individual::isalive() const
{
    return alive;
}

void Individual::disperse()
{
    habitat = habitat == 0u ? 1u : 0u;
}

void Individual::feed(const vecDbl &food)
{

    fitness = feeding[0u] * food[0u] + feeding[1u] * food[1u];
    ecotype = feeding[1u] * food[1u] > feeding[0u] * food[0u];

    assert(fitness >= 0.0);
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
    if (matepref >= 0.0)
       prob = calcAssortProb(matepref, ecotrait, x, p.sexsel);
    else
       prob = calcDisassortProb(matepref, ecotrait, x, p.sexsel);
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
