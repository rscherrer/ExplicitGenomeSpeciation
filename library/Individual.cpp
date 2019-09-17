#include "Individual.h"
#include "Population.h"
#include "utils.h"
#include "Random.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>


struct Locus;

/// Function to calculate feeding rates
std::vector<double> calcFeedingRates(const double &sel, const double &trait,
 const double &maxi)
{
    const double rate1 = maxi * exp(- sel * sqr(trait + 1.0));
    const double rate2 = maxi * exp(- sel * sqr(trait - 1.0));
    assert(rate1 >= 0.0);
    assert(rate2 >= 0.0);
    return { rate1, rate2 };
}


/// Constructor with randomly generated genome
Individual::Individual(const Genome &genome, const MultiNet &networks,
 const double &snpfreq, const double &scaleE) :
    sequence(makeSequence(genome.nloci, snpfreq)),
    genexp(zeros(genome.nloci)),
    locusGeneticValues(zeros(genome.nloci)),
    isFemale(determineSex(genome.femgamy)),
    geneticValues(develop(genome, networks)),
    ecoTrait(geneticValues[0u] + rnd::normal(0.0, scaleE)),
    matePref(geneticValues[1u] + rnd::normal(0.0, scaleE)),
    neutral(geneticValues[2u] + rnd::normal(0.0, scaleE)),
    traitvalues({ ecoTrait, matePref, neutral }),
    fitness(1.0),
    feedingRates(calcFeedingRates(1.0, ecoTrait)),
    ecotype(0u)
{
    assert(sequence.size() == 2u);
    for (size_t strain = 0u; strain < 2u; ++strain)
        assert(sequence[strain].size() == genome.nloci);
    assert(genexp.size() == genome.nloci);
    assert(traitvalues.size() == 3u);
    assert(fitness > 0.0);
    for (size_t res = 0u; res < 2u; ++res)
        assert(feedingRates[res] > 0.0);
}


/// Constructor that inherits a parental genome
Individual::Individual(const Genome &genome,
 const MultiNet &networks, const Haplotype &egg, const Haplotype &sperm,
  const double &scaleE) :
    sequence(fecundate(egg, sperm)),
    genexp(zeros(genome.nloci)),
    locusGeneticValues(zeros(genome.nloci)),
    isFemale(determineSex(genome.femgamy)),
    geneticValues(develop(genome, networks)),
    ecoTrait(geneticValues[0u] + rnd::normal(0.0, scaleE)),
    matePref(geneticValues[1u] + rnd::normal(0.0, scaleE)),
    neutral(geneticValues[2u] + rnd::normal(0.0, scaleE)),
    traitvalues({ ecoTrait, matePref, neutral }),
    fitness(1.0),
    feedingRates(calcFeedingRates(1.0, ecoTrait)),
    ecotype(0u)
{

    assert(sequence.size() == 2u);
    for (size_t strain = 0u; strain < 2u; ++strain)
        assert(sequence[strain].size() == genome.nloci);
    assert(genexp.size() == genome.nloci);
    assert(traitvalues.size() == 3u);
    assert(fitness > 0.0);
    for (size_t res = 0u; res < 2u; ++res)
        assert(feedingRates[res] > 0.0);

}


bool Individual::determineSex(const bool &femheterogamy)
{
    const bool hetero = sequence[0u][0u] != sequence[1u][0u];
    const bool isZW = hetero && femheterogamy;
    const bool isXX = !hetero && !femheterogamy;
    const bool isfemale = isZW || isXX;
    return isfemale;
}


/// Generate a diploid allele sequence
Diplotype Individual::makeSequence(const size_t &nloci, const double &prob)
{

    Diplotype sequences;

    for (size_t strain = 0u; strain < 2u; ++strain) {

        // Generate a random genetic sequence of alleles
        Haplotype haplotype;
        for (size_t locus = 0u; locus < nloci; ++locus)
            haplotype.push_back(rnd::bernoulli(prob));
        sequences.push_back(haplotype);

    }

    assert(sequences[0u].size() == sequences[1u].size());

    return sequences;
}


/// Fecundation
Diplotype Individual::fecundate(const Haplotype &egg, const Haplotype &sperm)
{
    Diplotype zygote(2u);
    zygote[0u] = egg;
    zygote[1u] = sperm;
    assert(zygote[0u].size() == zygote[1u].size());
    return zygote;
}


/// Development
std::vector<double> Individual::develop(const Genome &genome,
 const MultiNet &networks)
{

    // Development reads the genome and computes trait values
    // Loop throughout the genome
    // Each gene contributes to the trait value
    // But each gene has a certain allele in a certain individual
    // And each gene is diploid
    // And there is dominance
    // And there are multiple traits
    // And there is epistasis...

    vecDbl phenotypes {0.0, 0.0, 0.0};

    for (size_t locus = 0u; locus < genome.nloci; ++locus) {

        // Determine genotype
        size_t genotype = 0u;
        for (size_t hap = 0u; hap < 2u; ++hap)
            genotype += sequence[hap][locus];

        // Determine gene expression
        double expression;
        const double dominance = genome.dominances[locus];
        switch(genotype) {
            case 1u : expression = dominance; break; // Aa
            case 2u : expression = 1.0; break; // AA
            default : expression = -1.0; break; // aa
        }

        assert(expression >= -1.0);
        assert(expression <= 1.0);

        genexp[locus] = expression; // record gene expression

        // Determine the encoded trait
        const size_t trait = genome.traits[locus];

        // Contribute to trait
        locusGeneticValues[locus] = genome.effects[locus] * expression;
        phenotypes[trait] += locusGeneticValues[locus];

    }

    // For epistasis,
    // For each trait,
    // loop through the edges of the network
    // record the expression level of both partners and multiply them
    // multiply the result by the weight of the interaction
    // add the result to the phenotype

    for (size_t trait = 0u; trait < 3u; ++trait) {

        for (size_t e = 0u; e < networks[trait].nedges; ++e) {

            assert(networks[trait].edges.size() > 0u);

            // Level of expression of an interaction
            const Edge edge = networks[trait].edges[e];
            const double intexp = genexp[edge.first] * genexp[edge.second];

            assert(intexp >= -1.0);
            assert(intexp <= 1.0);

            const double interaction = intexp * networks[trait].weights[e];
            locusGeneticValues[edge.first] += 0.5 * interaction;
            locusGeneticValues[edge.second] += 0.5 * interaction;
            phenotypes[trait] += interaction;

        }
    }

    // Normalize!

    return phenotypes;
}


/// Feed and get a fitness value
void Individual::feed(const vecDbl &food)
{
    fitness = feedingRates[0u] * food[0u] + feedingRates[1u] * food[1u];
    assert(fitness >= 0.0);
}


/// Function to calculate mating probability under homogamy
double calcAssortProb(const double &y, const double &xi,
 const double &xj, const double &alpha)
{
    const double d = xi - xj;
    return exp(- 0.5 * alpha * sqr(y * d));
}


/// Function to calculate mating probability under heterogamy
double calcDisassortProb(const double &y, const double &xi,
 const double &xj, const double &alpha)
{
    return 1.0 - sqr(sqr(y)) * calcAssortProb(y, xi, xj, alpha);
}


/// Function to evaluate a potential mate
bool Individual::acceptMate(const double &xj, const double &strength) const
{

    const double tiny = 0.00000001;

    // Calculate the probability of mating
    double mateProb = matePref >= 0.0 ?
     calcAssortProb(matePref, ecoTrait, xj, strength) :
      calcDisassortProb(matePref, ecoTrait, xj, strength);

    if (mateProb < tiny) mateProb = 0.0;
    if (mateProb > 1.0 - tiny) mateProb = 1.0;

    assert(mateProb >= 0.0);
    assert(mateProb <= 1.0);

    // Sample mating event
    return rnd::bernoulli(mateProb);

}


/// Meiosis to produce a gamete
Haplotype Individual::recombine(const vecDbl &locations,
 const vecDbl &chromosomes, const double &rate)
{
    Haplotype gamete;

    // Choose a random haplotype
    // Loop through loci along this haplotype
    // Add this locus to the inherited gamete
    // Yes but there are crossing overs
    // Upon a crossover the strain changes
    // The crossover point is a location along the genome, not a specific locus
    // There can be several crossover points
    // The rate of recombination can be provided
    // A rate of 3 means three crossovers are expected across the genome
    // which means that the genome size is equivalent to 300cM
    // 1cM = 1% change recombination
    // But wait, there is free recombination between the chromosomes!

    const size_t nloci = sequence[0u].size();

    size_t locus = 0u;
    size_t chrom = 0u;

    double crossover = rnd::exponential(rate);
    double position = locations[0u];
    double chromend = chromosomes[0u];

    size_t hap = 0u;

    while (locus < nloci) {

        // What is the next thing coming up next?
        vecDbl closest = { crossover, chromend, position };
        size_t next = argmin(closest);

        switch (next) {

        // Crossover point
        case 0u:
            hap = hap ? 0u : 1u;
            crossover += rnd::exponential(rate);
            break;

        // Free recombination point
        case 1u:
            hap = rnd::random(2u);
            ++chrom;
            chromend = chromosomes[chrom];
            break;

        // Gene
        default:
            gamete.push_back(sequence[hap][locus]);
            ++locus;
            //assert(locus < locations.size());
            position = locations[locus];
            break;

        }
    }

    assert(locus == nloci);
    assert(chrom == chromosomes.size() - 1u);
    assert(gamete.size() == nloci);

    return gamete;
}


/// Mutation
void Individual::mutate(Haplotype &gamete, const double &rate)
{

    // Sample a number of mutations from a poisson
    // Sample the mutated targets
    // Flip the alleles

    size_t nmut = rnd::poisson(rate * gamete.size());

    while (nmut) {

        size_t target = rnd::random(gamete.size());
        gamete[target] = gamete[target] ? 0u : 1u;

        --nmut;
    }

}


void Individual::setEcotype(const double &mean)
{

    ecotype = ecoTrait > mean;

}

size_t Individual::getZygosity(const size_t &locus)
{
    return sequence[0u][locus] + sequence[1u][locus];
}

double Individual::getLocusGenValue(const size_t &locus)
{
    return locusGeneticValues[locus];
}

