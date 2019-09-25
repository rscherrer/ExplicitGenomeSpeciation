#include "Individual.h"
#include "Deme.h"
#include "Utilities.h"
#include "Random.h"
#include <cmath>
#include <cassert>
#include <iostream>
#include <algorithm>


size_t Individual::getAlleleSum(const size_t &hap)
{
    return genome[hap].count();
}

bool Individual::checkIndividual(const size_t &nLoci)
{

    assert(genome.size() == 2u);
    for (size_t hap = 0u; hap < 2u; ++hap) {
        assert(genome[hap].size() == nLoci);
    }
    assert(transcriptome.size() == nLoci);
    assert(traitvalues.size() == 3u);
    assert(fitness > 0.0);
    for (size_t res = 0u; res < 2u; ++res)
        assert(feedingRates[res] > 0.0);

    // For the sake of not breaking in release mode
    const double out = nLoci;
    return nLoci == out;
}

vecDbl Individual::calcFeedingRates(const double &sel, const double &trait,
 const double &maxi)
{
    const double rate1 = maxi * exp(- sel * utl::sqr(trait + 1.0));
    const double rate2 = maxi * exp(- sel * utl::sqr(trait - 1.0));
    assert(rate1 >= 0.0);
    assert(rate2 >= 0.0);
    return { rate1, rate2 };
}

bool Individual::determineSex(const bool &femheterogamy)
{

    // The first locus determines the sex
    // If the individual is heterozygous
    // and if females are the heterozygous sex, then ZW
    // otherwise, if males are the heterozygous sex, then XY

    const bool hetero = genome[0u][0u] != genome[1u][0u];
    const bool isZW = hetero && femheterogamy;
    const bool isXX = !hetero && !femheterogamy;
    const bool isfemale = isZW || isXX;
    return isfemale;

}

Genome Individual::generateGenome(const GenArch &arch, double prob)
{

    Genome sequences;

    // SNP frequency is read in architecture if not provided explicitly
    if (prob < 0.0) prob = arch.snpFreq;

    assert(prob <= 1.0);
    assert(prob >= 0.0);

    // For each haplotype
    for (size_t hap = 0u; hap < 2u; ++hap) {

        // Generate a genetic sequence of N zero-alleles
        Haplotype haplotype(arch.nLoci);

        // Throw mutations here and there
        if (prob < 0.5) {

            // Binomial distribution
            size_t nmut = rnd::binomial(arch.nLoci, prob);
            vecDbl probs = utl::ones(arch.nLoci);
            while (nmut) {
                size_t mutant = rnd::sample(probs);
                haplotype.set(mutant);
                probs[mutant] = 0.0; // without replacement
                --nmut;
            }

        } else {

            // Bernoulli events
            for (size_t locus = 0u; locus < arch.nLoci; ++locus)
                if (rnd::bernoulli(prob))
                    haplotype.set(locus);

        }

        sequences.push_back(haplotype);

    }

    assert(sequences[0u].size() == sequences[1u].size());

    return sequences;
}

Genome Individual::fecundate(const Haplotype &egg, const Haplotype &sperm)
{
    Genome zygote(2u); // diploid genome
    zygote[0u] = egg;
    zygote[1u] = sperm;
    assert(zygote[0u].size() == zygote[1u].size());
    return zygote;
}

void Individual::develop(const GenArch &arch)
{

    // Development reads the genome and computes trait values
    // Loop throughout the genome
    // Each gene contributes to the trait value
    // But each gene has a certain allele in a certain individual
    // And each gene is diploid
    // And there is dominance
    // And there are multiple traits
    // And there is epistasis...

    for (size_t locus = 0u; locus < arch.nLoci; ++locus) {

        // Determine the encoded trait
        const size_t trait = arch.traits[locus];

        // Determine genotype
        size_t genotype = 0u;
        for (size_t hap = 0u; hap < 2u; ++hap)
            genotype += genome[hap][locus];

        // Determine gene expression
        double expression;
        const double dominance = arch.dominances[locus];
        switch(genotype) {
            case 1u : expression = arch.scaleD[trait] * dominance; break;
            case 2u : expression = 1.0; break;
            default : expression = -1.0; break;
        }

        assert(expression >= -1.0);
        assert(expression <= 1.0);

        transcriptome[locus] = expression; // record gene expression

        // Contribute to trait
        double locuseffect = arch.effects[locus] * expression;
        locuseffect *= arch.scaleA[trait];
        locivalues[locus] = locuseffect;
        genvalues[trait] += locivalues[locus];

    }

    // For epistasis,
    // For each trait,
    // loop through the edges of the network
    // record the expression level of both partners and multiply them
    // multiply the result by the weight of the interaction
    // add the result to the phenotype

    for (size_t trait = 0u; trait < 3u; ++trait) {

        for (size_t e = 0u; e < arch.networks[trait].nedges; ++e) {

            assert(arch.networks[trait].edges.size() > 0u);

            // Level of expression of an interaction
            const Edge edge = arch.networks[trait].edges[e];
            double intexp = transcriptome[edge.first];
            intexp *= transcriptome[edge.second];

            assert(intexp >= -1.0);
            assert(intexp <= 1.0);

            double interaction = intexp * arch.networks[trait].weights[e];
            interaction *= arch.scaleI[trait];
            locivalues[edge.first] += 0.5 * interaction;
            locivalues[edge.second] += 0.5 * interaction;
            genvalues[trait] += interaction;

        }
    }

    for (size_t trait = 0u; trait < 3u; ++trait) {
        const double envnoise = rnd::normal(0.0, arch.scaleE[trait]);
        traitvalues[trait] = genvalues[trait] + envnoise;
    }

    ecotrait = traitvalues[0u];
    matepref = traitvalues[1u];
    neutrait = traitvalues[2u];
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
    return exp(- 0.5 * alpha * utl::sqr(y * d));
}


/// Function to calculate mating probability under heterogamy
double calcDisassortProb(const double &y, const double &xi,
 const double &xj, const double &alpha)
{
    return 1.0 - utl::sqr(utl::sqr(y)) * calcAssortProb(y, xi, xj, alpha);
}


/// Function to evaluate a potential mate
bool Individual::acceptMate(const double &xj, const double &sexsel) const
{

    const double tiny = 1E-15;

    // Calculate the probability of mating
    double mateProb = matepref >= 0.0 ?
     calcAssortProb(matepref, ecotrait, xj, sexsel) :
      calcDisassortProb(matepref, ecotrait, xj, sexsel);

    if (mateProb < tiny) mateProb = 0.0;
    if (mateProb > 1.0 - tiny) mateProb = 1.0;

    assert(mateProb >= 0.0);
    assert(mateProb <= 1.0);

    // Sample mating event
    return rnd::bernoulli(mateProb);

}


Haplotype Individual::recombine(const GenArch &arch)
{

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
    // 1cM = 1% chance recombination
    // But wait, there is free recombination between the chromosomes!

    size_t locus = 0u;
    size_t chrom = 0u;

    double crossover = rnd::exponential(arch.recombinationRate);
    double position = arch.locations[0u];
    double chromend = arch.chromosomes[0u];

    Haplotype gamete = genome[0u];

    size_t hap = rnd::bernoulli(0.5);

    while (locus < arch.nLoci) {

        // What is the thing coming up next?
        vecDbl closest = { crossover, chromend, position };
        size_t next = utl::argmin(closest);

        switch (next) {

        // Crossover point
        case 0u:
            hap = hap ? 0u : 1u;
            crossover += rnd::exponential(arch.recombinationRate);
            break;

        // Free recombination point
        case 1u:
            hap = rnd::random(2u);
            ++chrom;
            chromend = arch.chromosomes[chrom];
            break;

        // Gene
        default:
            if (hap == 1u)
                if (genome[0u].test(locus) != genome[1u].test(locus))
                    gamete.flip(locus);
            ++locus;
            position = arch.locations[locus];
            break;

        }
    }

    assert(locus == arch.nLoci);
    assert(chrom == arch.nChromosomes - 1u);
    assert(gamete.size() == arch.nLoci);

    return gamete;
}


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

    ecotype = ecotrait > mean;

}

size_t Individual::getZygosity(const size_t &locus)
{
    return genome[0u][locus] + genome[1u][locus];
}

double Individual::getLocusValue(const size_t &locus)
{
    return locivalues[locus];
}

void Individual::setEcoTrait(const double &value, const double &sel,
 const double &max)
{
    ecotrait = value;
    traitvalues[0u] = value;
    feedingRates = calcFeedingRates(sel, value, max);
}
void Individual::setMatePref(const double &value)
{
    matepref = value;
    traitvalues[1u] = value;
}

void Individual::setGender(const bool &sex)
{
    isFemale = sex;
}

