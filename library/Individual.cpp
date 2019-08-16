#include "Individual.h"
#include "Population.h"
#include "utils.h"
#include "Random.h"
#include <cmath>
#include <cassert>
#include <iostream>

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


/// Constructor
Individual::Individual(const Genome &genome,
 const MultiNet &networks) :
    sequence(makeSequence(genome.effects.size())),
    genexp(zeros(genome.nloci)),
    isFemale(rnd::bernoulli(0.5)),
    traits(develop(genome, networks)),
    ecoTrait(traits[0u]),
    matePref(traits[1u]),
    neutral(traits[2u]),
    fitness(1.0),
    feedingRates(calcFeedingRates(1.0, ecoTrait))
{

    assert(sequence.size() == 2u);
    for (size_t strain = 0u; strain < 2u; ++strain)
        assert(sequence[strain].size() == genome.nloci);
    assert(genexp.size() == genome.nloci);
    assert(traits.size() == 3u);
    assert(fitness > 0.0);
    for (size_t res = 0u; res < 2u; ++res)
        assert(feedingRates[res] > 0.0);
}


/// Generate a diploid allele sequence
std::vector<std::vector<bool> > Individual::makeSequence(const size_t &nloci)
{

    std::vector<std::vector<bool> > sequences;

    for (size_t strain = 0u; strain < 2u; ++strain) {

        // Generate a random genetic sequence of alleles
        std::vector<bool> haplotype;
        for (size_t locus = 0u; locus < nloci; ++locus)
            haplotype.push_back(rnd::bernoulli(0.5));
        sequences.push_back(haplotype);

    }

    assert(sequences[0u].size() == sequences[1u].size());

    return sequences;
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

    std::vector<double> phenotypes {0.0, 0.0, 0.0};

    for (size_t locus = 0u; locus < sequence[0u].size(); ++locus) {

        // Determine genotype
        size_t genotype = 0u;
        for (size_t strain = 0u; strain < 2u; ++strain)
            if (sequence[strain][locus])
                ++genotype;

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
        phenotypes[trait] += genome.effects[locus] * expression;

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

            phenotypes[trait] += intexp * networks[trait].weights[e];

        }
    }

    return phenotypes;
}


/// Feed and get a fitness value
void Individual::feed(const std::vector<double> &food)
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

/*
// Accessory functions
double calcAssortProb(const double &matePreference, const double &matingTrait,
 const double &ecoTraitDistance)
{
    return exp(- matePreference * matingTrait * ecoTraitDistance * 0.5);
}

double calcDisassortProb(const double &matePreference,
 const double &matingTrait, const double &ecoTraitDistance)
{
    return 1.0 - sqr(sqr(matingTrait)) * exp(- matePreference * matingTrait *
     ecoTraitDistance * 0.5);
}
*/


/*
// Constructors
Individual::Individual(const ParameterSet& parameters,
 const GeneticArchitecture &geneticArchitecture) :
isHeterogamous(rnd::bernoulli(0.5)), ecotype(0u), habitat(0u)
{

    setGenomeSequence(parameters.nBits, parameters.freqSNP);
    mutate(parameters);
    develop(parameters, geneticArchitecture);

}

Individual::Individual(const std::vector<bool>& sequence,
 const ParameterSet& parameters,
  const GeneticArchitecture &geneticArchitecture) : genomeSequence(sequence),
   isHeterogamous(rnd::bernoulli(0.5)), habitat(0u), ecotype(0u)
{
    mutate(parameters);
    develop(parameters, geneticArchitecture);
}

Individual::Individual(Individual const * const mother,
 Individual const * const father, const ParameterSet& parameters,
  const GeneticArchitecture &geneticArchitecture) :
        isHeterogamous(false), habitat(mother->habitat)
{

    // Recombination and transmission of genes
    inheritGamete(mother, parameters, geneticArchitecture);
    inheritGamete(father, parameters, geneticArchitecture);

    // Mutation and development
    mutate(parameters);
    develop(parameters, geneticArchitecture);
}
*/

/*
/// Function to reset one's habitat
void Individual::resetHabitat(const size_t &newHabitat) const
{
    habitat = newHabitat;
}
*/


// Getters

/*
bool Individual::isFemale(const bool &isFemaleHeterogamy) const
{
    return isHeterogamous == isFemaleHeterogamy;
}

Individual::Locus Individual::getLocus(const size_t &locus) const
{
    return genotypes[locus];
}
*/


// Setters

/*
void Individual::disperse(const size_t &nHabitat) const
{
    habitat = (habitat + 1u) % nHabitat;
}

void Individual::setFitness (const std::pair<double, double> &resources) const
{
    fitness = attackRates.first * resources.first + attackRates.second *
     resources.second;
}

void Individual::setBurninFitness(const std::pair<double, double> &resources,
 const double &ecoSelCoeff) const
{
    setFitness(resources);
    fitness *= exp(-ecoSelCoeff * sqr(phenotypes[1u]));
}

void Individual::setAttackRates(const double &ecoSelCoeff)
{
    attackRates.first  = exp(-ecoSelCoeff * sqr(phenotypes[0u] + 1.0));
    attackRates.second = exp(-ecoSelCoeff * sqr(phenotypes[0u] - 1.0));
}

void Individual::setMatePreference(const double &matePreferenceStrength)
{
    matePreference = matePreferenceStrength * phenotypes[1u];
}

void Individual::chooseMates(const double &matingSeasonEnd,
        std::discrete_distribution<size_t> &maleMarket,
        const std::vector<PInd> &males,
        const ParameterSet &parameters) const
{
    // Loop through offspring and through the mating season
    for (size_t t = 0u; nOffspring > 0u && t < matingSeasonEnd; ++t) {

        // Sample a male
        const size_t idFoundMale = maleMarket(rnd::rng);

        // Mate choice
        bool isMating = acceptMate(males[idFoundMale], parameters);

        // Birth
        if (isMating) {

            mates.push_back(idFoundMale);
            --nOffspring;

        }
    }
}

double Individual::assessMatingProb(const double &ecoTraitDistance,
 const double &tiny) const
{
    double matingProb;

    // Assortative or disassortative mating
    if (matePreference >= 0) {
        matingProb = calcAssortProb(matePreference, phenotypes[1u],
         ecoTraitDistance);
    }
    else {
        matingProb = calcDisassortProb(matePreference, phenotypes[1u],
         ecoTraitDistance);
    }

    // Normalize
    matingProb = matingProb < tiny ? 0.0 : matingProb;
    matingProb = matingProb > 1.0 - tiny ? 1.0 : matingProb;
    assert(matingProb >= 0.0 || matingProb <= 1.0);

    return matingProb;

}

bool Individual::acceptMate(Individual const * const male,
 const ParameterSet& parameters) const
{

    // In case of random mating, mate
    if(matePreference == 0.0) return true;

    // Observed male
    const double maleEcoTrait = male->phenotypes[0u];
    const double ecoTraitDistance = sqr(phenotypes[0u] - maleEcoTrait);

    // Mating probability
    double matingProb = assessMatingProb(ecoTraitDistance, parameters.tiny);

    bool isAccepted = rnd::bernoulli(matingProb);
    return isAccepted;

}

size_t Individual::sampleClutchSize(const double &birthRate) const
{
    return rnd::poisson(birthRate * fitness);
}

bool Individual::survive(const double &survivalProb) const
{
    return rnd::bernoulli(survivalProb);
}

void Individual::setEcotype(const std::pair<double, double> &ecotypeBoundary)
 const
{
    ecotype = compareAlongTradeOff(attackRates, ecotypeBoundary) ? 1u : 0u;
}

void Individual::mutate(const ParameterSet& parameters)
{

    // Sample mutations
    size_t nMutations = rnd::poisson(parameters.nBits *
     parameters.mutationRate);

    // Distribute them across the genome
    while (nMutations) {
        size_t i = rnd::random_int(parameters.nBits);
        genomeSequence[i] = !genomeSequence[i];
        --nMutations;
    }
}


// Function to reserve a certain length for vector of genotypes along the genome
void Individual::initializeSizeGenotypeVector(const size_t &nLoci) {

    genotypes.reserve(nLoci);

}

void Individual::initializeSizeIndivTraitSpecificMetrics(const size_t &nTraits)
{

    geneticValues.reserve(nTraits);
    envirValues.reserve(nTraits);
    phenotypes.reserve(nTraits);

}


void Individual::develop(const ParameterSet& parameters,
 const GeneticArchitecture &geneticArchitecture)
{

    initializeSizeGenotypeVector(parameters.nLoci);

    // Set genetic value of all loci
    for (size_t i = 0u; i < parameters.nLoci; ++i) {

        setLocusGeneticValue(i, geneticArchitecture, parameters);

    }

    // Accumulate phenotypic contributions and add environmental effect

    initializeSizeIndivTraitSpecificMetrics(parameters.nTraits);

    for (size_t trait = 0u; trait < parameters.nTraits; ++trait) {
        setGeneticValue(trait, geneticArchitecture);
        setEnvirValue(trait, parameters.scaleE[trait]);
        setPhenotype(trait);
    }

    // Compute attack rates
    setAttackRates(parameters.ecoSelCoeff);

    // Compute mate preference (useful for females only)
    setMatePreference(parameters.matePreferenceStrength);

}

void Individual::expressGene(const size_t &nucleotide, 
        const size_t &locus, 
        const double &scaleD, 
        const double &dominanceCoeff)
{

    // Determine genotype and local effect on the phenotype
    bool isHomozygous = genomeSequence[nucleotide] ==
     genomeSequence[nucleotide + 1u];
    if (isHomozygous) {

        // Homozygote AA
        if (genomeSequence[nucleotide]) {
            genotypes[locus].alleleCount = 2u;
            genotypes[locus].expression = 1.0;
        }

        // Homozygote aa
        else {
            genotypes[locus].alleleCount = 0u;
            genotypes[locus].expression = -1.0;
        }
    }

    else {

        // Heterozygote Aa
        genotypes[locus].alleleCount = 1u;
        genotypes[locus].expression = scaleD * dominanceCoeff;
    }
}

void Individual::setAdditiveValue(const size_t &i, const double &scaleA,
 const double &effectSize)
{
    genotypes[i].locusGeneticValue = scaleA * effectSize *
     genotypes[i].expression;
}

void Individual::setEpistaticValue(const size_t &i, const double &scaleI, const std::list<std::pair<size_t, double> > &neighbors)
{
    // For each interaction
    for (std::pair<size_t, double> edge : neighbors) {

        size_t j = edge.first;

        // Compute interaction strength and distribute phenotypic effect over
        // contributing loci
        double epistaticEffect = 0.5 * scaleI * edge.second *
         genotypes[i].expression * genotypes[j].expression;
        genotypes[i].locusGeneticValue += epistaticEffect;
        genotypes[j].locusGeneticValue += epistaticEffect;
    }
}

void Individual::setPhenotype(const size_t &trait)
{
    phenotypes[trait] = geneticValues[trait] + envirValues[trait];
}

void Individual::setEnvirValue(const size_t &trait, const double &scaleE)
{
    // Add environmental effect
    envirValues[trait] = rnd::normal(0.0, scaleE);
}

void Individual::setGeneticValue(const size_t &trait,
 const GeneticArchitecture &geneticArchitecture)
{
    geneticValues[trait] = 0.0;

    // Accumulate genetic contributions
    for (size_t i : geneticArchitecture.networkVertices[trait]) {
        geneticValues[trait] += genotypes[i].locusGeneticValue;
    }
}

void Individual::setLocusGeneticValue(const size_t &locus,
 const GeneticArchitecture &geneticArchitecture, const ParameterSet &parameters)
{

    size_t nucleotidePos = locus << 1u;  // Nucleotide position
    size_t trait = geneticArchitecture.locusConstants[locus].trait;

    // Express the gene
    expressGene(nucleotidePos, locus, parameters.scaleD[trait], geneticArchitecture.locusConstants[locus].dominanceCoeff);

    // Compute local non-epistatic genetic value
    setAdditiveValue(locus, parameters.scaleA[trait], geneticArchitecture.locusConstants[locus].effectSize);

    // Compute local epistatic value
    setEpistaticValue(locus, parameters.scaleI[trait], geneticArchitecture.locusConstants[locus].neighbors);

}

// Function to set a random genome sequence from a SNP frequency
void Individual::setGenomeSequence(const size_t &nBits, const double &freqSNP)
{
    for (size_t nucleotide = 0u; nucleotide < nBits; nucleotide += 2u) {
        genomeSequence.push_back(nucleotide % 4u == 0u);
        genomeSequence.push_back(nucleotide % 4u == 0u);
        if (rnd::uniform() < freqSNP) {
            genomeSequence[nucleotide] = !genomeSequence[nucleotide];
        }
        if (rnd::uniform() < freqSNP) {
            genomeSequence[nucleotide + 1u] = !genomeSequence[nucleotide + 1u];
        }
    }
}

void Individual::recombineFreely(size_t &haplotype,
        size_t &chromosome,
        const size_t &nChromosomes,
        const double &chromosomeSize,
        double &freeRecombinationPoint) const 
{
    // Recombine by switching to random haplotype
    haplotype = rnd::bernoulli(0.5) ? 0u : 1u;

    // Set next free recombination point
    if (chromosome < nChromosomes - 1u) {
        freeRecombinationPoint = chromosomeSize;
        ++chromosome;
    }
    else {
        freeRecombinationPoint = 1.0;
    }
}

void Individual::crossOver(size_t &haplotype, const double &recombinationRate,
 const double &mapLength, double &crossOverPoint) const
{
    // Cross over to the opposite haplotype
    haplotype = (haplotype + 1u) % 2u;

    // Set next cross-over point
    crossOverPoint += rnd::exponential(recombinationRate * mapLength);
}

void Individual::inheritLocus(Individual const * const parent,
 const bool &isMother, const size_t &locus, const size_t &haplotype)
{
    size_t genomePosition = isMother ? locus << 1 : (locus << 1) + 1u;
    genomeSequence[genomePosition] =
     parent->genomeSequence[(locus << 1u) + haplotype];
}

void Individual::determineSex(const bool &isMother,
 const bool &isFemaleHeterogamy, const size_t &haplotype)
{
    if (isMother) {
        isHeterogamous = haplotype == 0u && isFemaleHeterogamy;
    }
    else {
        isHeterogamous = haplotype == 1u && !isFemaleHeterogamy;
    }
}

void Individual::inheritGamete(Individual const * const parent,
 const ParameterSet &parameters, const GeneticArchitecture &geneticArchitecture)
{

    double freeRecombinationPoint = 0.0;
    double crossOverPoint = 0.0;

    // Loop through loci
    for (size_t locus = 0u, chromosome = 0u, haplotype = 0u;
     locus < parameters.nLoci; ++locus) {

        // Interchromosomal recombination
        bool isFreeRecombination =
         geneticArchitecture.locusConstants[locus].location >
          freeRecombinationPoint;
        if (isFreeRecombination) {
            recombineFreely(haplotype, chromosome, 
                    parameters.nChromosomes, 
                    geneticArchitecture.chromosomeSizes[chromosome], 
                    freeRecombinationPoint);
        }

        // Intrachromosomal recombination
        bool isCrossOver = geneticArchitecture.locusConstants[locus].location >
         crossOverPoint;
        if (isCrossOver) {
            crossOver(haplotype, parameters.recombinationRate,
             parameters.genomeLength, crossOverPoint);
        }

        // Inherit parental haplotype
        bool isMother = parent->isFemale(parameters.isFemaleHeterogamy);
        inheritLocus(parent, isMother, locus, haplotype);

        // Sex determination locus
        const bool isSexDeterminationLocus = locus == 0u;
        if (isSexDeterminationLocus) {
            determineSex(isMother, parameters.isFemaleHeterogamy, haplotype);
        }
    }
}

*/





