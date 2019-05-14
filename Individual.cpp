/*==================================================================================================================================
                                                     individual.cpp
====================================================================================================================================

C++-code accompanying:	
		 
		(ms. in prep).

Written by:
        G. Sander van Doorn
       	Centre for Ecological and Evolutionary Studies - Theoretical Biology Group
        University of Groningen
        the Netherlands

Program version
		xx/xx/2018	:	

Instructions for compiling and running the program
		
	Versions of this program were compiled and run on Windows and Mac, using Microsoft Visual C++
	2010 and XCode, respectively. The code is written in standard C++ and should be compatible with 
	other compilers. 

=================================================================================================================================*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <string>
#include <algorithm>
#include <set>
#include <list>
#include <cassert>
#include "Individual.h"
#include "random.h"
#include "ParameterSet.h"
#include "Population.h"
#include "GeneticArchitecture.h"

/*=======================================================================================================
                                         member functions
========================================================================================================*/

// Constructors

Individual::Individual(const ParameterSet& parameters, const GeneticArchitecture &geneticArchitecture) :
isHeteroGamous(rnd::bernoulli(0.5)), habitat(0u), ecotype(0u)
{
    setGenomeSequence(parameters.nBits, parameters.freqSNP);
    mutate(parameters);
    develop(parameters, geneticArchitecture);
}

Individual::Individual(const std::vector<bool>& sequence, const ParameterSet& parameters, const GeneticArchitecture &geneticArchitecture) :
genomeSequence(sequence), isHeteroGamous(rnd::bernoulli(0.5)), habitat(0u), ecotype(0u)
{
    mutate(parameters);
    develop(parameters, geneticArchitecture);
}

Individual::Individual(Individual const * const mother, Individual const * const father, const ParameterSet& parameters, const GeneticArchitecture &geneticArchitecture) :
        isHeteroGamous(false), habitat(mother->habitat)
{

    // Recombination and transmission of genes
    inheritGamete(mother, parameters, geneticArchitecture);
    inheritGamete(father, parameters, geneticArchitecture);

    // Mutation and development
    mutate(parameters);
    develop(parameters, geneticArchitecture);
}


// Getter

bool Individual::isFemale(const bool &isFemaleHeterogamety) const
{
    return isHeteroGamous == isFemaleHeterogamety;
}


// Setters

void Individual::disperse(const size_t &nHabitat)
{
    habitat = (habitat + 1u) % nHabitat;
}

void Individual::setFitness (const std::pair<double, double> &resources)
{
    fitness = attackRates.first * resources.first + attackRates.second * resources.second;
}

void Individual::setBurninFitness(const std::pair<double, double> &resources, const double &ecoSelCoeff)
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


void Individual::chooseMates(const double &matingSeasonEnd, std::discrete_distribution<size_t> &maleMarket, std::vector<PInd> &males, const ParameterSet &parameters)
{
    // Loop through offspring and through the mating season
    for (size_t t = 0u; nOffspring && t < matingSeasonEnd; ++t) {

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

double Individual::assessMatingProb(const double &ecoTraitDistance, const double &tiny)
{
    double matingProb;

    // Assortative or disassortative mating
    if (matePreference >= 0) {
        matingProb = calcAssortProb(matePreference, phenotypes[1u], ecoTraitDistance);
    }
    else {
        matingProb = calcDisassortProb(matePreference, phenotypes[1u], ecoTraitDistance);
    }

    // Normalize
    matingProb = matingProb < tiny ? 0.0 : matingProb;
    matingProb = matingProb > 1.0 - tiny ? 1.0 : matingProb;
    assert(matingProb >= 0.0 || matingProb <= 1.0);
}

bool Individual::acceptMate(Individual const * const male, const ParameterSet& parameters)
{

    if(matePreference == 0.0) return true;

    // Observed male
    const double maleEcoTrait = male->phenotypes[0u];
    const double ecoTraitDistance = sqr(phenotypes[0u] - maleEcoTrait);

    // Mating probability
    double matingProb = assessMatingProb(ecoTraitDistance, parameters.tiny);

    bool isAccepted = rnd::bernoulli(matingProb);
    return isAccepted;

}

size_t Individual::sampleClutchSize(const double &birthRate)
{
    return rnd::poisson(birthRate * fitness);
}

bool Individual::survive(const double &survivalProb)
{
    return rnd::bernoulli(survivalProb);
}

void Individual::mutate(const ParameterSet& parameters)
{

    // Sample mutations
    size_t nMutations = rnd::poisson(parameters.nBits * parameters.mutationRate);

    // Distribute them across the genome
    while (nMutations) {
        size_t i = rnd::random_int(parameters.nBits);
        genomeSequence[i] = !genomeSequence[i];
        --nMutations;
    }
}


void Individual::develop(const ParameterSet& parameters, const GeneticArchitecture &geneticArchitecture)
{

    // Set genetic value of all loci
    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        setLocusGeneticValue(i, geneticArchitecture, parameters);
    }

    // Accumulate phenotypic contributions and add environmental effect
    for (size_t trait = 0u; trait < parameters.nCharacter; ++trait) {
        setGeneticValue(trait, geneticArchitecture);
        setEnvirValue(trait, parameters.scaleE[trait]);
        setPhenotypes(trait);
    }

    // Compute viability
    setViability(parameters.costIncompat, geneticArchitecture);

    // Compute attack rates
    setAttackRates(parameters.ecoSelCoeff);

    // Compute mate preference (useful for females only)
    setMatePreference(parameters.matePreferenceStrength);

}

void Individual::expressGene(const size_t &nt, const size_t &i, const double &scaleD, const double &dominanceCoeff)
{
    // Determine genotype and local effect on the phenotype
    bool isHomozygous = genomeSequence[nt] == genomeSequence[nt + 1u];
    if (isHomozygous) {

        // Homozygote AA
        if (genomeSequence[nt]) {
            genotypes[i].alleleCount = 2u;
            genotypes[i].expression = 1.0;
        }

            // Homozygote aa
        else {
            genotypes[i].alleleCount = 0u;
            genotypes[i].expression = -1.0;
        }
    }

        // Heterozygote Aa
    else {
        genotypes[i].alleleCount = 1u;
        genotypes[i].expression = scaleD * dominanceCoeff;
    }
}

void Individual::setAdditiveValue(const size_t &i, const double &scaleA, const double &effectSize)
{
    genotypes[i].locusGeneticValue = scaleA * effectSize * genotypes[i].expression;
}

void Individual::setEpistaticValue(const size_t &i, const double &scaleI, const std::list<std::pair<size_t, double> > &neighbors)
{
    // For each interaction
    for (std::pair<size_t, double> edge : neighbors) {

        size_t j = edge.first;

        // Compute interaction strength and distribute phenotypic effect over contributing loci
        double epistaticEffect = 0.5 * scaleI * edge.second * genotypes[i].expression * genotypes[j].expression;
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

void Individual::setGeneticValue(const size_t &trait, const GeneticArchitecture &geneticArchitecture)
{
    geneticValues[trait] = 0.0;

    // Accumulate genetic contributions
    for (size_t i : geneticArchitecture.networkVertices[trait]) {
        geneticValues[trait] += genotypes[i].locusGeneticValue;
    }
}

void Individual::setLocusGeneticValue(const size_t &i,
                                      const GeneticArchitecture &geneticArchitecture,
                                      const ParameterSet &parameters)
{
    size_t nt = i << 1u;  // nucleotide position
    size_t trait = geneticArchitecture.locusConstants[i].character;

    // Express the gene
    expressGene(nt, i, parameters.scaleD[trait], geneticArchitecture.locusConstants[i].dominanceCoeff);

    // Compute local non-epistatic genetic value
    setAdditiveValue(i, parameters.scaleA[trait], geneticArchitecture.locusConstants[i].effectSize);

    // Compute local epistatic value
    setEpistaticValue(i, parameters.scaleI[trait], geneticArchitecture.locusConstants[i].neighbors);

}

void Individual::setGenomeSequence(const size_t &nBits, const double &freqSNP)
{
    for (size_t i = 0u; i < nBits; i += 2u) {
        genomeSequence[i] = genomeSequence[i + 1u] = (i % 4u == 0u);
        if (rnd::uniform() < freqSNP) {
            genomeSequence[i] = !genomeSequence[i];
        }
        if (rnd::uniform() < freqSNP) {
            genomeSequence[i + 1u] = !genomeSequence[i + 1u];
        }
    }
}

void Individual::recombineFreely(size_t &haplotype,
        size_t &lg,
        const size_t &nChromosomes,
        const double &chromosomeSize,
        double &freeRecombinationPoint)
{
    // Recombine by switching to random haplotype
    haplotype = rnd::bernoulli(0.5) ? 0u : 1u;

    // Set next free recombination point
    if (lg < nChromosomes - 1u) {
        freeRecombinationPoint = chromosomeSize;
        ++lg;
    }
    else {
        freeRecombinationPoint = 1.0;
    }
}

void Individual::crossOver(size_t &haplotype, const double &recombinationRate, const double &mapLength, double &crossOverPoint)
{
    // Cross over to the opposite haplotype
    haplotype = (haplotype + 1u) % 2u;

    // Set next cross-over point
    crossOverPoint += rnd::exponential(recombinationRate * mapLength);
}

void Individual::inheritLocus(Individual const * const parent, const bool &isMother, const size_t &locus, const size_t &haplotype)
{
    size_t genomePosition = isMother ? locus << 1 : (locus << 1) + 1u;
    genomeSequence[genomePosition] = parent->genomeSequence[(locus << 1u) + haplotype];
}

void Individual::determineSex(const bool &isMother, const bool &isFemaleHeteroGamety, const size_t &haplotype)
{
    if (isMother) {
        isHeteroGamous = haplotype == 0u && isFemaleHeteroGamety ? true : false;
    }
    else {
        isHeteroGamous = haplotype == 1u && !isFemaleHeteroGamety ? true : false;
    }
}

void Individual::inheritGamete(Individual const * const parent, const ParameterSet &parameters, const GeneticArchitecture &geneticArchitecture)
{

    double freeRecombinationPoint = 0.0;
    double crossOverPoint = 0.0;

    // Loop through loci
    for (size_t i = 0u, lg = 0u, haplotype = 0u; i < parameters.nLoci; ++i) {

        // Interchromosomal recombination
        bool isFreeRecombination = geneticArchitecture.locusConstants[i].location > freeRecombinationPoint;
        if (isFreeRecombination) {
            recombineFreely(haplotype, lg, parameters.nChromosomes, geneticArchitecture.chromosomeSizes[lg], freeRecombinationPoint)
        }

        // Intrachromosomal recombination
        bool isCrossOver = geneticArchitecture.locusConstants[i].location > crossOverPoint;
        if (isCrossOver) {
            crossOver(haplotype, parameters.recombinationRate, parameters.mapLength, crossOverPoint);
        }

        // Inherit parental haplotype
        bool isMother = parent->isFemale(parameters.isFemaleHeteroGamety);
        inheritLocus(parent, isMother, i, haplotype);

        // Sex determination locus
        const bool isSexDeterminationLocus = i == 0u;
        if (isSexDeterminationLocus) {
            determineSex(isMother, parameters.isFemaleHeteroGamety, haplotype);
        }
    }
}


// Accessory functions

bool tradeOffCompare (const std::pair<double, double> &x, const std::pair<double, double> &y)
{
    bool yOnLeft = y.first < y.second;
    if(x.first < x.second) {
        if(yOnLeft) return (x.first < y.first);
        else return true;
    }
    else {
        if(yOnLeft) return false;
        else return (x.second < y.second);
    }
}

double calcAssortProb(const double &matePreference, const double &matingTrait, const double &ecoTraitDistance)
{
    return exp(- matePreference * matingTrait * ecoTraitDistance * 0.5);
}

double calcDisassortProb(const double &matePreference, const double &matingTrait, const double &ecoTraitDistance)
{
    return 1.0 - sqr(sqr(matingTrait)) * exp(- matePreference * matingTrait * ecoTraitDistance * 0.5);
}


// To be taken care of

void Individual::setEcotype(const std::pair<double, double> &threshold) const
{
    ecotype = tradeOffCompare(attackRates, threshold) ? 2u : 1u;
}





