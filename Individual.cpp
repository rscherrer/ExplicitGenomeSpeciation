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
#include "Individual.h"
#include "random.h"
#include "ParameterSet.h"
#include "Population.h"
#include "GeneticArchitecture.h"

/*=======================================================================================================
                                         member functions
========================================================================================================*/

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
    }
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

void Individual::setGeneticValue(const size_t &trait, const GeneticArchitecture &geneticArchitecture)
{
    geneticValues[trait] = 0.0;

    // Accumulate genetic contributions
    for (size_t i : geneticArchitecture.networkVertices[trait]) {
        geneticValues[trait] += genotypes[i].locusGeneticValue;
    }
}

void Individual::setEnvirValue(const size_t &trait, const double &scaleE)
{
    // Add environmental effect
    envirValues[trait] = rnd::normal(0.0, scaleE);
}

void Individual::setPhenotype(const size_t &trait)
{
    phenotypes[trait] = geneticValues[trait] + envirValues[trait];
}

void Individual::setViability(const double &costIncompat, const GeneticArchitecture &geneticArchitecture)
{
    if (costIncompat > 0.0) {

        // Initialize the number of incompatibilities
        size_t nIncompatibilities = 0;

        // For each locus underlying trait the neutral trait
        for (size_t i : geneticArchitecture.networkVertices[2u]) {

            // For each interaction
            for (std::pair<size_t, double> edge : geneticArchitecture.locusConstants[i].neighbors) {

                // Record expression of both interacting genes
                size_t j = edge.first;
                double ei = genotypes[i].expression;
                double ej = genotypes[j].expression;

                // If both expression levels are negative, there is an incompatibility
                if (ei < 0.0 && ej < 0.0) {
                    ++nIncompatibilities;
                }
            }
        }

        // Viability is related to the number of incompatibilities
        viability = exp(- nIncompatibilities * costIncompat);

    }
    else {
        viability = 1.0;
    }
}

void Individual::setAttackRates(const double &ecoSelCoeff)
{
    attackRates.first  = exp(-ecoSelCoeff * sqr(phenotypes[0u] + 1.0));
    attackRates.second = exp(-ecoSelCoeff * sqr(phenotypes[0u] - 1.0));
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
        setPhenotype(trait);
    }

    // Compute viability
    setViability(parameters.costIncompat, geneticArchitecture);

    // Compute attack rates
    setAttackRates(parameters.ecoSelCoeff);

}

void Individual::prepareChoice() const
{
    double xi = phenotypes[0u];
    if(!obs.empty()) obs.clear();
    obs.push_back(0.0);
    xsum = xi;
    xxsum = xi * xi;
}

bool Individual::acceptMate(Individual const * const male, const ParameterSet& parameters) const
{

    double scale = parameters.matePreferenceStrength * traitP[1u];
    if(scale == 0.0) return true;

    // observed male
    const double xj = male->traitP[0u];
    const double dij = sqr(traitP[0u] - xj);

    if(parameters.isTypeIIMateChoice) {

        double matingProb = scale >= 0 ? exp(- parameters.matePreferenceStrength * sqr(traitP[1u]) * dij / 2.0) : 1.0 - sqr(sqr(traitP[1u])) * exp(- parameters.matePreferenceStrength * sqr(traitP[1u]) * dij / 2.0);
        matingProb = matingProb < parameters.tiny ? 0.0 : matingProb;
        matingProb = matingProb > 1.0 - parameters.tiny ? 1.0 : matingProb;
        if(matingProb < 0.0 || matingProb > 1.0) throw std::logic_error("mating probability out of bounds");

        return(rnd::bernoulli(matingProb));

    } else {

        // insert observation in sorted list
        std::list<double>::iterator it = std::upper_bound (obs.begin(), obs.end(), dij);
        obs.insert(it, dij);
        const size_t n = obs.size();

        // update statistics
        xsum += xj;
        xxsum += xj * xj;
        const double mu = xsum / n;
        const double var = (xxsum - n * mu * mu) / (n - 1u);

        if(scale < 0.0) {
            // preference towards higher ecotype distance (disassortative mating)

            // reject male if there is no variation in the sample
            if(var < parameters.tiny) return false;
            scale /= var;

            // determine threshold for mate acceptance
            double delta, sum = 0.0, dxy0 = obs.back(), dxy1;
            size_t k = 0u;
            for (std::list<double>::reverse_iterator rit = obs.rbegin();;) {
                ++rit; ++k;  // the first term can be skipped because it has zero weight in the integral
                if(rit != obs.rend()) {
                    dxy1 = *rit;
                    delta = scale * k * (dxy1 - dxy0);
                    if(sum + delta < n) {
                        sum += delta;
                        dxy0 = dxy1;
                    }
                    else break;
                }
                else return true;
            }
            double aux = (n - sum) / delta;
            const double threshold = (1.0 - aux) * dxy0 + aux * dxy1;

            return (dij > threshold);
        }
        else {
            // preference towards lower ecotype distance (assortative mating)

            // accept male if there is no variation in the sample
            if(var < parameters.tiny) return true;
            scale /= var;

            // determine threshold for mate acceptance
            double delta, sum = 0.0, dxy0 = obs.front(), dxy1;
            size_t k = 0u;
            for (std::list<double>::iterator it = obs.begin();;) {
                ++it; ++k;  // the first term can be skipped because it has zero weight in the integral
                if(it != obs.end()) {
                    dxy1 = *it;
                    delta = scale * k * (dxy1 - dxy0);
                    if(sum + delta < n) {
                        sum += delta;
                        dxy0 = dxy1;
                    }
                    else break;
                }
                else return true;
            }
            double aux = (n - sum) / delta;
            const double threshold = (1.0 - aux) * dxy0 + aux * dxy1;

            return (dij < threshold);
        }

    }


}

size_t Individual::setEcotype(const TradeOffPt &threshold) const
{
    return ecotype = tradeOffCompare(attackRate, threshold) ? 2u : 1u;
}

double Individual::getBurnInRpSc(double selCoeff) const
{
    return exp(-selCoeff * sqr(traitP[1u]));
}


