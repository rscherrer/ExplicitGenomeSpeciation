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

Individual::Individual(const ParameterSet& parameters) :
isHeteroGamous(rnd::bernoulli(0.5)), habitat(0u), ecotype(0u)
{
    // Generate a genotype
    for (size_t i = 0u; i < parameters.nBits; i += 2u) {
        genomeSequence[i] = genomeSequence[i + 1u] = (i % 4u == 0u);
        if (rnd::uniform() < parameters.freqSNP) {
            genomeSequence[i] = !genomeSequence[i];
        }
        if (rnd::uniform() < parameters.freqSNP) {
            genomeSequence[i + 1u] = !genomeSequence[i + 1u];
        }
    }
    mutate(parameters);
    develop(parameters, genome);
}

Individual::Individual(const std::vector<bool>& sequence, const ParameterSet& parameters) :
genomeSequence(sequence), isHeteroGamous(rnd::bernoulli(0.5)), habitat(0u), ecotype(0u)
{
    mutate(parameters);
    develop(parameters, genome);
}

Individual::Individual(Individual const * const mother, Individual const * const father, const ParameterSet& parameters, const GeneticArchitecture &geneticArchitecture) :
    isHeteroGamous(false), habitat(mother->habitat)
{

    // Inheritance from mom
    double freeRecombinationPoint = 0.0;
    double crossOverPoint = 0.0;

    // Loop through loci
    for (size_t i = 0u, lg = 0u, haplotype = 0u; i < parameters.nLoci; ++i) {

        // Free interchromosomal recombination
        if (geneticArchitecture.locusConstants[i].location > freeRecombinationPoint) {

            // Recombinate by switching to random haplotype
            haplotype = rnd::bernoulli(0.5) ? 0u : 1u;

            // Set next free recombination point
            if (lg < parameters.nChromosomes - 1u) {
                freeRecombinationPoint = geneticArchitecture.chromosomeSizes[lg];
                ++lg;
            }
            else {
                freeRecombinationPoint = 1.0;
            }
        }

        // Intrachromosomal recombination
        if (geneticArchitecture.locusConstants[i].location > crossOverPoint) {

            // Cross over to the opposite haplotype
            haplotype = (haplotype + 1u) % 2u;

            // Set next cross-over point
            crossOverPoint += rnd::exponential(0.01 * parameters.mapLength);
        }

        // Inherit maternal haplotype
        genomeSequence[i << 1] = mother->genomeSequence[(i << 1u) + haplotype];
        if (i == 0u && haplotype == 0u && parameters.isFemaleHeteroGamety) {
            isHeteroGamous = true;
        }
    }
    
    // Inheritance from dad
    freeRecombinationPoint = 0.0;
    crossOverPoint = 0.0;

    // Loop through loci
    for (size_t i = 0u, lg = 0u, haplotype = 0u; i < parameters.nLoci; ++i) {

        // Free interchromosomal recombination
        if(geneticArchitecture.locusConstants[i].location > freeRecombinationPoint) {

            // Recombinate by switching to random haplotype
            haplotype = (rnd::bernoulli(0.5) ? 0u : 1u);

            // Set next free recombination point
            if(lg < parameters.nChromosomes - 1u) {
                freeRecombinationPoint = geneticArchitecture.chromosomeSizes[lg];
                ++lg;
            }

            else freeRecombinationPoint = 1.0;
        }

        // Intrachromosomal recombination
        if (geneticArchitecture.locusConstants[i].location > crossOverPoint) {

            // Cross over to the opposite haplotype
            haplotype = (haplotype + 1u) % 2u;

            // Set next cross-over point
            crossOverPoint += rnd::exponential(0.01 * parameters.mapLength);
        }

        // Inherit paternal haplotype
        genomeSequence[(i << 1) + 1u] = father->genomeSequence[(i << 1u) + haplotype];
        if(i == 0u && haplotype == 1u && !parameters.isFemaleHeteroGamety) isHeteroGamous = true;
    }

    mutate(parameters);
    develop(parameters, genome);
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

    // Non-epistatic component (additive and dominance) of each locus

    for (size_t i = 0u; i < parameters.nLoci; ++i) {
        
        size_t k = i << 1u;
        size_t crctr = geneticArchitecture.locusConstants[i].character;

        expressGene();

        // Determine genotype and local effect on the phenotype
        if (genomeSequence[k] == genomeSequence[k + 1u]) {

            // Homozygote AA
            if (genomeSequence[k]) {
                traitLocus[i].alleleCount = 2u;
                traitLocus[i].expression = 1.0;
            }

            // Homozygote aa
            else {
                traitLocus[i].alleleCount = 0u;
                traitLocus[i].expression = -1.0;
            }
        }

        // Heterozygote Aa
        else {
            traitLocus[i].alleleCount = 1u;
            traitLocus[i].expression =
                    parameters.scaleD[crctr] * geneticArchitecture.locusConstants[i].dominanceCoeff;
        }

        // Compute non-epistatic local genetic value
        computeGeneNonEpistaticValue();

        traitLocus[i].geneticValue = parameters.scaleA[crctr] *
                geneticArchitecture.locusConstants[i].effectSize * traitLocus[i].expression;
    }

    // Epistatic component of each locus
    for (size_t i = 0u; i < parameters.nLoci; ++i) {

        size_t crctr = geneticArchitecture.locusConstants[i].character;

        computeGeneEpistaticValue();

        // For each interaction
        for (std::pair<size_t, double> edge : geneticArchitecture.locusConstants[i].neighbors) {

            size_t j = edge.first;
            
            // Compute interaction strength and distribute phenotypic effect over contributing loci
            double epistaticEffect = 0.5 * parameters.scaleI[crctr] * edge.second *
                traitLocus[i].expression * traitLocus[j].expression;
            traitLocus[i].geneticValue += epistaticEffect;
            traitLocus[j].geneticValue += epistaticEffect;
        }
    }

    // Accumulate phenotypic contributions and add environmental effect
    for (size_t crctr = 0u; crctr < parameters.nCharacter; ++crctr) {

        setGeneticValue();

        traitG[crctr] = 0.0;

        // Accumulate genetic contributions
        for (size_t i : geneticArchitecture.networkVertices[crctr]) {
            traitG[crctr] += traitLocus[i].geneticValue;
        }

        setEnvirValue();

        // Add environmental effect
        traitE[crctr] = rnd::normal(0.0, parameters.scaleE[crctr]);

        setPhenotypeValue();
        traitP[crctr] = traitG[crctr] + traitE[crctr];
        
    }

    // Compute viability
    computeViability();

    if (parameters.costIncompat > 0.0) {

        // Initialize the number of incompatibilities
        size_t nIncompatibilities = 0;

        // For each locus underlying trait the neutral trait
        for (size_t i : geneticArchitecture.networkVertices[2u]) {

            // For each interaction
            for (std::pair<size_t, double> edge : geneticArchitecture.locusConstants[i].neighbors) {

                // Record expression of both interacting genes
                size_t j = edge.first;
                double ei = traitLocus[i].expression;
                double ej = traitLocus[j].expression;

                // If both expression levels are negative, there is an incompatibility
                if (ei < 0.0 && ej < 0.0) {
                    ++nIncompatibilities;
                }
            }
        }

        // Viability is related to the number of incompatibilities
        viability = exp(- nIncompatibilities * parameters.costIncompat);

    }
    else {
        viability = 1.0;
    }
    
    // Compute attack rate
    setAttackRates();

    attackRate.first  = exp(-parameters.ecoSelCoeff * sqr(traitP[0u] + 1.0));
    attackRate.second = exp(-parameters.ecoSelCoeff * sqr(traitP[0u] - 1.0));

}

void Individual::prepareChoice() const
{
    double xi = traitP[0u];
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


