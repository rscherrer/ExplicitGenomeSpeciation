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
#include "individual.h"
#include "random.h"

extern double mutationRate, mapLength, ecoSelCoeff, matePreferenceStrength, costIncompat;
extern bool isTypeIIMateChoice;
extern std::array<double, nCharacter> scaleA, scaleD, scaleI, scaleE;

/*=======================================================================================================
                                         member functions
========================================================================================================*/

Individual::Individual(double freqSNP) :
isHeteroGamous(rnd::bernoulli(0.5)), habitat(0u), ecotype(0u)
// default constructor; called on initialisation
{
    // initial genotype
    for(size_t i = 0u; i < nBits; i += 2u) {
        genome[i] = genome[i + 1u] = (i % 4u == 0u);
        if(rnd::uniform() < freqSNP) genome.flip(i);
        if(rnd::uniform() < freqSNP) genome.flip(i + 1u);
    }
    mutate();
    develop();
}

Individual::Individual(const std::string &sequence) :
genome(sequence), isHeteroGamous(rnd::bernoulli(0.5)), habitat(0u), ecotype(0u)
// default constructor; called on initialisation
{
    mutate();
    develop();
}

Individual::Individual(Individual const * const mother, Individual const * const father) :
    isHeteroGamous(false), habitat(mother->habitat)
// constructor implementing sexual reproduction
{
    // transmission of genes from mother
    double freeRecombinationPoint = 0.0, crossOverPoint = 0.0;
    for(size_t i = 0u, lnkgr = 0u, hpltp = 0u; i < nLoci; ++i) {
        if(characterLocus[i].location > freeRecombinationPoint) {
            // switch to random haplotype
            hpltp = (rnd::bernoulli(0.5) ? 0u : 1u);
            // set next free recombination point
            if(lnkgr < nChromosomes - 1u) {
                freeRecombinationPoint = chromosomeSize[lnkgr];
                ++lnkgr;
            }
            else freeRecombinationPoint = 1.0;
        }
        if(characterLocus[i].location > crossOverPoint) {
            // cross over to other haplotype
            hpltp = (hpltp + 1u) % 2u;
            // set next cross-over point
            crossOverPoint += rnd::exponential(0.01 * mapLength);
        }
        genome[i << 1] = mother->genome[(i << 1u) + hpltp];
        if(i == 0u && hpltp == 0u && isFemaleHeteroGamety) isHeteroGamous = true;
    }
    
    // transmission of genes from father
    freeRecombinationPoint = crossOverPoint = 0.0;
    for(size_t i = 0u, lnkgr = 0u, hpltp = 0u; i < nLoci; ++i) {
        if(characterLocus[i].location > freeRecombinationPoint) {
            // switch to random haplotype
            hpltp = (rnd::bernoulli(0.5) ? 0u : 1u);
            // set next free recombination point
            if(lnkgr < nChromosomes - 1u) {
                freeRecombinationPoint = chromosomeSize[lnkgr];
                ++lnkgr;
            }
            else freeRecombinationPoint = 1.0;
        }
        if(characterLocus[i].location > crossOverPoint) {
            // cross over to other haplotype
            hpltp = (hpltp + 1u) % 2u;
            // set next cross-over point
            crossOverPoint += rnd::exponential(0.01 * mapLength);
        }
        genome[(i << 1) + 1u] = father->genome[(i << 1u) + hpltp];
        if(i == 0u && hpltp == 1u && !isFemaleHeteroGamety) isHeteroGamous = true;
    }
    mutate();
    develop();
}

void Individual::mutate()
// implements mutation
{
    size_t k = rnd::poisson(nBits * mutationRate);
    while(k) {
        size_t i = rnd::random_int(nBits);
        genome.flip(i);
        --k;
    }
}

void Individual::develop()
// implements genotype->phenotype map
{
    // additive component
    for(size_t i = 0u; i < nLoci; ++i) {
        
        size_t k = i << 1u;
        size_t crctr = characterLocus[i].character;
        
        // determine genotype and locus phenotypic effect
        if(genome[k] == genome[k + 1u]) {
            if(genome[k]) {         // homozygote AA
                traitLocus[i].alleleCount = 2u;
                traitLocus[i].expression = +1.0;
            }
            else {                  // homozygote aa
                traitLocus[i].alleleCount = 0u;
                traitLocus[i].expression = -1.0;
            }
        }
        else {                      // heterozygote Aa
            traitLocus[i].alleleCount = 1u;
            traitLocus[i].expression =
                scaleD[crctr] * characterLocus[i].dominanceCoeff;
        }
        traitLocus[i].geneticValue = scaleA[crctr] *
            characterLocus[i].effectSize * traitLocus[i].expression;
    }
    // epistatic interactions
    for(size_t i = 0u; i < nLoci; ++i) {
        size_t crctr = characterLocus[i].character;
        for(std::pair<size_t, double> edge : characterLocus[i].edges) {
            size_t j = edge.first;
            
            // compute interaction strength and distribute phenotypic effect over contributing loci
            double Iij = 0.5 * scaleI[crctr] * edge.second *
                traitLocus[i].expression * traitLocus[j].expression;
            traitLocus[i].geneticValue += Iij;
            traitLocus[j].geneticValue += Iij;
        }
    }
    // accumulate phenotypic contributions and add environmental effect
    for(size_t crctr = 0u; crctr < nCharacter; ++crctr) {
        traitG[crctr] = 0.0;
        for(size_t i : vertices[crctr])
            traitG[crctr] += traitLocus[i].geneticValue;
        traitE[crctr] = rnd::normal(0.0, scaleE[crctr]);
        traitP[crctr] = traitG[crctr] + traitE[crctr];
    }

    // compute viability
    if(costIncompat > 0.0) {
        // initialize the number of incompatibilities
        size_t nIncompatibilities = 0;
        // for each vertex underlying trait Z
        for(size_t i : vertices[2u]) {
            // for each of its edges
            for(std::pair<size_t, double> edge : characterLocus[i].edges) {
                // record expression level of i and j
                size_t j = edge.first;
                double ei = traitLocus[i].expression;
                double ej = traitLocus[j].expression;
                // if both expression levels are negative, there is an incompatibility
                if(ei < 0.0 && ej < 0.0) {
                    ++nIncompatibilities;
                }
            }
        }
        // compute viability
        viability = exp(- nIncompatibilities * costIncompat);
    }
    else {
        viability = 1.0;
    }
    
    // compute attack rate
    attackRate.first  = exp(-ecoSelCoeff * sqr(traitP[0u] + 1.0));
    attackRate.second = exp(-ecoSelCoeff * sqr(traitP[0u] - 1.0));

}

void Individual::prepareChoice() const
{
    double xi = traitP[0u];
    if(!obs.empty()) obs.clear();
    obs.push_back(0.0);
    xsum = xi;
    xxsum = xi * xi;
}

bool Individual::acceptMate(Individual const * const male) const
{

    double scale = matePreferenceStrength * traitP[1u];
    if(scale == 0.0) return true;

    // observed male
    const double xj = male->traitP[0u];
    const double dij = sqr(traitP[0u] - xj);

    if(isTypeIIMateChoice) {

        double matingProb = scale >= 0 ? exp(- matePreferenceStrength * sqr(traitP[1u]) * dij / 2.0) : 1.0 - sqr(sqr(traitP[1u])) * exp(- matePreferenceStrength * sqr(traitP[1u]) * dij / 2.0);
        matingProb = matingProb < tiny ? 0.0 : matingProb;
        matingProb = matingProb > 1.0 - tiny ? 1.0 : matingProb;
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
            if(var < tiny) return false;
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
            if(var < tiny) return true;
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

size_t Individual::setEcotype(const Individual::TradeOffPt &threshold) const
{
    return ecotype = tradeOffCompare(attackRate, threshold) ? 2u : 1u;
}

double Individual::getBurnInRpSc(double selCoeff) const
{
    return exp(-selCoeff * sqr(traitP[1u]));
}


