#ifndef EXPLICITGENOMESPECIATION_PARAM_H
#define EXPLICITGENOMESPECIATION_PARAM_H

#include "Types.h"
//#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>
//#include <sstream>
#include <cassert>
#include <cstdint>

class Param {

public:

    Param() : seed(makeDefaultSeed())
    {
        capEdges();
        checkParams();
    }

    /// Getters
    size_t getNChromosomes() const { return nChromosomes; }
    size_t getNLoci() const { return nLoci; }
    vecUns getNLociPerTrait() const { return nLociPerTrait; }
    vecUns getNEdgesPerTrait() const { return nEdgesPerTrait; }
    vecDbl getSkewnesses() const { return skewnesses; }
    double getEffectSizeShape() const { return effectSizeShape; }
    double getEffectSizeScale() const { return effectSizeScale; }
    double getInteractionWeightShape() const { return interactionWeightShape; }
    double getInteractionWeightScale() const { return interactionWeightScale; }
    double getDominanceVariance() const { return dominanceVariance; }
    double getSNPFreq() const { return freqSNP; }
    size_t getSeed() const { return seed; }
    int getTEndSim() const { return tEndSim; }
    int getTSave() const { return tSave; }
    int getTBurnIn() const { return tBurnIn; }
    size_t getInitialPopSize() const { return initialPopSize; }
    double getDispersalRate() const { return dispersalRate; }
    double getSurvivalProb() const { return survivalProb; }
    double getBirthRate() const { return birthRate; }
    double getEcoSelCoeff() const { return ecoSelCoeff; }
    double getMatePreferenceStrength() const { return matePreferenceStrength; }
    double getMateEvaluationCost() const { return mateEvaluationCost; }
    double getMaxResourceCapacity() const { return maxResourceCapacity; }
    double getMaxResourceGrowth() const { return maxResourceGrowth; }
    double getHabitatSymmetry() const { return habitatSymmetry; }
    double getRecombinationRate() const { return recombinationRate; }
    double getMutationRate() const { return mutationRate; }
    double getMaxFeedingRate() const { return maxFeedingRate; }
    bool getIsFemaleHeterogamy() const { return isFemaleHeterogamy; }
    bool getRecord() const { return record; }
    vecDbl getScaleA() const { return scaleA; }
    vecDbl getScaleD() const { return scaleD; }
    vecDbl getScaleI() const { return scaleI; }
    vecDbl getScaleE() const { return scaleE; }

    /// Setters
    void readParams(std::ifstream&);
    void setTEndSim(const size_t &t) { tEndSim = t; }
    void setTSave(const size_t &t) { tSave = t; }
    void setInitialPopSize(const size_t &n) { initialPopSize = n; }
    void setDispersalRate(const double &d) { dispersalRate = d; }
    void setSurvivalProb(const double &p) { survivalProb = p; }
    void setBirthRate(const double &b) { birthRate = b; }
    void setMatePreferenceStrength(const double &);
    void setNChromosomes(const size_t &nchrom) { nChromosomes = nchrom; }
    void setNLociPerTrait(const vecUns&);
    void setNEdgesPerTrait(const vecUns&);
    void setSkewnesses(const vecDbl &skews) { skewnesses = skews; }
    void setDominanceVariance(const double &x) { dominanceVariance = x; }
    void setNLoci(const size_t &nloci) { nLoci = nloci; }
    void setSeed(const size_t &number) { seed = number; }
    void capEdges();
    void setEffectSizeScale(const double &x) { effectSizeScale = x; }
    void setInteractionWeightScale(const double&);

    /// Makers
    size_t makeDefaultSeed();

    /// Checkers
    void checkParams();

private:

    // Ecological parameters
    double maxResourceCapacity     = 500.0;
    double maxResourceGrowth       = 1.0;
    double habitatSymmetry         = 1.0;
    double ecoSelCoeff             = 1.0;
    size_t initialPopSize          = 100u;
    double dispersalRate           = 1.0e-3;
    double birthRate               = 2.0;
    double survivalProb            = 0.6;
    double matePreferenceStrength  = 10.0;
    double mateEvaluationCost      = 0.01;
    double maxFeedingRate          = 4.0E-4;

    // Genetic parameters
    size_t nEcoLoci = 10;
    size_t nMatLoci = 10u;
    size_t nNtrLoci = 10u;
    size_t nEcoEdges = 0u;
    size_t nMatEdges = 0u;
    size_t nNtrEdges = 0u;
    size_t nChromosomes = 3u;

    size_t nLoci = nEcoLoci + nMatLoci + nNtrLoci;
    vecUns nLociPerTrait = { nEcoLoci, nMatLoci, nNtrLoci };
    vecUns nEdgesPerTrait = { nEcoEdges, nMatEdges, nNtrEdges };

    double  mutationRate            = 1.0e-5;
    double  recombinationRate       = 0.01;
    double  freqSNP                 = 0.5;
    bool    isFemaleHeterogamy      = false;

    // Genotype-phenotype map
    vecDbl scaleA = {1.0, 1.0, 1.0};
    vecDbl scaleD = {0.0, 0.0, 0.0};
    vecDbl scaleI = {0.0, 0.0, 0.0};
    vecDbl scaleE = {0.0, 0.0, 0.0};
    vecDbl skewnesses = { 1.0, 1.0, 1.0 };
    double effectSizeShape = 2.0;
    double effectSizeScale = 1.0;
    double interactionWeightShape = 5.0;
    double interactionWeightScale = 1.0;
    double dominanceVariance = 1.0;

    // Simulation parameters
    int  tBurnIn = 10;
    int  tEndSim = 10;
    int  tSave = 10;
    bool record = true;
    size_t seed;

};

#endif
