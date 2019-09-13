#ifndef EXPLICITGENOMESPECIATION_PARAMETERSET_H
#define EXPLICITGENOMESPECIATION_PARAMETERSET_H

#include <vector>
#include <fstream>
#include <stddef.h>

typedef std::vector<size_t> vecUns;
typedef std::vector<double> vecDbl;

class ParameterSet {

public:

    /// Constructor
    ParameterSet();

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
    size_t getSeed() const { return seed; }
    size_t getTEndSim() const { return tEndSim; }
    size_t getTSave() const { return tSave; }
    size_t getTBurnIn() const { return tBurnIn; }
    size_t getInitialPopSize() const { return initialPopSize; }
    double getDispersalRate() const { return dispersalRate; }
    double getSurvivalProb() const { return survivalProb; }
    double getBirthRate() const { return birthRate; }
    double getMatePreferenceStrength() const { return matePreferenceStrength; }
    double getMateEvaluationCost() const { return mateEvaluationCost; }
    double getMaxResourceCapacity() const { return maxResourceCapacity; }
    double getMaxResourceGrowth() const { return maxResourceGrowth; }
    double getHabitatSymmetry() const { return habitatSymmetry; }
    bool getIsFemaleHeterogamy() const { return isFemaleHeterogamy; }
    bool getRecord() const { return record; }

    /// Setters
    void readParams(std::ifstream&);
    void setTEndSim(const size_t &t) { tEndSim = t; }
    void setTSave(const size_t &t) { tSave = t; }
    void setInitialPopSize(const size_t &n) { initialPopSize = n; }
    void setDispersalRate(const double &d) { dispersalRate = d; }
    void setSurvivalProb(const double &p) { survivalProb = p; }
    void setBirthRate(const double &b) { birthRate = b; }
    void setMatePreferenceStrength(const double &alpha)
    {
        matePreferenceStrength = alpha;
    }
    void setNChromosomes(const size_t &nchrom) { nChromosomes = nchrom; }
    void setNLociPerTrait(const vecUns &locipertrait)
    {
        nLociPerTrait = locipertrait;
        nLoci = nLociPerTrait[0u] + nLociPerTrait[1u] + nLociPerTrait[2u];
        capEdges();
    }
    void setNEdgesPerTrait(const vecUns &edgespertrait)
    {
        nEdgesPerTrait = edgespertrait;
        capEdges();
    }
    void setSkewnesses(const vecDbl &skews) { skewnesses = skews; }
    void setDominanceVariance(const double &x) { dominanceVariance = x; }
    void setNLoci(const size_t &nloci) { nLoci = nloci; }
    void setSeed(const size_t &number) { seed = number; }
    void capEdges();
    void setEffectSizeScale(const double &x) { effectSizeScale = x; }
    void setInteractionWeightScale(const double &x)
    {
        interactionWeightScale = x;
    }

    /// Makers
    size_t makeDefaultSeed();

    /// Checkers
    void checkParams();


private:

    // Ecological parameters
    double  maxResourceCapacity     = 100.0;
    double  maxResourceGrowth       = 1.0;
    double  habitatSymmetry         = 1.0;
    double  ecoSelCoeff             = 1.0;
    size_t  initialPopSize          = 100u;
    double  dispersalRate           = 1.0e-3;
    double  birthRate               = 4.0;
    double  survivalProb            = 0.8;
    double  matePreferenceStrength  = 10.0;
    double  mateEvaluationCost      = 0.01;

    // Genetic parameters
    size_t nEcoLoci = 400u;
    size_t nMatLoci = 200u;
    size_t nNtrLoci = 400u;
    size_t nEcoEdges = 1000u;
    size_t nMatEdges = 500u;
    size_t nNtrEdges = 0u;
    size_t nChromosomes = 3u;

    size_t nLoci = nEcoLoci + nMatLoci + nNtrLoci;
    vecUns nLociPerTrait = { nEcoLoci, nMatLoci, nNtrLoci };
    vecUns nEdgesPerTrait = { nEcoEdges, nMatEdges, nNtrEdges };

    double  mutationRate            = 1.0e-5;
    double  recombinationRate       = 0.01;
    double  genomeLength            = 300.0;
    double  freqSNP                 = 0.02;
    bool    isFemaleHeterogamy      = false;

    // Genotype-phenotype map
    bool isGenerateArchitecture = true;
    std::string architectureFileName = "";
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
    size_t  tBurnIn                 = 1;
    size_t  tEndSim                 = 100;
    size_t  tSave                   = 1;
    double tiny                  = 1.0e-12;
    size_t seed;
    bool record = true;

};

#endif
