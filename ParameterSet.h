#ifndef EXPLICITGENOMESPECIATION_PARAMETERSET_H
#define EXPLICITGENOMESPECIATION_PARAMETERSET_H

#include <vector>
#include <string>

class ParameterSet {

public:

    // Setters
    void setDefaultSeed();
    void readSeed(std::ifstream&, std::string&);
    void readIsArchitecture(std::ifstream&, std::string&);
    void setScale(std::vector<double>&, std::ifstream&);
    void readInput(const std::string&, std::ifstream&);
    void readParameters(const std::string&);
    void newArchitectureFileName();
    void setArchitectureFileName(const std::string&);

    // Ecological parameters
    size_t  initialPopSize          = 100u;
    double  dispersalRate           = 1.0e-3;
    double  birthRate               = 4.0;
    double  habitatAsymmetry        = 0.5;
    double  survivalProb            = 0.8;
    double  ecoSelCoeff             = 1.0;
    double  matePreferenceStrength  = 10.0;
    double  mateEvaluationCost      = 0.01;
    double  maxResourceCapacity     = 5000.0;
    double  maxResourceGrowth       = 1.0;

    // Genetic parameters
    size_t nEcoLoci         = 400u;
    size_t nMatLoci         = 200u;
    size_t nNtrLoci         = 400u;
    size_t nEcoInteractions = 1000u;
    size_t nMatInteractions = 500u;
    size_t nNtrInteractions = 0u;
    size_t nChromosomes     = 3u;
    size_t nLoci = nEcoLoci + nMatLoci + nNtrLoci;
    size_t nBits = 2u * nLoci;
    std::string strInitialSequence;
    std::vector<bool> initialSequence;
    double  freqSNP                 = 0.02;
    double  mutationRate            = 1.0e-5;
    double  genomeLength            = 300.0;
    bool    isFemaleHeterogamy      = false;
    double  recombinationRate       = 0.01;

    // Genotype-phenotype map
    bool isGenerateArchitecture = true;
    std::string architectureFileName;
    double  networkSkewness = 1.0;
    std::vector<double> scaleA {1.0, 1.0, 1.0};
    std::vector<double> scaleD {0.0, 0.0, 0.0};
    std::vector<double> scaleI {0.0, 0.0, 0.0};
    std::vector<double> scaleE {0.0, 0.0, 0.0};
    double shapeEffectSizes = 2.0;
    double scaleEffectSizes = 1.0;
    double shapeInteractionWeights = 5.0;
    double scaleInteractionWeights = 1.0;

    // Simulation parameters
    int  tBurnIn                 = 1;
    int  tEndSim                 = 5;
    int  tGetDat                 = 1;
    int  tSavDat                 = 1;
    double tiny                  = 1.0e-12;    // for clipping towards zero
    size_t seed;

    // Unsettables
    size_t  nHabitats               = 2u;
    size_t  nTraits                 = 3u;

private:

};

#endif
