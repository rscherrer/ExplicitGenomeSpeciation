#ifndef EXPLICITGENOMESPECIATION_PARAMETERSET_H
#define EXPLICITGENOMESPECIATION_PARAMETERSET_H

#include <vector>
#include <string>

class ParameterSet {

public:

    // Ecological parameters
    size_t  nIndividualInit         = 100u;
    double  dispersalRate           = 1.0e-3;
    double  alpha                   = 5.0e-3;
    double  birthRate               = 4.0;
    double  habitatAsymmetry        = 0.5;
    double  survivalProb            = 0.8;
    double  ecoSelCoeff             = 1.0;
    double  matePreferenceStrength  = 10.0;
    double  mateEvaluationCost      = 0.01;
    size_t  nHabitat                = 2u;
    double  maxResourceCapacity     = 5000.0;
    double  maxResourceGrowth       = 1.0;

    bool isTypeIIResourceUtilisation = true;
    bool isTypeIIMateChoice = true;

    // Genetic parameters
    size_t nEcoLoci         = 400u;
    size_t nMatLoci         = 200u;
    size_t nNtrLoci         = 400u;
    size_t nEcoInteractions = 1000u;
    size_t nMatInteractions = 500u;
    size_t nNtrInteractions = 0u;
    size_t nChromosomes     = 3u;
    size_t nCharacter       = 3u;
    size_t nLoci = nEcoLoci + nMatLoci + nNtrLoci;
    size_t nBits = 2u * nLoci;
    std::string sequenceString;
    std::vector<bool> sequence;
    double  freqSNP                 = 0.02;
    double  mutationRate            = 1.0e-5;
    double  mapLength               = 300.0;
    bool    isFemaleHeterogamety    = false;
    double  recombinationRate       = 0.01;

    // Genotype-phenotype map
    bool isGenerateArchitecture;
    std::string architectureFilename;
    double  networkSkewness = 1.0;
    double  costIncompat = 0.0;
    std::vector<double> scaleA {1.0, 1.0, 1.0};
    std::vector<double> scaleD {0.0, 0.0, 0.0};
    std::vector<double> scaleI {0.0, 0.0, 0.0};
    std::vector<double> scaleE {0.0, 0.0, 0.0};
    double alphaAdditive = 2.0;
    double alphaInteraction = 5.0;

    // Simulation parameters
    int  tBurnIn                 = 1;
    int  tEndSim                 = 5;
    int  tGetDat                 = 1;
    int  tSavDat                 = 1;
    double tiny                  = 1.0e-12;    // for clipping towards zero
    size_t seed;

    // Member functions
    void readParameters(const std::string&);
    void writeParameters(std::ofstream&, const char = ' ');

};


#endif
