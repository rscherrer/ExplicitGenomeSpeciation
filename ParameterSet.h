//
// Created by p278834 on 7-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_PARAMETERSET_H
#define EXPLICITGENOMESPECIATION_PARAMETERSET_H

#include <array>

class ParameterSet {

public:

    const size_t nIndividualInit = 100u;

    double freqSNP = 0.02;

    std::array<double, 3u> scaleA {1.0, 1.0, 1.0};
    std::array<double, 3u> scaleD {0.0, 0.0, 0.0};
    std::array<double, 3u> scaleI {0.0, 0.0, 0.0};
    std::array<double, 3u> scaleE {0.0, 0.0, 0.0};

    double  mutationRate            = 1.0e-5;
    double  mapLength               = 300.0;
    bool    isFemaleHeteroGamety    = false;

    double  dispersalRate           = 1.0e-3;
    double  alpha                   = 5.0e-3;
    double  beta                    = 4.0;
    double  habitatAsymmetry        = 0.5;
    double  survivalProb            = 0.8;
    double  ecoSelCoeff             = 1.0;
    double  matePreferenceStrength  = 10.0;
    double  mateEvaluationCost      = 0.01;
    double  costIncompat            = 0.0;
    double  networkSkewness         = 1.0;

    bool isTypeIIResourceUtilisation = true;
    bool isTypeIIMateChoice = true;

    int  tBurnIn                 = 1;
    int  tEndSim                 = 5;
    int  tGetDat                 = 1;
    int  tSavDat                 = 1;

    unsigned int seed;
    bool generateArchitecture;
    std::string architecture, sequence;

    const size_t nEcoLoci         = 400u;
    const size_t nMatLoci         = 200u;
    const size_t nNtrLoci         = 400u;
    const size_t nEcoInteractions = 1000u;
    const size_t nMatInteractions = 500u;
    const size_t nNtrInteractions = 0u;
    const size_t nChromosomes     = 3u;
    const size_t nHabitat         = 2u;
    const size_t nCharacter       = 3u;
    const double tiny             = 1.0e-12;    // for clipping towards zero
    const size_t nLoci = nEcoLoci + nMatLoci + nNtrLoci;
    const size_t nBits = 2u * nLoci;

//private:
    //ParameterSet(const ParameterSet&) = delete;
};


#endif //EXPLICITGENOMESPECIATION_PARAMETERSET_H
