#include "ParameterSet.h"
#include <iostream>
#include <chrono>
#include <sstream>

/// Constructor and member initialization list
ParameterSet::ParameterSet() : initialPopSize(100u),
                                dispersalRate(1.0e-3),
                                birthRate(4.0),
                                habitatAsymmetry(0.5),
                                survivalProb(0.8),
                                ecoSelCoeff(1.0),
                                matePreferenceStrength(10.0),
                                mateEvaluationCost(0.01),
                                maxResourceCapacity(5000.0),
                                maxResourceGrowth(1.0),
                                nEcoLoci(400u),
                                nMatLoci(200u),
                                nNtrLoci(400u),
                                nEcoInteractions(1000u),
                                nMatInteractions(500u),
                                nNtrInteractions(0u),
                                nChromosomes(3u),
                                nLoci(nEcoLoci + nMatLoci + nNtrLoci),
                                nLociPerTrait({nEcoLoci, nMatLoci, nNtrLoci}),
                                nBits(2u * nLoci),
                                freqSNP(0.02),
                                mutationRate(1.0e-5),
                                genomeLength(300u),
                                isFemaleHeterogamy(false),
                                recombinationRate(0.01),
                                isGenerateArchitecture(true),
                                networkSkewness(1.0),
                                scaleA({1.0, 1.0, 1.0}),
                                scaleD({0.0, 0.0, 0.0}),
                                scaleI({0.0, 0.0, 0.0}),
                                scaleE({0.0, 0.0, 0.0}),
                                shapeEffectSizes(2.0),
                                scaleEffectSizes(1.0),
                                shapeInteractionWeights(5.0),
                                scaleInteractionWeights(1.0),
                                tBurnIn(1),
                                tEndSim(5),
                                tGetDat(1),
                                tSavDat(1),
                                tiny(1.0e-12),
                                nHabitats(2u),
                                nTraits(3u)
{}

template <class T>
void ParameterSet::setParameter(T &parameter, std::ifstream &inputFile)
{
    inputFile >> parameter;
}

// Setters

void ParameterSet::setDefaultSeed()
{
    seed = static_cast<size_t>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
}

std::vector<bool> string2bits(const std::string &stringSequence)
{
    std::vector<bool> bitSequence;
    for (auto a : stringSequence) {
        bitSequence.push_back(a == '1');
    }
    return bitSequence;
}

void ParameterSet::setScale(std::vector<double> &scaleVector, std::ifstream &inputFile)
{
    for (size_t trait = 0u; trait < nTraits; ++trait) {
        inputFile >> scaleVector[trait];
    }
}

void ParameterSet::readInput(const std::string &input, std::ifstream &inputFile) {
    bool isSeedProvided = false;
    if (input == "rng_seed_user") {
        inputFile >> seed;
        isSeedProvided = true;
    } else if (input == "architecture_generate")
    {
        isGenerateArchitecture = true;
    } else if (input == "architecture_load") {
        isGenerateArchitecture = false;
        inputFile >> architectureFileName;
    }
    if (input == "initialPopSize") {
        inputFile >> initialPopSize;
    } else if (input == "dispersalRate") {
        inputFile >> dispersalRate;
    } else if (input == "birthRate") {
        inputFile >> birthRate;
    } else if (input == "habitatAsymmetry") {
        inputFile >> habitatAsymmetry;
    } else if (input == "survivalProb") {
        inputFile >> survivalProb;
    } else if (input == "ecoSelCoeff") {
        inputFile >> ecoSelCoeff;
    } else if (input == "matePreferenceStrength") {
        inputFile >> matePreferenceStrength;
    } else if (input == "mateEvaluationCost") {
        inputFile >> mateEvaluationCost;
    } else if (input == "maxResourceCapacity") {
        inputFile >> maxResourceCapacity;
    } else if (input == "maxResourceGrowth") {
        inputFile >> maxResourceGrowth;
    } else if (input == "nEcoLoci") {
        inputFile >> nEcoLoci;
    } else if (input == "nMatLoci") {
        inputFile >> nMatLoci;
    } else if (input == "nNtrLoci") {
        inputFile >> nNtrLoci;
    } else if (input == "nEcoInteractions") {
        inputFile >> nEcoInteractions;
    } else if (input == "nMatInteractions") {
        inputFile >> nMatInteractions;
    } else if (input == "nNtrInteractions") {
        inputFile >> nNtrInteractions;
    } else if (input == "nChromosomes") {
        inputFile >> nChromosomes;
    } else if (input == "initialSequence") {
        inputFile >> strInitialSequence;
        initialSequence = string2bits(strInitialSequence);
    } else if (input == "freqSNP") {
        inputFile >> freqSNP;
    } else if (input == "mutationRate") {
        inputFile >> mutationRate;
    } else if (input == "genomeLength") {
        inputFile >> genomeLength;
    } else if (input == "isFemaleHeterogamy") {
        inputFile >> isFemaleHeterogamy;
    } else if (input == "recombinationRate") {
        inputFile >> recombinationRate;
    } else if (input == "networkSkewness") {
        inputFile >> networkSkewness;
    } else if (input == "scaleA") {
        setScale(scaleA, inputFile);
    } else if (input == "scaleD") {
        setScale(scaleD, inputFile);
    } else if (input == "scaleI") {
        setScale(scaleI, inputFile);
    } else if (input == "scaleE") {
        setScale(scaleE, inputFile);
    } else if (input == "shapeEffectSizes") {
        inputFile >> shapeEffectSizes;
    } else if (input == "scaleEffectSizes") {
        inputFile >> scaleEffectSizes;
    } else if (input == "shapeInteractionWeights") {
        inputFile >> shapeInteractionWeights;
    } else if (input == "scaleInteractionWeights") {
        inputFile >> scaleInteractionWeights;
    } else if (input == "tBurnIn") {
        inputFile >> tBurnIn;
    } else if (input == "tEndSim") {
        inputFile >> tEndSim;
    } else if (input == "tGetDat") {
        inputFile >> tGetDat;
    } else if (input == "tSavDat") {
        inputFile >> tSavDat;
    } else if (input == "tiny") {
        inputFile >> tiny;
    } else {
        throw std::runtime_error("Unknown parameter " + input);
    }
    if (isSeedProvided) {
        std::clog << "Custom seed was provided by the user\n";
    }
    else
    {
        std::clog << "No custom seed was provided, using clock seed instead\n";
    }
}

// Function to read parameters from a parameter file
void ParameterSet::readParameters(std::ifstream& inputFile)
{
    std::string input;
    while (inputFile >> input) {
        readInput(input, inputFile);
    }
}

void ParameterSet::newArchitectureFileName()
{
    std::ostringstream oss;
    oss << "architecture_" << seed << ".txt";
    setArchitectureFileName(oss.str());
}


void ParameterSet::setArchitectureFileName(const std::string &filename)
{
    architectureFileName = filename;
}