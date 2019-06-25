#include "ParameterSet.h"
#include <iostream>
#include <chrono>

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

void ParameterSet::readSeed(std::ifstream &inputFile, std::string &input)
{
    inputFile >> input;
    if (input == "rng_seed_clock") {
        setDefaultSeed();
    }
    else if (input == "rng_seed_user") {
        inputFile >> seed;
    }
    else {
        throw std::logic_error("\'rng_seed_clock\' or \'rng_seed_user <arg>\' expected at first line of parameterfile\n");
    }
}

void ParameterSet::readIsArchitecture(std::ifstream &inputFile, std::string &input)
{
    inputFile >> input;
    if (input == "architecture_generate") {
        isGenerateArchitecture = true;
    }
    else if (input == "architecture_load") {
        isGenerateArchitecture = false;
        inputFile >> architectureFileName;
    }
    else {
        throw std::logic_error("\'architecture_generate\' or \'architecture_load <arg>\' expected at second line of parameterfile\n");
    }
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

void ParameterSet::readInput(const std::string &input, std::ifstream &inputFile)
{
    if (input == "initialPopSize") {
        inputFile >> initialPopSize;
    }
    else if (input == "dispersalRate") {
        inputFile >> dispersalRate;
    }
    else if (input == "birthRate") {
        inputFile >> birthRate;
    }
    else if (input == "habitatAsymmetry") {
        inputFile >> habitatAsymmetry;
    }
    else if (input == "survivalProb") {
        inputFile >> survivalProb;
    }
    else if (input == "ecoSelCoeff") {
        inputFile >> ecoSelCoeff;
    }
    else if (input == "matePreferenceStrength") {
        inputFile >> matePreferenceStrength;
    }
    else if (input == "mateEvaluationCost") {
        inputFile >> mateEvaluationCost;
    }
    else if (input == "maxResourceCapacity") {
        inputFile >> maxResourceCapacity;
    }
    else if (input == "maxResourceGrowth") {
        inputFile >> maxResourceGrowth;
    }
    else if (input == "nEcoLoci") {
        inputFile >> nEcoLoci;
    }
    else if (input == "nMatLoci") {
        inputFile >> nMatLoci;
    }
    else if (input == "nNtrLoci") {
        inputFile >> nNtrLoci;
    }
    else if (input == "nEcoInteractions") {
        inputFile >> nEcoInteractions;
    }
    else if (input == "nMatInteractions") {
        inputFile >> nMatInteractions;
    }
    else if (input == "nNtrInteractions") {
        inputFile >> nNtrInteractions;
    }
    else if (input == "nChromosomes") {
        inputFile >> nChromosomes;
    }
    else if (input == "initialSequence") {
        inputFile >> strInitialSequence;
        initialSequence = string2bits(strInitialSequence);
    }
    else if (input == "freqSNP") {
        inputFile >> freqSNP;
    }
    else if (input == "mutationRate") {
        inputFile >> mutationRate;
    }
    else if (input == "genomeLength") {
        inputFile >> genomeLength;
    }
    else if (input == "isFemaleHeterogamy") {
        inputFile >> isFemaleHeterogamy;
    }
    else if (input == "recombinationRate") {
        inputFile >> recombinationRate;
    }
    else if (input == "networkSkewness") {
        inputFile >> networkSkewness;
    }
    else if (input == "scaleA") {
        setScale(scaleA, inputFile);
    }
    else if (input == "scaleD") {
        setScale(scaleD, inputFile);
    }
    else if (input == "scaleI") {
        setScale(scaleI, inputFile);
    }
    else if (input == "scaleE") {
        setScale(scaleE, inputFile);
    }
    else if (input == "shapeEffectSizes") {
        inputFile >> shapeEffectSizes;
    }
    else if (input == "scaleEffectSizes") {
        inputFile >> scaleEffectSizes;
    }
    else if (input == "shapeInteractionWeights") {
        inputFile >> shapeInteractionWeights;
    }
    else if (input == "scaleInteractionWeights") {
        inputFile >> scaleInteractionWeights;
    }
    else if (input == "tBurnIn") {
        inputFile >> tBurnIn;
    }
    else if (input == "tEndSim") {
        inputFile >> tEndSim;
    }
    else if (input == "tGetDat") {
        inputFile >> tGetDat;
    }
    else if (input == "tSavDat") {
        inputFile >> tSavDat;
    }
    else if (input == "tiny") {
        inputFile >> tiny;
    }
    else {
        throw std::runtime_error("Unknown parameter " + input);
    }
}

// Function to read parameters from a parameter file
void ParameterSet::readParameters(const std::string& filename)
{
    std::clog << "Reading parameters from file " << filename << '\n';

    // Open parameter file
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        throw std::runtime_error("Unable to open parameter file in readParameters()");
    }

    // Read input
    std::string input;
    readSeed(inputFile, input);
    readIsArchitecture(inputFile, input);

    while (inputFile >> input) {
        readInput(input, inputFile);
    }

    std::clog << "Parameters were read in successfully\n";
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