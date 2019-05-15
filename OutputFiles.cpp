//
// Created by p278834 on 15-5-2019.
//

#include "OutputFiles.h"
#include "ParameterSet.h"
#include <sstream>

void OutputFiles::openAll(const size_t &seed)
{
    std::ostringstream oss;
    oss << "simulation_" << seed;
    logFile.open(oss.str() + ".log");
    datFile.open(oss.str() + ".dat");
    arcFile.open(oss.str() + "_fossil_record.txt");
    if (!(logFile.is_open() && datFile.is_open() && arcFile.is_open())) {
        throw std::runtime_error("Unable to open output files in main()");
    }
}

template <class T>
void OutputFiles::writeLogLine(std::string &name, T &value)
{
    logFile << name << ' ' << value << '\n';
}

void OutputFiles::writeLogLineVector(std::string &name, std::vector<double> &values)
{
    logFile << name << ' ';
    for (double element : values) {
        logFile << element << ' ';
    }
    logFile << '\n';
}

void OutputFiles::writeParameters(const ParameterSet &parameters)
{

    writeLogLine("seed", parameters.seed);
    writeLogLine("isGenerateArchitecture", parameters.isGenerateArchitecture);
    if (parameters.isGenerateArchitecture) {
        writeLogLine("architectureFileName", parameters.architectureFileName);
    }
    writeLogLine("initialPopSize", parameters.initialPopSize);
    writeLogLine("dispersalRate", parameters.dispersalRate);
    writeLogLine("birthRate", parameters.birthRate);
    writeLogLine("habitatAsymmetry", parameters.habitatAsymmetry);
    writeLogLine("survivalProb", parameters.survivalProb);
    writeLogLine("ecoSelCoeff", parameters.ecoSelCoeff);
    writeLogLine("matePreferenceStrength", parameters.matePreferenceStrength);
    writeLogLine("mateEvaluationCost", parameters.mateEvaluationCost);
    writeLogLine("maxResourceCapacity", parameters.maxResourceCapacity);
    writeLogLine("maxResourceGrowth", parameters.maxResourceGrowth);
    writeLogLine("nEcoloci", parameters.nEcoLoci);
    writeLogLine("nMatLoci", parameters.nMatLoci);
    writeLogLine("nNtrLoci", parameters.nNtrLoci);
    writeLogLine("nEcoInteractions", parameters.nEcoInteractions);
    writeLogLine("nMatInteractions", parameters.nMatInteractions);
    writeLogLine("nNtrInteractions", parameters.nNtrInteractions);
    writeLogLine("nChromosomes", , parameters.nChromosomes);
    writeLogLine("initialSequence", parameters.strInitialSequence);
    writeLogLine("freqSNP", parameters.freqSNP);
    writeLogLine("mutationRate", parameters.mutationRate);
    writeLogLine("genomeLength", parameters.genomeLength);
    writeLogLine("isFemaleHeterogamy", parameters.isFemaleHeterogamy);
    writeLogLine("recombinationRate", parameters.recombinationRate);
    writeLogLine("networkSkewness", parameters.networkSkewness);
    writeLogLineVector("scaleA", parameters.scaleA);
    writeLogLineVector("scaleD", parameters.scaleD);
    writeLogLineVector("scaleI", parameters.scaleI);
    writeLogLineVector("scaleE", parameters.scaleE);
    writeLogLine("shapeEffectSizes", parameters.shapeEffectSizes);
    writeLogLine("scaleEffectSizes", parameters.scaleEffectSizes);
    writeLogLine("shapeInteractionWeights", parameters.shapeInteractionWeights);
    writeLogLine("scaleInteractionWeights", parameters.scaleInteractionWeights);
    writeLogLine("tBurnIn", parameters.tBurnIn);
    writeLogLine("tEndSim", parameters.tEndSim);
    writeLogLine("tGetDat", parameters.tGetDat);
    writeLogLine("tSavDat", parameters.tSavDat);
    writeLogLine("tiny", parameters.tiny);
}
