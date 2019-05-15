#include "OutputFile.h"
#include "ParameterSet.h"
#include <sstream>

void OutputFile::open(const size_t &seed, const std::string &extension)
{
    std::ostringstream oss;
    oss << "simulation_" << seed;
    file.open(oss.str() + extension);
}



template <class T>
void OutputFile::writeLine(std::string &name, T &value)
{
    file << name << ' ' << value << '\n';
}

void OutputFile::writeLineVector(std::string &name, std::vector<double> &values)
{
    file << name << ' ';
    for (double element : values) {
        file << element << ' ';
    }
    file << '\n';
}

void OutputFile::writeParameters(const ParameterSet &parameters)
{

    writeLine("seed", parameters.seed);
    writeLine("isGenerateArchitecture", parameters.isGenerateArchitecture);
    if (parameters.isGenerateArchitecture) {
        writeLine("architectureFileName", parameters.architectureFileName);
    }
    writeLine("initialPopSize", parameters.initialPopSize);
    writeLine("dispersalRate", parameters.dispersalRate);
    writeLine("birthRate", parameters.birthRate);
    writeLine("habitatAsymmetry", parameters.habitatAsymmetry);
    writeLine("survivalProb", parameters.survivalProb);
    writeLine("ecoSelCoeff", parameters.ecoSelCoeff);
    writeLine("matePreferenceStrength", parameters.matePreferenceStrength);
    writeLine("mateEvaluationCost", parameters.mateEvaluationCost);
    writeLine("maxResourceCapacity", parameters.maxResourceCapacity);
    writeLine("maxResourceGrowth", parameters.maxResourceGrowth);
    writeLine("nEcoloci", parameters.nEcoLoci);
    writeLine("nMatLoci", parameters.nMatLoci);
    writeLine("nNtrLoci", parameters.nNtrLoci);
    writeLine("nEcoInteractions", parameters.nEcoInteractions);
    writeLine("nMatInteractions", parameters.nMatInteractions);
    writeLine("nNtrInteractions", parameters.nNtrInteractions);
    writeLine("nChromosomes", , parameters.nChromosomes);
    writeLine("initialSequence", parameters.strInitialSequence);
    writeLine("freqSNP", parameters.freqSNP);
    writeLine("mutationRate", parameters.mutationRate);
    writeLine("genomeLength", parameters.genomeLength);
    writeLine("isFemaleHeterogamy", parameters.isFemaleHeterogamy);
    writeLine("recombinationRate", parameters.recombinationRate);
    writeLine("networkSkewness", parameters.networkSkewness);
    writeLineVector("scaleA", parameters.scaleA);
    writeLineVector("scaleD", parameters.scaleD);
    writeLineVector("scaleI", parameters.scaleI);
    writeLineVector("scaleE", parameters.scaleE);
    writeLine("shapeEffectSizes", parameters.shapeEffectSizes);
    writeLine("scaleEffectSizes", parameters.scaleEffectSizes);
    writeLine("shapeInteractionWeights", parameters.shapeInteractionWeights);
    writeLine("scaleInteractionWeights", parameters.scaleInteractionWeights);
    writeLine("tBurnIn", parameters.tBurnIn);
    writeLine("tEndSim", parameters.tEndSim);
    writeLine("tGetDat", parameters.tGetDat);
    writeLine("tSavDat", parameters.tSavDat);
    writeLine("tiny", parameters.tiny);
}

void addColumn(std::ofstream &file, std::string &name)
{
    file << '\t' << name;
}

void OutputFile::writeHeader()
{
    addColumn(file, "popSize");
    addColumn(file, "nFemales");
    addColumn(file, "nMales");
    addColumn(file, "popSizeE0H0");
    addColumn(file, "popSizeE1H0");
    addColumn(file, "popSizeE0H1");
    addColumn(file, "popSizeE1H1");
    addColumn(file, "meanAttackRateE1H0");
    addColumn(file, "meanAttackRateE2H0");
    addColumn(file, "meanAttackRateE1H1");
    addColumn(file, "meanAttackRateE2H1");
    addColumn(file, "resource0H0");
    addColumn(file, "resource1H0");
    addColumn(file, "resource0H1");
    addColumn(file, "resource1H1");
    addColumn(file, "resource1H0");

    addColumn(file, "meanEcoTraitE0H0");
    addColumn(file, "meanEcoTraitE1H0");
    addColumn(file, "meanEcoTraitE0H1");
    addColumn(file, "meanEcoTraitE1H1");
    addColumn(file, "varPEcoTraitE0H0");
    addColumn(file, "varPEcoTraitE1H0");
    addColumn(file, "varPEcoTraitE0H1");
    addColumn(file, "varPEcoTraitE1H1");
    addColumn(file, "varGEcoTraitE0H0");
    addColumn(file, "varGEcoTraitE1H0");
    addColumn(file, "varGEcoTraitE0H1");
    addColumn(file, "varGEcoTraitE1H1");
    addColumn(file, "varDEcoTraitE0H0");
    addColumn(file, "varDEcoTraitE1H0");
    addColumn(file, "varDEcoTraitE0H1");
    addColumn(file, "varDEcoTraitE1H1");
    addColumn(file, "varIEcoTraitE0H0");
    addColumn(file, "varIEcoTraitE1H0");
    addColumn(file, "varIEcoTraitE0H1");
    addColumn(file, "varIEcoTraitE1H1");
    addColumn(file, "FstEcoTraitE1H1");
    addColumn(file, "GstEcoTraitE1H1");
    addColumn(file, "QstEcoTraitE1H1");
    addColumn(file, "CstEcoTraitE1H1");

    addColumn(file, "meanMatTraitE0H0");
    addColumn(file, "meanMatTraitE1H0");
    addColumn(file, "meanMatTraitE0H1");
    addColumn(file, "meanMatTraitE1H1");
    addColumn(file, "varPMatTraitE0H0");
    addColumn(file, "varPMatTraitE1H0");
    addColumn(file, "varPMatTraitE0H1");
    addColumn(file, "varPMatTraitE1H1");
    addColumn(file, "varGMatTraitE0H0");
    addColumn(file, "varGMatTraitE1H0");
    addColumn(file, "varGMatTraitE0H1");
    addColumn(file, "varGMatTraitE1H1");
    addColumn(file, "varDMatTraitE0H0");
    addColumn(file, "varDMatTraitE1H0");
    addColumn(file, "varDMatTraitE0H1");
    addColumn(file, "varDMatTraitE1H1");
    addColumn(file, "varIMatTraitE0H0");
    addColumn(file, "varIMatTraitE1H0");
    addColumn(file, "varIMatTraitE0H1");
    addColumn(file, "varIMatTraitE1H1");
    addColumn(file, "FstMatTraitE1H1");
    addColumn(file, "GstMatTraitE1H1");
    addColumn(file, "QstMatTraitE1H1");
    addColumn(file, "CstMatTraitE1H1");

    addColumn(file, "meanNtrTraitE0H0");
    addColumn(file, "meanNtrTraitE1H0");
    addColumn(file, "meanNtrTraitE0H1");
    addColumn(file, "meanNtrTraitE1H1");
    addColumn(file, "varPNtrTraitE0H0");
    addColumn(file, "varPNtrTraitE1H0");
    addColumn(file, "varPNtrTraitE0H1");
    addColumn(file, "varPNtrTraitE1H1");
    addColumn(file, "varGNtrTraitE0H0");
    addColumn(file, "varGNtrTraitE1H0");
    addColumn(file, "varGNtrTraitE0H1");
    addColumn(file, "varGNtrTraitE1H1");
    addColumn(file, "varDNtrTraitE0H0");
    addColumn(file, "varDNtrTraitE1H0");
    addColumn(file, "varDNtrTraitE0H1");
    addColumn(file, "varDNtrTraitE1H1");
    addColumn(file, "varINtrTraitE0H0");
    addColumn(file, "varINtrTraitE1H0");
    addColumn(file, "varINtrTraitE0H1");
    addColumn(file, "varINtrTraitE1H1");
    addColumn(file, "FstNtrTraitE1H1");
    addColumn(file, "GstNtrTraitE1H1");
    addColumn(file, "QstNtrTraitE1H1");
    addColumn(file, "CstNtrTraitE1H1");

    addColumn(file, "spatialIsolation");
    addColumn(file, "ecologicalIsolation");
    addColumn(file, "matingIsolation");

    file << '\n';
}