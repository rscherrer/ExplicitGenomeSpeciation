#include "OutputFile.h"
#include "ParameterSet.h"
#include <sstream>
#include <string>
#include <vector>

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

void OutputFile::addColumn(std::string &name)
{
    file << '\t' << name;
}

void OutputFile::writeHeader()
{
    addColumn("popSize");
    addColumn("nFemales");
    addColumn("nMales");
    addColumn("popSizeE0H0");
    addColumn("popSizeE1H0");
    addColumn("popSizeE0H1");
    addColumn("popSizeE1H1");
    addColumn("meanAttackRateE1H0");
    addColumn("meanAttackRateE2H0");
    addColumn("meanAttackRateE1H1");
    addColumn("meanAttackRateE2H1");
    addColumn("resource0H0");
    addColumn("resource1H0");
    addColumn("resource0H1");
    addColumn("resource1H1");
    addColumn("resource1H0");

    addColumn("meanEcoTraitE0H0");
    addColumn("meanEcoTraitE1H0");
    addColumn("meanEcoTraitE0H1");
    addColumn("meanEcoTraitE1H1");
    addColumn("varPEcoTraitE0H0");
    addColumn("varPEcoTraitE1H0");
    addColumn("varPEcoTraitE0H1");
    addColumn("varPEcoTraitE1H1");
    addColumn("varGEcoTraitE0H0");
    addColumn("varGEcoTraitE1H0");
    addColumn("varGEcoTraitE0H1");
    addColumn("varGEcoTraitE1H1");
    addColumn("varDEcoTraitE0H0");
    addColumn("varDEcoTraitE1H0");
    addColumn("varDEcoTraitE0H1");
    addColumn("varDEcoTraitE1H1");
    addColumn("varIEcoTraitE0H0");
    addColumn("varIEcoTraitE1H0");
    addColumn("varIEcoTraitE0H1");
    addColumn("varIEcoTraitE1H1");
    addColumn("FstEcoTraitE1H1");
    addColumn("GstEcoTraitE1H1");
    addColumn("QstEcoTraitE1H1");
    addColumn("CstEcoTraitE1H1");

    addColumn("meanMatTraitE0H0");
    addColumn("meanMatTraitE1H0");
    addColumn("meanMatTraitE0H1");
    addColumn("meanMatTraitE1H1");
    addColumn("varPMatTraitE0H0");
    addColumn("varPMatTraitE1H0");
    addColumn("varPMatTraitE0H1");
    addColumn("varPMatTraitE1H1");
    addColumn("varGMatTraitE0H0");
    addColumn("varGMatTraitE1H0");
    addColumn("varGMatTraitE0H1");
    addColumn("varGMatTraitE1H1");
    addColumn("varDMatTraitE0H0");
    addColumn("varDMatTraitE1H0");
    addColumn("varDMatTraitE0H1");
    addColumn("varDMatTraitE1H1");
    addColumn("varIMatTraitE0H0");
    addColumn("varIMatTraitE1H0");
    addColumn("varIMatTraitE0H1");
    addColumn("varIMatTraitE1H1");
    addColumn("FstMatTraitE1H1");
    addColumn("GstMatTraitE1H1");
    addColumn("QstMatTraitE1H1");
    addColumn("CstMatTraitE1H1");

    addColumn("meanNtrTraitE0H0");
    addColumn("meanNtrTraitE1H0");
    addColumn("meanNtrTraitE0H1");
    addColumn("meanNtrTraitE1H1");
    addColumn("varPNtrTraitE0H0");
    addColumn("varPNtrTraitE1H0");
    addColumn("varPNtrTraitE0H1");
    addColumn("varPNtrTraitE1H1");
    addColumn("varGNtrTraitE0H0");
    addColumn("varGNtrTraitE1H0");
    addColumn("varGNtrTraitE0H1");
    addColumn("varGNtrTraitE1H1");
    addColumn("varDNtrTraitE0H0");
    addColumn("varDNtrTraitE1H0");
    addColumn("varDNtrTraitE0H1");
    addColumn("varDNtrTraitE1H1");
    addColumn("varINtrTraitE0H0");
    addColumn("varINtrTraitE1H0");
    addColumn("varINtrTraitE0H1");
    addColumn("varINtrTraitE1H1");
    addColumn("FstNtrTraitE1H1");
    addColumn("GstNtrTraitE1H1");
    addColumn("QstNtrTraitE1H1");
    addColumn("CstNtrTraitE1H1");

    addColumn("spatialIsolation");
    addColumn("ecologicalIsolation");
    addColumn("matingIsolation");

    file << '\n';
}