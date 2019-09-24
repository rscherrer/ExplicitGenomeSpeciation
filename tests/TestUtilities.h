#ifndef EXPLICITGENOMESPECIATION_TESTUTILITIES_H
#define EXPLICITGENOMESPECIATION_TESTUTILITIES_H

#include <boost/test/included/unit_test.hpp>
#include <vector>

void makeValidParamFile()
{
    std::ofstream out;
    out.open("valid_paramfile_test.txt");
    if (!out.is_open())
        std::cout << "Unable to open valid parameter test file.\n";

    out << "maxResourceCapacity" << '\t' << 100.0 << '\n'
        << "maxResourceGrowth" << '\t' << 1.0 << '\n'
        << "habitatSymmetry" << '\t' << 1.0 << '\n'
        << "ecoSelCoeff" << '\t' << 1.0 << '\n'
        << "initialPopSize" << '\t' << 100 << '\n'
        << "dispersalRate" << '\t' << 1.0e-3 << '\n'
        << "birthRate" << '\t' << 4.0 << '\n'
        << "survivalProb" << '\t' << 0.8 << '\n'
        << "matePreferenceStrength" << '\t' << 10.0 << '\n'
        << "mateEvalutationCost" << '\t' << 0.01 << '\n'
        << "nEcoLoci" << '\t' << 400 << '\n'
        << "nMatLoci" << '\t' << 200 << '\n'
        << "nNtrLoci" << '\t' << 400 << '\n'
        << "nEcoEdges" << '\t' << 1000 << '\n'
        << "nMatEdges" << '\t' << 500 << '\n'
        << "nNtrEdges" << '\t' << 0 << '\n'
        << "nChromosomes" << '\t' << 3 << '\n'
        << "mutationRate" << '\t' << 1.0e-5 << '\n'
        << "recombinationRate" << '\t' << 0.01 << '\n'
        << "freqSNP" << '\t' << 0.02 << '\n'
        << "isFemaleHeterogamy" << '\t' << 0 << '\n'
        // << "isGeneticArchitecture" << '\t' << 0 << '\n'
        // << "architectureFileName" << '\t' << "architecture.txt" << '\n'
        << "scaleA" << '\t' << 1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
        << "scaleD" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
        << "scaleI" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
        << "scaleE" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
        << "skewnesses" << '\t' << 1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
        << "effectSizeShape" << '\t' << 2.0 << '\n'
        << "effectSizeScale" << '\t' << 1.0 << '\n'
        << "interactionWeightShape" << '\t' << 5.0 << '\n'
        << "interactionWeightScale" << '\t' << 1.0 << '\n'
        << "dominanceVariance" << '\t' << 1.0 << '\n'
        << "tBurnIn" << '\t' << 1 << '\n'
        << "tEndSim" << '\t' << 5 << '\n'
        << "tSave" << '\t' << 1 << '\n'
        << "seed" << '\t' << 42 << '\n'
        << "record" << '\t' << 1 << '\n';

    out.close();
}

void makeInvalidParamName()
{
    std::ofstream out;
    out.open("invalid_paramname_test.txt");
    if (!out.is_open())
        std::cout << "Unable to open invalid parameter name test file.\n";
    out << "nonsense" << '\t' << 3.0 << '\n';
    out.close();
}

void makeInvalidParamValue()
{
    std::ofstream out;
    out.open("invalid_paramvalue_test.txt");
    if (!out.is_open())
        std::cout << "Unable to open invalid parameter value test file.\n";

    out << "maxResourceCapacity" << '\t' << -1.0 << '\n'
        << "maxResourceGrowth" << '\t' << -1.0 << '\n'
        << "habitatSymmetry" << '\t' << -1.0 << '\n'
        << "ecoSelCoeff" << '\t' << -1.0 << '\n'
        << "initialPopSize" << '\t' << -100 << '\n'
        << "dispersalRate" << '\t' << 10 << '\n'
        << "birthRate" << '\t' << -4.0 << '\n'
        << "survivalProb" << '\t' << -0.8 << '\n'
        << "matePreferenceStrength" << '\t' << -10.0 << '\n'
        << "mateEvalutationCost" << '\t' << -0.01 << '\n'
        << "nEcoLoci" << '\t' << 1 << '\n'
        << "nMatLoci" << '\t' << 1 << '\n'
        << "nNtrLoci" << '\t' << 1 << '\n'
        << "nChromosomes" << '\t' << 0 << '\n'
        << "mutationRate" << '\t' << -1.0e-5 << '\n'
        << "recombinationRate" << '\t' << -0.01 << '\n'
        << "freqSNP" << '\t' << -0.02 << '\n'
        << "scaleA" << '\t' << -1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
        << "scaleD" << '\t' << -1.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
        << "scaleI" << '\t' << -1.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
        << "scaleE" << '\t' << -1.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
        << "skewnesses" << '\t' << -1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
        << "effectSizeShape" << '\t' << -2.0 << '\n'
        << "effectSizeScale" << '\t' << -1.0 << '\n'
        << "interactionWeightShape" << '\t' << -5.0 << '\n'
        << "interactionWeightScale" << '\t' << -1.0 << '\n'
        << "dominanceVariance" << '\t' << -1.0 << '\n';

    out.close();
}

void makeInvalidParamValue2()
{
    std::ofstream out;
    out.open("invalid_paramvalue_test.txt");
    if (!out.is_open())
        std::cout << "Unable to open invalid parameter value test file.\n";

    out << "dispersalRate" << '\t' << -1.0e-3 << '\n'
        << "survivalProb" << '\t' << 2.0 << '\n'
        << "mutationRate" << '\t' << 10 << '\n'
        << "habitatSymmetry" << '\t' << 2.0 << '\n'
        << "freqSNP" << '\t' << 10 << '\n';

    out.close();
}



#endif
