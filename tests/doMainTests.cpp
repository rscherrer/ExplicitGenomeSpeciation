#define BOOST_TEST_MAIN

#include "library/doMain.h"
#include "library/Population.h"
#include "library/MetaPop.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

typedef std::vector<Network> MultiNet;

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
        << "genomeLength" << '\t' << 300.0 << '\n'
        << "freqSNP" << '\t' << 0.02 << '\n'
        << "isFemaleHeterogamy" << '\t' << 0 << '\n'
        << "isGeneticArchitecture" << '\t' << 0 << '\n'
        << "architectureFileName" << '\t' << "architecture.txt" << '\n'
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
        << "tiny" << '\t' << 1.0e-12 << '\n'
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
        << "genomeLength" << '\t' << -300.0 << '\n'
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
        << "dominanceVariance" << '\t' << -1.0 << '\n'
        << "tiny" << '\t' << -1.0e-12 << '\n';

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

// Check that the program can run without arguments
BOOST_AUTO_TEST_CASE(testUseNoArgs)
{
    std::cout << "Testing that the main runs without arguments...\n";
    BOOST_CHECK_EQUAL(doMain({ "EGS_test" }), 0);
}


// Check that the program cannot run with more than one argument
BOOST_AUTO_TEST_CASE(testAbuseTooManyArgs)
{
    std::cout << "Testing providing too many arguments to the main...\n";
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "arg1", "arg2" }), 1);
}

BOOST_AUTO_TEST_CASE(testAbuseInvalidFilename)
{
    std::cout << "Testing invalid parameter file name...\n";
    makeValidParamFile();
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "nonsense.txt" }), 1);
}


BOOST_AUTO_TEST_CASE(testUseValidFilename)
{
    std::cout << "Testing valid parameter file name...\n";
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "valid_paramfile_test.txt" }), 0);
}


BOOST_AUTO_TEST_CASE(testAbuseInvalidParamName)
{
    std::cout << "Testing invalid parameter name...\n";
    makeInvalidParamName();
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "invalid_paramname_test.txt" }), 1);
}

BOOST_AUTO_TEST_CASE(testAbuseInvalidParamValue)
{
    std::cout << "Testing invalid parameter value...\n";
    makeInvalidParamValue();
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "invalid_paramvalue_test.txt" }), 1);
    makeInvalidParamValue2();
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "invalid_paramvalue_test.txt" }), 1);
}


// Simulation should reach tmax in the absence of mortality
BOOST_AUTO_TEST_CASE(checkImmortalPopulation)
{

    std::cout << "Testing an immortal population...\n";

    ParameterSet pars;
    pars.setTEndSim(100u);
    pars.setTSave(1u);
    pars.setInitialPopSize(10u);
    pars.setDispersalRate(0.0);
    pars.setSurvivalProb(1.0);
    pars.setBirthRate(0.0);
    pars.setMatePreferenceStrength(0.0);

    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();

    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(pars.getInitialPopSize(), genome, networks);

    MetaPop meta = MetaPop({ pop1, pop2 }, pars);

    int t = meta.evolve(genome, networks);

    BOOST_CHECK_EQUAL(t, pars.getTEndSim());
    BOOST_CHECK(meta.getPops()[0u].getPopSize() > 0u);
    BOOST_CHECK(meta.getPops()[1u].getPopSize() > 0u);
}


// Simulation should end prematurely with high mortality
BOOST_AUTO_TEST_CASE(checkProgressiveExtinction)
{

    std::cout << "Testing progressive extinction...\n";

    ParameterSet pars;
    pars.setTEndSim(100u);
    pars.setTSave(1u);
    pars.setInitialPopSize(10u);
    pars.setDispersalRate(0.0);
    pars.setSurvivalProb(0.1);
    pars.setBirthRate(0.0);
    pars.setMatePreferenceStrength(0.0);

    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();

    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(pars.getInitialPopSize(), genome, networks);

    MetaPop meta = MetaPop({ pop1, pop2 }, pars);

    int t = meta.evolve(genome, networks);

    BOOST_CHECK(t < pars.getTEndSim());
    BOOST_CHECK(meta.getPops()[0u].getPopSize() == 0u);
    BOOST_CHECK(meta.getPops()[1u].getPopSize() == 0u);
}

// Test case: a population with ecological isolation = 1
BOOST_AUTO_TEST_CASE(fullEcologicalIsolation)
{

    std::cout << "Testing full ecological isolation...\n";

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();
    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(pars.getInitialPopSize(), genome, networks);
    pop1.resetEcoTraits(-1.0, 1.0);
    pop2.resetEcoTraits(1.0, 1.0);
    MetaPop meta = MetaPop({ pop1, pop2 }, pars);
    meta.analyze(genome.nloci, genome.traits, pars.getScaleE());
    BOOST_CHECK_EQUAL(meta.getEcoIsolation(), 1.0);

}

// Test case: a population with spatial isolation = 1
BOOST_AUTO_TEST_CASE(fullSpatialIsolation)
{

    std::cout << "Testing full spatial isolation...\n";

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();
    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(pars.getInitialPopSize(), genome, networks);
    pop1.resetEcoTraits(-1.0, 1.0);
    pop2.resetEcoTraits(1.0, 1.0);
    MetaPop meta = MetaPop({ pop1, pop2 }, pars);
    meta.analyze(genome.nloci, genome.traits, pars.getScaleE());
    BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 1.0);

}

// Test case: a population with mating isolation = 1
BOOST_AUTO_TEST_CASE(fullMatingIsolation)
{

    std::cout << "Testing full mating isolation...\n";

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();
    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(pars.getInitialPopSize(), genome, networks);
    pop1.resetEcoTraits(-1.0, 1.0);
    pop2.resetEcoTraits(1.0, 1.0);
    pop1.resetMatePrefs(1.0);
    pop2.resetMatePrefs(1.0);
    MetaPop meta = MetaPop({ pop1, pop2 }, pars);
    meta.analyze(genome.nloci, genome.traits, pars.getScaleE());
    BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 1.0);

}

BOOST_AUTO_TEST_CASE(abuseSpatialIsolationOnePop)
{

    std::cout << "Testing spatial isolation with only one population...\n";

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();
    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(0u, genome, networks);
    pop1.resetEcoTraits(-1.0, 1.0);
    MetaPop meta = MetaPop({ pop1, pop2 }, pars);
    meta.analyze(genome.nloci, genome.traits, pars.getScaleE());
    BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 0.0);
}

BOOST_AUTO_TEST_CASE(abuseSpatialIsolationOneEcotype)
{

    std::cout << "Testing spatial isolation with only one ecotype...\n";

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();
    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(pars.getInitialPopSize(), genome, networks);
    pop1.resetEcotypes(1u);
    pop2.resetEcotypes(1u);
    MetaPop meta = MetaPop({ pop1, pop2 }, pars);
    meta.analyze(genome.nloci, genome.traits, pars.getScaleE());
    BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 0.0);
}

BOOST_AUTO_TEST_CASE(abuseMatingIsolationOneSex)
{

    std::cout << "Testing mating isolation when only one sex...\n";

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();
    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(0u, genome, networks);
    pop1.resetGenders(true); // only females
    MetaPop meta = MetaPop({ pop1, pop2 }, pars);
    meta.analyze(genome.nloci, genome.traits, pars.getScaleE());
    BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 0.0);
}

BOOST_AUTO_TEST_CASE(abuseMatingIsolationOneEcotype)
{

    std::cout << "Testing mating isolation when one reproducing ecotype...\n";

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();
    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(0u, genome, networks);
    pop1.resetEcotypes(1u);
    MetaPop meta = MetaPop({ pop1, pop2 }, pars);
    meta.analyze(genome.nloci, genome.traits, pars.getScaleE());
    BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 0.0);
}
