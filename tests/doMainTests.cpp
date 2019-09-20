#define BOOST_TEST_MAIN

#include "library/doMain.h"
#include "library/Population.h"
#include "library/MetaPop.h"
#include "tests/GenFixture.h"
#include "tests/testUtilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

typedef std::vector<Network> MultiNet;

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

    const size_t n0 = pars.getInitialPopSize();
    const double s = pars.getEcoSelCoeff();
    const double max = pars.getMaxFeedingRate();
    const vecDbl k = rep(pars.getMaxResourceCapacity(), 2u);
    const vecDbl r = rep(pars.getMaxResourceGrowth(), 2u);

    Population pop1 = Population(n0, s, max, k, r, arch);
    Population pop2 = Population(n0, s, max, k, r, arch);

    MetaPop meta = MetaPop({ pop1, pop2 }, pars);

    int t = meta.evolve(arch);

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

    const size_t n0 = pars.getInitialPopSize();
    const double s = pars.getEcoSelCoeff();
    const double max = pars.getMaxFeedingRate();
    const vecDbl k = rep(pars.getMaxResourceCapacity(), 2u);
    const vecDbl r = rep(pars.getMaxResourceGrowth(), 2u);

    Population pop1 = Population(n0, s, max, k, r, arch);
    Population pop2 = Population(n0, s, max, k, r, arch);

    MetaPop meta = MetaPop({ pop1, pop2 }, pars);

    int t = meta.evolve(arch);

    BOOST_CHECK(t < pars.getTEndSim());
    BOOST_CHECK(meta.getPops()[0u].getPopSize() == 0u);
    BOOST_CHECK(meta.getPops()[1u].getPopSize() == 0u);
}

BOOST_FIXTURE_TEST_SUITE(analysisTestSuite, GenFixture)

    // Test case: a population with ecological isolation = 1
    BOOST_AUTO_TEST_CASE(fullEcologicalIsolation)
    {

        std::cout << "Testing full ecological isolation...\n";

        const size_t n0 = pars.getInitialPopSize();
        const double s = pars.getEcoSelCoeff();
        const double max = pars.getMaxFeedingRate();
        const vecDbl k = rep(pars.getMaxResourceCapacity(), 2u);
        const vecDbl r = rep(pars.getMaxResourceGrowth(), 2u);

        Population pop1 = Population(n0, s, max, k, r, arch);
        Population pop2 = Population(n0, s, max, k, r, arch);

        pop1.resetEcoTraits(-1.0, 1.0, 4.0E-4);
        pop2.resetEcoTraits(1.0, 1.0, 4.0E-4);
        MetaPop meta = MetaPop({ pop1, pop2 }, pars);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getEcoIsolation(), 1.0);

    }

    // Test case: a population with spatial isolation = 1
    BOOST_AUTO_TEST_CASE(fullSpatialIsolation)
    {

        std::cout << "Testing full spatial isolation...\n";

        const size_t n0 = pars.getInitialPopSize();
        const double s = pars.getEcoSelCoeff();
        const double max = pars.getMaxFeedingRate();
        const vecDbl k = rep(pars.getMaxResourceCapacity(), 2u);
        const vecDbl r = rep(pars.getMaxResourceGrowth(), 2u);

        Population pop1 = Population(n0, s, max, k, r, arch);
        Population pop2 = Population(n0, s, max, k, r, arch);
        pop1.resetEcoTraits(-1.0, 1.0, 4.0E-4);
        pop2.resetEcoTraits(1.0, 1.0, 4.0E-4);
        MetaPop meta = MetaPop({ pop1, pop2 }, pars);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 1.0);

    }

    // Test case: a population with mating isolation = 1
    BOOST_AUTO_TEST_CASE(fullMatingIsolation)
    {

        std::cout << "Testing full mating isolation...\n";

        const size_t n0 = pars.getInitialPopSize();
        const double s = pars.getEcoSelCoeff();
        const double max = pars.getMaxFeedingRate();
        const vecDbl k = rep(pars.getMaxResourceCapacity(), 2u);
        const vecDbl r = rep(pars.getMaxResourceGrowth(), 2u);

        Population pop1 = Population(n0, s, max, k, r, arch);
        Population pop2 = Population(n0, s, max, k, r, arch);

        pop1.resetEcoTraits(-1.0, 1.0, 4.0E-4);
        pop2.resetEcoTraits(1.0, 1.0, 4.0E-4);
        pop1.resetMatePrefs(1.0);
        pop2.resetMatePrefs(1.0);
        MetaPop meta = MetaPop({ pop1, pop2 }, pars);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 1.0);

    }

    BOOST_AUTO_TEST_CASE(abuseSpatialIsolationOnePop)
    {

        std::cout << "Testing spatial isolation with only one population...\n";

        const size_t n0 = pars.getInitialPopSize();
        const double s = pars.getEcoSelCoeff();
        const double max = pars.getMaxFeedingRate();
        const vecDbl k = rep(pars.getMaxResourceCapacity(), 2u);
        const vecDbl r = rep(pars.getMaxResourceGrowth(), 2u);

        Population pop1 = Population(n0, s, max, k, r, arch);
        Population pop2 = Population(0u, s, max, k, r, arch);
        pop1.resetEcoTraits(-1.0, 1.0, 4.0E-4);
        MetaPop meta = MetaPop({ pop1, pop2 }, pars);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 0.0);
    }

    BOOST_AUTO_TEST_CASE(abuseSpatialIsolationOneEcotype)
    {

        std::cout << "Testing spatial isolation with only one ecotype...\n";

        const size_t n0 = pars.getInitialPopSize();
        const double s = pars.getEcoSelCoeff();
        const double max = pars.getMaxFeedingRate();
        const vecDbl k = rep(pars.getMaxResourceCapacity(), 2u);
        const vecDbl r = rep(pars.getMaxResourceGrowth(), 2u);

        Population pop1 = Population(n0, s, max, k, r, arch);
        Population pop2 = Population(n0, s, max, k, r, arch);
        pop1.resetEcotypes(1u);
        pop2.resetEcotypes(1u);
        MetaPop meta = MetaPop({ pop1, pop2 }, pars);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 0.0);
    }

    BOOST_AUTO_TEST_CASE(abuseMatingIsolationOneSex)
    {

        std::cout << "Testing mating isolation when only one sex...\n";

        const size_t n0 = pars.getInitialPopSize();
        const double s = pars.getEcoSelCoeff();
        const double max = pars.getMaxFeedingRate();
        const vecDbl k = rep(pars.getMaxResourceCapacity(), 2u);
        const vecDbl r = rep(pars.getMaxResourceGrowth(), 2u);

        Population pop1 = Population(n0, s, max, k, r, arch);
        Population pop2 = Population(0u, s, max, k, r, arch);
        pop1.resetGenders(true); // only females
        MetaPop meta = MetaPop({ pop1, pop2 }, pars);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 0.0);
    }

    BOOST_AUTO_TEST_CASE(abuseMatingIsolationOneEcotype)
    {

        std::cout << "Testing mating isolation when one ecotype...\n";

        const size_t n0 = pars.getInitialPopSize();
        const double s = pars.getEcoSelCoeff();
        const double max = pars.getMaxFeedingRate();
        const vecDbl k = rep(pars.getMaxResourceCapacity(), 2u);
        const vecDbl r = rep(pars.getMaxResourceGrowth(), 2u);

        Population pop1 = Population(n0, s, max, k, r, arch);
        Population pop2 = Population(0u, s, max, k, r, arch);
        pop1.resetEcotypes(1u);
        MetaPop meta = MetaPop({ pop1, pop2 }, pars);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 0.0);
    }


BOOST_AUTO_TEST_SUITE_END()

