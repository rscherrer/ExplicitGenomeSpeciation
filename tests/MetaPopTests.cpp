#include "library/MetaPop.h"
#include "library/Utilities.h"
#include "tests/PopFixture.h"
#include <boost/test/unit_test.hpp>

// Simulation should reach tmax in the absence of mortality
BOOST_AUTO_TEST_CASE(checkImmortalPopulation)
{

    std::cout << "Testing an immortal population...\n";

    Param pars;
    pars.setTEndSim(100u);
    pars.setTSave(1u);
    pars.setInitialPopSize(10u);
    pars.setDispersalRate(0.0);
    pars.setSurvivalProb(1.0);
    pars.setBirthRate(0.0);
    pars.setMatePreferenceStrength(0.0);

    GenArch arch = GenArch(pars);
    const size_t n0 = pars.getInitialPopSize();
    MetaPop meta = MetaPop(utl::repUns(n0, 2u), pars, arch);
    int t = meta.evolve(arch);

    BOOST_CHECK_EQUAL(t, pars.getTEndSim());
    BOOST_CHECK(meta.getPopSize(0u) > 0u);
    BOOST_CHECK(meta.getPopSize(1u) > 0u);
}


// Simulation should end prematurely with high mortality
BOOST_AUTO_TEST_CASE(checkProgressiveExtinction)
{

    std::cout << "Testing progressive extinction...\n";

    Param pars;
    pars.setTEndSim(100u);
    pars.setTSave(1u);
    pars.setInitialPopSize(10u);
    pars.setDispersalRate(0.0);
    pars.setSurvivalProb(0.1);
    pars.setBirthRate(0.0);
    pars.setMatePreferenceStrength(0.0);

    GenArch arch = GenArch(pars);
    const size_t n0 = pars.getInitialPopSize();
    MetaPop meta = MetaPop(utl::repUns(n0, 2u), pars, arch);
    int t = meta.evolve(arch);

    BOOST_CHECK(t < pars.getTEndSim());
    BOOST_CHECK(meta.getPopSize(0u) == 0u);
    BOOST_CHECK(meta.getPopSize(1u) == 0u);
}

BOOST_FIXTURE_TEST_SUITE(analysisTestSuite, PopFixture)

    // Test case: a population with ecological isolation = 1
    BOOST_AUTO_TEST_CASE(fullEcologicalIsolation)
    {

        std::cout << "Testing full ecological isolation...\n";
        MetaPop meta = MetaPop(utl::repUns(n0, 2u), pars, arch);
        meta.resetEcoTraits(0u, -1.0);
        meta.resetEcoTraits(1u, 1.0);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getEcoIsolation(), 1.0);

    }

    // Test case: a population with spatial isolation = 1
    BOOST_AUTO_TEST_CASE(fullSpatialIsolation)
    {

        std::cout << "Testing full spatial isolation...\n";
        MetaPop meta = MetaPop(utl::repUns(n0, 2u), pars, arch);
        meta.resetEcoTraits(0u, -1.0);
        meta.resetEcoTraits(1u, 1.0);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 1.0);

    }

    // Test case: a population with mating isolation = 1
    BOOST_AUTO_TEST_CASE(fullMatingIsolation)
    {

        std::cout << "Testing full mating isolation...\n";
        MetaPop meta = MetaPop(utl::repUns(n0, 2u), pars, arch);
        meta.resetEcoTraits(0u, -1.0);
        meta.resetEcoTraits(1u, 1.0);
        meta.resetMatePrefs(0u, 1.0);
        meta.resetMatePrefs(1u, 1.0);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 1.0);

    }

    BOOST_AUTO_TEST_CASE(abuseSpatialIsolationOnePop)
    {

        std::cout << "Testing spatial isolation with only one population...\n";
        MetaPop meta = MetaPop({ n0, 0u }, pars, arch);
        meta.resetEcoTraits(0u, -1.0);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 0.0);
    }

    BOOST_AUTO_TEST_CASE(abuseSpatialIsolationOneEcotype)
    {

        std::cout << "Testing spatial isolation with only one ecotype...\n";
        MetaPop meta = MetaPop(utl::repUns(n0, 2u), pars, arch);
        meta.analyze(arch);
        meta.resetEcotypes(0u, 1u);
        meta.resetEcotypes(1u, 1u);
        BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 0.0);
    }

    BOOST_AUTO_TEST_CASE(abuseMatingIsolationOneSex)
    {

        std::cout << "Testing mating isolation when only one sex...\n";
        MetaPop meta = MetaPop({ n0, 0 }, pars, arch);
        meta.resetGenders(0u, true);
        meta.analyze(arch);
        BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 0.0);
    }

    BOOST_AUTO_TEST_CASE(abuseMatingIsolationOneEcotype)
    {

        std::cout << "Testing mating isolation when one ecotype...\n";
        MetaPop meta = MetaPop({ n0, 0u }, pars, arch);
        meta.analyze(arch);
        meta.resetEcotypes(0u, 1u);
        BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 0.0);
    }


BOOST_AUTO_TEST_SUITE_END()


