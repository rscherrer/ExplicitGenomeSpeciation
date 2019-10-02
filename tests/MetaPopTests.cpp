#include "library/MetaPop.h"
#include "library/Utilities.h"
#include "tests/DemeFixture.h"
#include <boost/test/unit_test.hpp>

// Simulation should reach tmax in the absence of mortality
BOOST_AUTO_TEST_CASE(checkImmortalPopulation)
{

    std::clog << "Testing an immortal population...\n";

    Param pars;
    pars.setTEndSim(10u);
    pars.setTSave(1u);
    pars.setInitialPopSize(10u);
    pars.setDispersalRate(0.0);
    pars.setSurvivalProb(1.0);
    pars.setBirthRate(0.0);
    pars.setMatePreferenceStrength(0.0);

    GenArch arch = GenArch(pars);
    const size_t n0 = pars.getInitialPopSize();
    MetaPop meta = MetaPop(utl::repUns(n0, 2u), pars, arch, true);
    int t = meta.evolve(arch);

    BOOST_CHECK_EQUAL(t, pars.getTEndSim());
    BOOST_CHECK(meta.getPopSize(0u) > 0u);
    BOOST_CHECK(meta.getPopSize(1u) > 0u);
}


// Simulation should end prematurely with high mortality
BOOST_AUTO_TEST_CASE(checkProgressiveExtinction)
{

    std::clog << "Testing progressive extinction...\n";

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
    MetaPop meta = MetaPop(utl::repUns(n0, 2u), pars, arch, true);
    int t = meta.evolve(arch);

    BOOST_CHECK(t < pars.getTEndSim());
    BOOST_CHECK(meta.getPopSize(0u) == 0u);
    BOOST_CHECK(meta.getPopSize(1u) == 0u);
}

// Test that habitats are initialized with only one resource if
// habitat symmetry is zero
BOOST_AUTO_TEST_CASE(habitatsHaveOneResourceIfCompleteAsymmetry)
{
    std::clog << "Testing habitat asymmetry...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    pars.setHabitatSymmetry(0.0); // full habitat asymmetry
    MetaPop meta = MetaPop(utl::repUns(100u, 2u), pars, arch, false);
    BOOST_CHECK_EQUAL(meta.getResource(0u, 0u), pars.getMaxResourceCapacity());
    BOOST_CHECK_EQUAL(meta.getResource(1u, 1u), pars.getMaxResourceCapacity());
    BOOST_CHECK_EQUAL(meta.getResource(0u, 1u), 0.0);
    BOOST_CHECK_EQUAL(meta.getResource(1u, 0u), 0.0);

}

// Test that Xst should be 1 if ecotypic variances are zero
BOOST_AUTO_TEST_CASE(XstDoesItsJob)
{

    std::clog << "Checking that variance is correctly partitioned...\n";
    BOOST_CHECK_EQUAL(utl::Xst(0.0, 0.0, 42.0, 10u, 10u, 20u), 1.0);
    BOOST_CHECK_EQUAL(utl::Xst(0.0, 0.0, 2.0, 10u, 10u, 20u), 1.0);
    BOOST_CHECK_EQUAL(utl::Xst(0.0, 0.0, 1.0, 10u, 10u, 20u), 1.0);

}

// Test that monomorphic ecotypes indeed have zero variance
BOOST_AUTO_TEST_CASE(fullEcologicalIsolation)
{
    std::clog << "Checking full ecological isolation...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    pars.setHabitatSymmetry(0.0); // full habitat asymmetry
    MetaPop meta = MetaPop(utl::repUns(100u, 2u), pars, arch, false);
    meta.resetEcoTraits(0u, -1.0); // only trait -1 in habitat 0
    meta.resetEcoTraits(1u, 1.0); // only trait 1 in habitat 1
    meta.consume();
    meta.analyze(arch);

    BOOST_CHECK_EQUAL(meta.getSumTrait(0u, 0u), -1.0 * meta.getPopSize(0u));
    BOOST_CHECK_EQUAL(meta.getSumTrait(0u, 1u), meta.getPopSize(1u));
    BOOST_CHECK_EQUAL(meta.getEcoIsolation(), 1.0);
    BOOST_CHECK_EQUAL(meta.getPst(0u), 1.0);
    BOOST_CHECK_EQUAL(meta.getVarP(0u, 0u), 0.0);
    BOOST_CHECK_EQUAL(meta.getVarP(0u, 1u), 0.0);
    BOOST_CHECK_EQUAL(meta.getSsqPhe(0u, 0u), meta.getEcoCount(0u));
    BOOST_CHECK_EQUAL(meta.getSsqPhe(0u, 1u), meta.getEcoCount(1u));
    BOOST_CHECK_EQUAL(meta.getSumPhe(0u, 0u), -1.0 * meta.getEcoCount(0u));
    BOOST_CHECK_EQUAL(meta.getSumPhe(0u, 1u), meta.getEcoCount(1u));
    BOOST_CHECK_EQUAL(meta.getSumEcotypes(0u), 0u);
    BOOST_CHECK_EQUAL(meta.getSumEcotypes(1u), meta.getPopSize(1u));

}

// Test case: a population with spatial isolation = 1
BOOST_AUTO_TEST_CASE(fullSpatialIsolation)
{
    std::clog << "Testing full spatial isolation...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    pars.setHabitatSymmetry(0.0); // full habitat asymmetry
    MetaPop meta = MetaPop(utl::repUns(100u, 2u), pars, arch, false);
    meta.resetEcoTraits(0u, -1.0); // only trait -1 in habitat 0
    meta.resetEcoTraits(1u, 1.0); // only trait 1 in habitat 1
    meta.consume();
    meta.analyze(arch);
    BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 1.0);

}

// Test case: a population with mating isolation = 1
BOOST_AUTO_TEST_CASE(fullMatingIsolation)
{
    std::clog << "Testing full mating isolation...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    pars.setHabitatSymmetry(0.0);
    MetaPop meta = MetaPop(utl::repUns(100u, 2u), pars, arch, false);
    meta.resetEcoTraits(0u, -1.0); // only trait -1 in habitat 0
    meta.resetEcoTraits(1u, 1.0); // only trait 1 in habitat 1
    meta.resetMatePrefs(0u, 1.0); // assortative mating everywhere
    meta.resetMatePrefs(1u, 1.0);
    meta.consume();
    meta.sortSexes();
    meta.analyze(arch);
    BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 1.0);

}

BOOST_AUTO_TEST_CASE(abuseSpatialIsolationOnePop)
{
    std::clog << "Testing spatial isolation with only one population...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    MetaPop meta = MetaPop({ 100u, 0u }, pars, arch, false);
    meta.resetEcoTraits(0u, -1.0);
    meta.consume();
    meta.analyze(arch);
    BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 0.0);
}

BOOST_AUTO_TEST_CASE(abuseSpatialIsolationOneEcotype)
{

    std::clog << "Testing spatial isolation with only one ecotype...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    MetaPop meta = MetaPop(utl::repUns(100u, 2u), pars, arch, false);
    meta.consume();
    meta.resetEcotypes(0u, 1u);
    meta.resetEcotypes(1u, 1u);
    meta.analyze(arch);
    BOOST_CHECK_EQUAL(meta.getSpatialIsolation(), 0.0);
}

BOOST_AUTO_TEST_CASE(abuseMatingIsolationOneSex)
{

    std::clog << "Testing mating isolation when only one sex...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    MetaPop meta = MetaPop({ 100u, 0 }, pars, arch, false);
    meta.consume();
    meta.resetGenders(0u, true);
    meta.analyze(arch);
    BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 0.0);
}

BOOST_AUTO_TEST_CASE(abuseMatingIsolationOneEcotype)
{

    std::clog << "Testing mating isolation when one ecotype...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    MetaPop meta = MetaPop({ 100u, 0u }, pars, arch, false);
    meta.consume();
    meta.resetEcotypes(0u, 1u);
    meta.analyze(arch);
    BOOST_CHECK_EQUAL(meta.getMatingIsolation(), 0.0);
}


