#include "library/Collector.h"
#include "library/Random.h"
#include "library/Utilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cassert>

// Test the analysis module

BOOST_AUTO_TEST_CASE(DivergenceIsFullIfWithinGroupVarianceIsZero)
{
    std::clog << "Testing full variance partitioning...\n";
    BOOST_CHECK_EQUAL(Xst({ 0.0, 0.0, 4.0 }, { 10u, 10u, 20u }), 1.0);
    BOOST_CHECK_EQUAL(Xst({ 0.0, 0.0, 2.0 }, { 10u, 10u, 20u }), 1.0);
    BOOST_CHECK_EQUAL(Xst({ 0.0, 0.0, 1.0 }, { 10u, 10u, 20u }), 1.0);
}

BOOST_AUTO_TEST_CASE(DivergenceIsZeroIfTotalVarianceIsZero)
{
    std::clog << "Testing zero variance partitioning...\n";
    BOOST_CHECK_EQUAL(Xst({ 0.0, 0.0, 0.0 }, { 10u, 10u, 20u }), 0.0);
}

// Test that monomorphic ecotypes indeed have zero variance
BOOST_AUTO_TEST_CASE(EcologicalIsolationIsOneIfEcotypesAreMonomorphic)
{
    std::clog << "Testing complete ecological differentiation...\n";
    Param pars;
    pars.hsymmetry = 0.0; // habitats are asymmetric in resources
    pars.demesizes = { 100u, 100u };
    pars.survival = 1.0;
    pars.birth = 0.0;
    pars.tburnin = 0u;
    pars.dispersal = 0.0;
    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    metapop.resetEcoTraits(0u, -1.0, pars); // only trait -1 in habitat 0
    metapop.resetEcoTraits(1u, 1.0, pars); // only trait 1 in habitat 1
    metapop.cycle(pars, arch);
    Collector collector = Collector(arch);
    collector.analyze(metapop, pars, arch);
    BOOST_CHECK_EQUAL(collector.getEI(), 1.0);    
}

// Test case: a population with spatial isolation = 1
BOOST_AUTO_TEST_CASE(SpatialIsolationIsOneIfEcotypesAreSeparated)
{
    std::clog << "Testing complete spatial differentiation...\n";
    Param pars;
    pars.hsymmetry = 0.0; // habitats are asymmetric in resources
    pars.demesizes = { 100u, 100u };
    pars.survival = 1.0;
    pars.birth = 0.0;
    pars.tburnin = 0u;
    pars.dispersal = 0.0;
    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    metapop.resetEcoTraits(0u, -1.0, pars); // only trait -1 in habitat 0
    metapop.resetEcoTraits(0u, 1.0, pars); // only trait 1 in habitat 1
    metapop.cycle(pars, arch);
    Collector collector = Collector(arch);
    collector.analyze(metapop, pars, arch);
    BOOST_CHECK_EQUAL(collector.getSI(), 1.0);
}

// Test case: a population with mating isolation = 1
BOOST_AUTO_TEST_CASE(MatingIsolationIsOneIfMatingIsAssortative)
{
    std::clog << "Testing complete mating differentiation...\n";
    Param pars;
    pars.hsymmetry = 0.0; // habitats are asymmetric in resources
    pars.demesizes = { 100u, 100u };
    pars.survival = 1.0;
    pars.birth = 0.0;
    pars.tburnin = 0u;
    pars.dispersal = 0.0;
    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    metapop.resetEcoTraits(0u, -1.0, pars); // only trait -1 in habitat 0
    metapop.resetEcoTraits(1u, 1.0, pars); // only trait 1 in habitat 1
    metapop.resetMatePrefs(1.0); // assortative mating
    metapop.cycle(pars, arch);
    Collector collector = Collector(arch);
    collector.analyze(metapop, pars, arch);
    BOOST_CHECK_EQUAL(utl::round(collector.getRI(), 4u), 1.0);
}

BOOST_AUTO_TEST_CASE(SpatialIsolationIsZeroIfOneHabitatIsEmpty)
{
    std::clog << "Testing zero spatial differentiation...\n";
    Param pars;
    pars.demesizes = { 100u, 0u };
    pars.hsymmetry = 1.0;
    pars.dispersal = 0.0;
    pars.survival = 1.0;
    pars.birth = 0.0;
    pars.tburnin = 0u;
    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    metapop.cycle(pars, arch);
    Collector collector = Collector(arch);
    collector.analyze(metapop, pars, arch);
    BOOST_CHECK_EQUAL(collector.getSI(), 0.0);
}

BOOST_AUTO_TEST_CASE(MatingIsolationIsZeroIfOnlyOneSex)
{
    std::clog << "Testing zero mating differentiation...\n";
    Param pars;
    pars.survival = 1.0;
    pars.birth = 0.0;
    pars.tburnin = 0u;
    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    metapop.resetGenders(true); // only females
    metapop.cycle(pars, arch);
    Collector collector = Collector(arch);
    collector.analyze(metapop, pars, arch);
    BOOST_CHECK_EQUAL(collector.getRI(), 0.0);
}

BOOST_AUTO_TEST_CASE(AllIsolationMetricsAreZeroIfOnlyOneEcotype)
{
    std::clog << "Testing zero differentiation whatsoever...\n";
    Param pars;
    pars.demesizes = { 100u, 0u };
    pars.hsymmetry = 0.0;
    pars.dispersal = 0.0;
    pars.survival = 1.0;
    pars.birth = 0.0;
    pars.tburnin = 0u;
    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    metapop.resetEcoTraits(-1.0, pars); // should produce only ecotype 0
    metapop.cycle(pars, arch);
    Collector collector = Collector(arch);
    collector.analyze(metapop, pars, arch);
    BOOST_CHECK_EQUAL(collector.getEI(), 0.0);
    BOOST_CHECK_EQUAL(collector.getSI(), 0.0);
    BOOST_CHECK_EQUAL(collector.getRI(), 0.0);
}


