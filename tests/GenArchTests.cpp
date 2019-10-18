#include "library/GenArch.h"
#include "library/Random.h"
#include "library/Utilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cassert>

// Test the correct generation of a genetic architecture

BOOST_AUTO_TEST_CASE(ChromosomesHaveEqualLength)
{
    std::clog << "Testing chromosome lengths...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.chromosomes[0u], 1.0 / 3.0);
    BOOST_CHECK_EQUAL(arch.chromosomes[1u], 2.0 / 3.0);
    BOOST_CHECK_EQUAL(arch.chromosomes[2u], 1.0);
}

BOOST_AUTO_TEST_CASE(LociEncodeTheRightTraits)
{
    std::clog << "Testing genetic encoding of traits...\n";
    Param pars;
    pars.nvertices = { 10u, 2u, 2u };
    pars.update();
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getSumTraits(), 6u);
}

BOOST_AUTO_TEST_CASE(EffectSizesAreZeroIfScaleParamIsZero)
{
    std::clog << "Testing locus effect sizes...\n";
    Param pars;
    pars.effectscale = 0.0;
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getSumEffects(), 0.0);
}

BOOST_AUTO_TEST_CASE(DominancesAreZeroIfVarianceIsZero)
{
    std::clog << "Testing dominance coefficients...\n";
    Param pars;
    pars.dominancevar = 0.0;
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getSumDominances(), 0.0);
}

BOOST_AUTO_TEST_CASE(NetworksAreEmptyIfNoEdges)
{
    std::clog << "Testing empty networks...\n";
    Param pars;
    pars.nvertices = utl::repUns(100u, 3u);
    pars.nedges = utl::uzeros(3u);
    pars.update();
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(0u), 0u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(1u), 0u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(2u), 0u);

}

BOOST_AUTO_TEST_CASE(NetworkWithOneEdgeConnectsNodesZeroAndOne)
{
    std::clog << "Testing single edge network...\n";
    Param pars;
    pars.nvertices = utl::repUns(100u, 3u);
    pars.nedges = { 1u, 0u, 0u };
    pars.update();
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(0u), 1u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(1u), 0u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(2u), 0u);
    BOOST_CHECK_EQUAL(arch.getEdge(0u, 0u).first, 0u);
    BOOST_CHECK_EQUAL(arch.getEdge(0u, 0u).second, 1u);
}

BOOST_AUTO_TEST_CASE(NetworkWithTooManyEdgesIsCappedAtMaxPossible)
{
    std::clog << "Testing too large networks...\n";
    Param pars;
    pars.nvertices = { 10u, 5u, 2u };
    pars.nedges = utl::repUns(1000u, 3u); // very large number of edges
    pars.update();
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(0u), 45u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(1u), 10u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(2u), 1u);
}

BOOST_AUTO_TEST_CASE(InteractionWeightsAreZeroIfScaleParamIsZero)
{
    std::clog << "Testing interaction weights...\n";
    Param pars;
    pars.interactionscale = 0.0;
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getSumWeights(0u), 0.0);
    BOOST_CHECK_EQUAL(arch.getSumWeights(1u), 0.0);
    BOOST_CHECK_EQUAL(arch.getSumWeights(2u), 0.0);
}

// Two simulations with the same architecture seed should have the exact
// same genetic architecture

BOOST_AUTO_TEST_CASE(SeedGeneticArchitecture)
{
    std::clog << "Testing that genetic architecture can be reused...\n";
    Param pars;
    rnd::rng.seed(pars.seed);
    pars.archseed = 42u;
    GenArch arch1 = GenArch(pars);
    GenArch arch2 = GenArch(pars);
    BOOST_CHECK_EQUAL(arch1.getSumEffects(), arch2.getSumEffects());
    BOOST_CHECK_EQUAL(arch1.getSumDominances(), arch2.getSumDominances());

    // Control that different seeds give different architectures
    pars.archseed = 42u;
    GenArch arch3 = GenArch(pars);
    pars.archseed = 43u;
    GenArch arch4 = GenArch(pars);
    BOOST_CHECK(arch3.getSumEffects() != arch4.getSumEffects());
    BOOST_CHECK(arch3.getSumDominances() != arch4.getSumDominances());
}
