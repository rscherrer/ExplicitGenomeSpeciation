#include "library/GenArch.h"
#include "library/Random.h"
#include "library/Utilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cassert>

// Test the correct generation of a genetic architecture

BOOST_AUTO_TEST_CASE(ChromosomesHaveEqualLength)
{
    Param pars;
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.chromosomes[0u], 1.0 / 3.0);
    BOOST_CHECK_EQUAL(arch.chromosomes[1u], 2.0 / 3.0);
    BOOST_CHECK_EQUAL(arch.chromosomes[2u], 1.0);
}

BOOST_AUTO_TEST_CASE(LociEncodeTheRightTraits)
{
    Param pars;
    pars.nvertices = { 10u, 2u, 2u };
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(utl::sumu(arch.traits), 6u);
}

BOOST_AUTO_TEST_CASE(EffectSizesAreZeroIfScaleParamIsZero)
{
    Param pars;
    pars.effectscale = 0.0;
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(utl::sum(arch.effects), 0.0);
}

BOOST_AUTO_TEST_CASE(DominancesAreZeroIfVarianceIsZero)
{
    Param pars;
    pars.dominancevar = 0.0;
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(utl::sum(arch.dominances), 0.0);
}

BOOST_AUTO_TEST_CASE(NetworksAreEmptyIfNoEdges)
{
    Param pars;
    pars.nvertices = utl::repUns(100u, 3u);
    pars.nloci = 300u;
    pars.nedges = utl::uzeros(3u);
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(0u), 0u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(1u), 0u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(2u), 0u);

}

BOOST_AUTO_TEST_CASE(NetworkWithOneEdgeConnectsNodesZeroAndOne)
{
    Param pars;
    pars.nvertices = utl::repUns(100u, 3u);
    pars.nloci = 300u;
    pars.nedges = { 1u, 0u, 0u };
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(0u), 1u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(1u), 0u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(2u), 0u);
    BOOST_CHECK_EQUAL(arch.getEdge(0u, 0u).first, 0u);
    BOOST_CHECK_EQUAL(arch.getEdge(0u, 0u).second, 1u);
}

BOOST_AUTO_TEST_CASE(NetworkWithTooManyEdgesIsCappedAtMaxPossible)
{
    Param pars;
    pars.nvertices = { 10u, 5u, 2u };
    pars.nloci = 17u;
    pars.nedges = utl::repUns(1000u, 3u); // very large number of edges
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(0u), 45u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(1u), 10u);
    BOOST_CHECK_EQUAL(arch.getNetworkSize(2u), 1u);
}

BOOST_AUTO_TEST_CASE(InteractionWeightsAreZeroIfScaleParamIsZero)
{
    Param pars;
    pars.interactionscale = 0.0;
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getSumWeights(0u), 0.0);
    BOOST_CHECK_EQUAL(arch.getSumWeights(1u), 0.0);
    BOOST_CHECK_EQUAL(arch.getSumWeights(2u), 0.0);
}
