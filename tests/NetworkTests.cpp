#include "library/Network.h"
#include "library/Utilities.h"
#include "library/Param.h"
#include "library/GenArch.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_CASE(checkNoNetworks)
{
    std::cout << "Testing empty networks...\n";
    Param pars;
    pars.setNLociPerTrait({ 100u, 100u, 100u });
    pars.setNEdgesPerTrait({ 0u, 0u, 0u });
    GenArch arch = GenArch(pars);
    for (size_t trait = 0u; trait < 3u; ++trait)
        BOOST_CHECK_EQUAL(arch.getNetwork(trait).map.size(), 0u);

}

BOOST_AUTO_TEST_CASE(checkSmallNetworks)
{
    std::cout << "Testing small networks...\n";
    Param pars;
    pars.setNLociPerTrait({ 100u, 100u, 100u });
    pars.setNEdgesPerTrait({ 1u, 0u, 0u });
    GenArch arch = GenArch(pars);
    Network net = arch.getNetwork(0u);
    BOOST_CHECK_EQUAL(net.map.size(), 1u);
    BOOST_CHECK_EQUAL(net.map[0u].first, 0u);
    BOOST_CHECK_EQUAL(net.map[0u].second, 1u);

}

BOOST_AUTO_TEST_CASE(checkTooBigNetworks)
{
    std::cout << "Testing that network size is capped...\n";
    Param pars;
    pars.setNLociPerTrait({ 10u, 5u, 2u });
    pars.setNEdgesPerTrait({ 1000u, 1000u, 1000u });
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.getNetwork(0u).map.size(), 45u);
    BOOST_CHECK_EQUAL(arch.getNetwork(1u).map.size(), 10u);
    BOOST_CHECK_EQUAL(arch.getNetwork(2u).map.size(), 1u);
}

BOOST_AUTO_TEST_CASE(checkInteractionWeights)
{
    std::cout << "Testing interaction weights...\n";
    Param pars;
    pars.setInteractionWeightScale(0.0);
    GenArch arch = GenArch(pars);
    for (size_t trait = 0u; trait < 3u; ++trait) {
        Network net = arch.getNetwork(trait);
        BOOST_CHECK_EQUAL(utl::sum(net.weights), 0.0);
    }
}
