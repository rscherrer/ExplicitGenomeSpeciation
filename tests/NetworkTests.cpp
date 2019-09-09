#include "library/Network.h"
#include "library/utils.h"
#include "tests/GenFixture.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_CASE(checkNoNetworks)
{
    std::cout << "Testing empty networks...\n";
    ParameterSet pars;
    pars.setNLociPerTrait({ 100u, 100u, 100u });
    pars.setNEdgesPerTrait({ 0u, 0u, 0u });
    GeneticArchitecture arch = GeneticArchitecture(pars);
    MultiNet networks = arch.getNetworks();
    for (auto net : networks)
        BOOST_CHECK_EQUAL(net.map.size(), 0u);

}

BOOST_AUTO_TEST_CASE(checkSmallNetworks)
{
    std::cout << "Testing small networks...\n";
    ParameterSet pars;
    pars.setNLociPerTrait({ 100u, 100u, 100u });
    pars.setNEdgesPerTrait({ 1u, 0u, 0u });
    GeneticArchitecture arch = GeneticArchitecture(pars);
    MultiNet networks = arch.getNetworks();
    BOOST_CHECK_EQUAL(networks[0u].map.size(), 1u);
    BOOST_CHECK_EQUAL(networks[0u].map[0u].first, 0u);
    BOOST_CHECK_EQUAL(networks[0u].map[0u].second, 1u);

}

BOOST_AUTO_TEST_CASE(checkTooBigNetworks)
{
    std::cout << "Testing that network size is capped...\n";
    ParameterSet pars;
    pars.setNLociPerTrait({ 10u, 5u, 2u });
    pars.setNEdgesPerTrait({ 1000u, 1000u, 1000u });
    GeneticArchitecture arch = GeneticArchitecture(pars);
    MultiNet networks = arch.getNetworks();
    BOOST_CHECK_EQUAL(networks[0u].map.size(), 45u);
    BOOST_CHECK_EQUAL(networks[1u].map.size(), 10u);
    BOOST_CHECK_EQUAL(networks[2u].map.size(), 1u);
}

BOOST_AUTO_TEST_CASE(checkInteractionWeights)
{
    std::cout << "Testing interaction weights...\n";
    ParameterSet pars;
    pars.setInteractionWeightScale(0.0);
    GeneticArchitecture arch = GeneticArchitecture(pars);
    MultiNet networks = arch.getNetworks();
    BOOST_CHECK_EQUAL(sum(networks[0u].weights), 0.0);
    BOOST_CHECK_EQUAL(sum(networks[1u].weights), 0.0);
    BOOST_CHECK_EQUAL(sum(networks[2u].weights), 0.0);
}
