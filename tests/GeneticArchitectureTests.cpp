#include "library/GeneticArchitecture.h"
#include "library/Random.h"
#include "library/utils.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cassert>

BOOST_AUTO_TEST_CASE(checkChromosomes)
{
    std::cout << "Testing chromosome lengths...\n";
    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    BOOST_CHECK_EQUAL(arch.chromosomes[0u], 1.0 / 3.0);
    BOOST_CHECK_EQUAL(arch.chromosomes[1u], 2.0 / 3.0);
    BOOST_CHECK_EQUAL(arch.chromosomes[2u], 1.0);
}

BOOST_AUTO_TEST_CASE(checkEncodedTraits)
{
    std::cout << "Testing genes underlying traits...\n";
    ParameterSet pars;
    pars.setNLociPerTrait({10u, 2u, 2u});
    GeneticArchitecture arch = GeneticArchitecture(pars);
    BOOST_CHECK_EQUAL(sumu(arch.traits), 6u);
}

BOOST_AUTO_TEST_CASE(checkEffectSizes)
{
    std::cout << "Testing gene effect sizes...\n";
    ParameterSet pars;
    pars.setEffectSizeScale(0.0);
    GeneticArchitecture arch = GeneticArchitecture(pars);
    BOOST_CHECK_EQUAL(sum(arch.effects), 0.0);
}

BOOST_AUTO_TEST_CASE(checkDominances)
{
    std::cout << "Testing gene dominance coefficients...\n";
    ParameterSet pars;
    pars.setDominanceVariance(0.0);
    GeneticArchitecture arch = GeneticArchitecture(pars);
    BOOST_CHECK_EQUAL(sum(arch.dominances), 0.0);
}
