#include "library/GenArch.h"
#include "library/Random.h"
#include "library/Utilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cassert>

BOOST_AUTO_TEST_CASE(checkChromosomes)
{
    std::cout << "Testing chromosome lengths...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(arch.chromosomes[0u], 1.0 / 3.0);
    BOOST_CHECK_EQUAL(arch.chromosomes[1u], 2.0 / 3.0);
    BOOST_CHECK_EQUAL(arch.chromosomes[2u], 1.0);
}

BOOST_AUTO_TEST_CASE(checkEncodedTraits)
{
    std::cout << "Testing genes underlying traits...\n";
    Param pars;
    pars.setNLociPerTrait({10u, 2u, 2u});
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(utl::sumu(arch.traits), 6u);
}

BOOST_AUTO_TEST_CASE(checkEffectSizes)
{
    std::cout << "Testing gene effect sizes...\n";
    Param pars;
    pars.setEffectSizeScale(0.0);
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(utl::sum(arch.effects), 0.0);
}

BOOST_AUTO_TEST_CASE(checkDominances)
{
    std::cout << "Testing gene dominance coefficients...\n";
    Param pars;
    pars.setDominanceVariance(0.0);
    GenArch arch = GenArch(pars);
    BOOST_CHECK_EQUAL(utl::sum(arch.dominances), 0.0);
}
