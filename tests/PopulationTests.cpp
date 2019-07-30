#include "library/Population.h"
#include <boost/test/unit_test.hpp>
#include <iostream>


// Check that zero survival leaves no survivors in the population
BOOST_AUTO_TEST_CASE(checkNoSurvivors)
{
    std::cout << "Testing population wipe-out...\n";
    Population pop = Population(10u);
    pop.survive(0.0);
    BOOST_CHECK_EQUAL(pop.getPopSize(), 0u);
}


// Check that no mortality leaves all individuals in the population
BOOST_AUTO_TEST_CASE(checkNoMortality)
{
    std::cout << "Testing absence of mortality...\n";
    Population pop = Population(10u);
    pop.survive(1.0);
    BOOST_CHECK_EQUAL(pop.getPopSize(), 10u);
}


// Check that a population has grown after reproduction
BOOST_AUTO_TEST_CASE(checkGrowth)
{
    std::cout << "Testing population growth...\n";
    Population pop = Population(10u);
    pop.reproduce(4.0);
    BOOST_CHECK(pop.getPopSize() > 10u);
}
