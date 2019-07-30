#include "library/Population.h"
#include <boost/test/unit_test.hpp>


// Check that zero survival leaves no survivors in the population
BOOST_AUTO_TEST_CASE(checkNoSurvivors)
{
    Population pop = Population(10u);
    pop.survive(0.0);
    BOOST_CHECK_EQUAL(pop.getPopSize(), 0u);
}


// Check that no mortality leaves all individuals in the population
