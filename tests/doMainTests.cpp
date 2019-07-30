#define BOOST_TEST_MAIN

#include "library/doMain.h"
#include "library/Population.h"
#include <boost/test/unit_test.hpp>
#include <iostream>


// Check that the program can run without arguments
BOOST_AUTO_TEST_CASE(testUseNoArgs)
{
    BOOST_CHECK_EQUAL(doMain({ "program" }), 0);
}


// Check that the program cannot run with more than one argument
BOOST_AUTO_TEST_CASE(testAbuseTooManyArgs)
{
    BOOST_CHECK_EQUAL(doMain({ "program", "arg1", "arg2" }), 1);
}


// After a simulation ends, the time should be tmax
BOOST_AUTO_TEST_CASE(checkTimeAfterSimul)
{
    size_t t = 0u;
    size_t tmax = 100u;

    ParameterSet pars;
    Population pop = Population(pars.getInitialPopSize());

    runSimulation(t, pop, tmax);

    BOOST_CHECK_EQUAL(t, tmax);

    // After a simulation, the population should still be there and
    // have individuals
    BOOST_CHECK(pop.getPopSize() > 0u);
}


