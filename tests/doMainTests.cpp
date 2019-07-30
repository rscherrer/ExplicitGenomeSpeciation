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
BOOST_AUTO_TEST_CASE(checkSimulReachesTMax)
{
    size_t t = 0u;
    const size_t tmax = 100u;
    const size_t initPopSize = 10u;
    const double survival = 1.0;

    ParameterSet pars;
    Population pop = Population(initPopSize);

    runSimulation(t, pop, tmax, survival);

    BOOST_CHECK_EQUAL(t, tmax);
}




