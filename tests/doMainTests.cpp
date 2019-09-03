#define BOOST_TEST_MAIN

#include "library/doMain.h"
#include "library/Population.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

typedef std::vector<Network> MultiNet;

// Check that the program can run without arguments
BOOST_AUTO_TEST_CASE(testUseNoArgs)
{
    std::cout << "Testing that the main runs without arguments...\n";
    BOOST_CHECK_EQUAL(doMain({ "program" }), 0);
}


// Check that the program cannot run with more than one argument
BOOST_AUTO_TEST_CASE(testAbuseTooManyArgs)
{
    std::cout << "Testing providing too many arguments to the main...\n";
    BOOST_CHECK_EQUAL(doMain({ "program", "arg1", "arg2" }), 1);
}


// Simulation should reach tmax in the absence of mortality
BOOST_AUTO_TEST_CASE(checkImmortalPopulation)
{

    std::cout << "Testing an immortal population...\n";

    const size_t tmax = 100u;
    const size_t initPopSize = 10u;
    const double dispersal = 0.0;
    const double survival = 1.0;
    const double birth = 0.0;
    const double mating = 0.0;

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();
    Population pop1 = Population(initPopSize, genome, networks);
    Population pop2 = Population(initPopSize, genome, networks);
    MetaPop metapop = {pop1, pop2};

    size_t t = runSimulation(metapop, tmax, dispersal, survival, birth, mating,
     genome, networks);

    BOOST_CHECK_EQUAL(t, tmax);
    BOOST_CHECK(metapop[0u].getPopSize() > 0u);
    BOOST_CHECK(metapop[1u].getPopSize() > 0u);
}


// Simulation should end prematurely with high mortality
BOOST_AUTO_TEST_CASE(checkProgressiveExtinction)
{

    std::cout << "Testing progressive extinction...\n";

    const size_t tmax = 100u;
    const size_t initPopSize = 10u;
    const double dispersal = 0.0;
    const double survival = 0.1;
    const double birth = 0.0;
    const double mating = 0.0;

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();
    Population pop1 = Population(initPopSize, genome, networks);
    Population pop2 = Population(initPopSize, genome, networks);
    MetaPop metapop = {pop1, pop2};

    size_t t = runSimulation(metapop, tmax, dispersal, survival, birth, mating,
     genome, networks);

    BOOST_CHECK(t < tmax);
    BOOST_CHECK(metapop[0u].getPopSize() == 0u);
    BOOST_CHECK(metapop[1u].getPopSize() == 0u);
}
