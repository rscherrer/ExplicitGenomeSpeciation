#define BOOST_TEST_MAIN

#include "library/doMain.h"
#include "library/Population.h"
#include "library/MetaPop.h"
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

BOOST_AUTO_TEST_CASE(testAbuseInvalidFilename)
{
    std::cout << "Testing invalid parameter file name...\n";
    BOOST_CHECK_EQUAL(doMain({ "program", "nonsense.txt" }), 1);
}


// Simulation should reach tmax in the absence of mortality
BOOST_AUTO_TEST_CASE(checkImmortalPopulation)
{

    std::cout << "Testing an immortal population...\n";

    ParameterSet pars;
    pars.setTEndSim(100u);
    pars.setTSave(1u);
    pars.setInitialPopSize(10u);
    pars.setDispersalRate(0.0);
    pars.setSurvivalProb(1.0);
    pars.setBirthRate(0.0);
    pars.setMatePreferenceStrength(0.0);

    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();

    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(pars.getInitialPopSize(), genome, networks);

    MetaPop meta = MetaPop({ pop1, pop2 }, pars);

    size_t t = meta.evolve(genome, networks);

    BOOST_CHECK_EQUAL(t, pars.getTEndSim());
    BOOST_CHECK(meta.getPops()[0u].getPopSize() > 0u);
    BOOST_CHECK(meta.getPops()[1u].getPopSize() > 0u);
}


// Simulation should end prematurely with high mortality
BOOST_AUTO_TEST_CASE(checkProgressiveExtinction)
{

    std::cout << "Testing progressive extinction...\n";

    ParameterSet pars;
    pars.setTEndSim(100u);
    pars.setTSave(1u);
    pars.setInitialPopSize(10u);
    pars.setDispersalRate(0.0);
    pars.setSurvivalProb(0.1);
    pars.setBirthRate(0.0);
    pars.setMatePreferenceStrength(0.0);

    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    MultiNet networks = arch.getNetworks();

    Population pop1 = Population(pars.getInitialPopSize(), genome, networks);
    Population pop2 = Population(pars.getInitialPopSize(), genome, networks);

    MetaPop meta = MetaPop({ pop1, pop2 }, pars);

    size_t t = meta.evolve(genome, networks);

    BOOST_CHECK(t < pars.getTEndSim());
    BOOST_CHECK(meta.getPops()[0u].getPopSize() == 0u);
    BOOST_CHECK(meta.getPops()[1u].getPopSize() == 0u);
}
