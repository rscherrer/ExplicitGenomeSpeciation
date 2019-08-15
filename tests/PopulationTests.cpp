#include "library/Population.h"
#include "tests/GenFixture.h"
#include <boost/test/unit_test.hpp>
#include <iostream>


BOOST_FIXTURE_TEST_SUITE(popTestSuite, GenFixture)

    // Check that zero survival leaves no survivors in the population
    BOOST_AUTO_TEST_CASE(checkNoSurvivors)
    {
        std::cout << "Testing population wipe-out...\n";
        Population pop = Population(10u, genome);
        pop.survive(0.0);
        BOOST_CHECK_EQUAL(pop.getPopSize(), 0u);
    }


    // Check that no mortality leaves all individuals in the population
    BOOST_AUTO_TEST_CASE(checkNoMortality)
    {
        std::cout << "Testing absence of mortality...\n";
        Population pop = Population(10u, genome);
        pop.survive(1.0);
        BOOST_CHECK_EQUAL(pop.getPopSize(), 10u);
    }


    // Check that a population has grown after reproduction
    BOOST_AUTO_TEST_CASE(checkGrowth)
    {
        std::cout << "Testing population growth...\n";
        Population pop = Population(10u, genome);
        pop.reproduce(4.0, 1.0, genome);
        pop.survive(1.0);
        BOOST_CHECK(pop.getPopSize() >= 10u);
    }


    // Newborns should not die
    BOOST_AUTO_TEST_CASE(checkNewbornsShouldNotDie)
    {
        std::cout << "Testing that newborns do not die...\n";
        Population pop = Population(10u, genome);
        pop.reproduce(1.0, 1.0, genome);
        pop.survive(0.0);
        BOOST_CHECK(pop.getPopSize() > 0u);
    }

    // After feeding, the resources should be depleted
    BOOST_AUTO_TEST_CASE(checkResourceIsDepleted)
    {
        std::cout << "Testing resource depletion...\n";
        Population pop = Population(10u, genome);
        std::vector<double> resourcesBefore = pop.getResources();
        pop.consume();
        std::vector<double> resourcesAfter = pop.getResources();
        for (size_t res = 0u; res < 2u; ++res)
            BOOST_CHECK(resourcesBefore[res] > resourcesAfter[res]);
    }

BOOST_AUTO_TEST_SUITE_END()
