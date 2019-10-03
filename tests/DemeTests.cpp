#include "library/Deme.h"
#include "tests/DemeFixture.h"
#include <boost/test/unit_test.hpp>
#include <iostream>


BOOST_FIXTURE_TEST_SUITE(popTestSuite, PopFixture)

    /*
    BOOST_AUTO_TEST_CASE(checkNoDispersal)
    {
        Deme pop = Deme(10u, s, max, k, r, arch);
        Crowd migrants = pop.emigrate(0.0);
        BOOST_CHECK_EQUAL(migrants.size(), 0u);
        BOOST_CHECK_EQUAL(pop.getPopSize(), 10u);
    }

    BOOST_AUTO_TEST_CASE(checkExodus)
    {
        Deme pop = Deme(10u, s, max, k, r, arch);
        Crowd migrants = pop.emigrate(1.0);
        BOOST_CHECK_EQUAL(migrants.size(), 10u);
        BOOST_CHECK_EQUAL(pop.getPopSize(), 0u);
    }

    BOOST_AUTO_TEST_CASE(checkBorderControl)
    {
        Deme pop = Deme(10u, s, max, k, r, arch);
        Crowd migrants = pop.emigrate(1.0);
        pop.immigrate(migrants);
        BOOST_CHECK_EQUAL(pop.getPopSize(), 10u);

    }
    */

    // Check that zero survival leaves no survivors in the population
    BOOST_AUTO_TEST_CASE(checkNoSurvivors)
    {
        Deme pop = Deme(10u, s, max, k, r, arch);
        pop.survive(0.0);
        BOOST_CHECK_EQUAL(pop.getPopSize(), 0u);
    }


    // Check that no mortality leaves all individuals in the population
    BOOST_AUTO_TEST_CASE(checkNoMortality)
    {
        Deme pop = Deme(10u, s, max, k, r, arch);
        pop.survive(1.0);
        BOOST_CHECK_EQUAL(pop.getPopSize(), 10u);
    }


    // Check that a population has grown after reproduction
    BOOST_AUTO_TEST_CASE(checkGrowth)
    {
        Deme pop = Deme(10u, s, max, k, r, arch);
        pop.reproduce(4.0, 0.0, 0.01, 1.0, 4.0E-4, arch);
        pop.survive(1.0);
        BOOST_CHECK(pop.getPopSize() >= 10u);
    }


    // Newborns should not die
    BOOST_AUTO_TEST_CASE(checkNewbornsShouldNotDie)
    {
        Deme pop = Deme(100u, s, max, k, r, arch);
        pop.reproduce(4.0, 0.0, 0.01, 1.0, 4.0E-4, arch);
        pop.survive(0.0); // kill all adults        
        BOOST_CHECK(pop.getPopSize() > 0u);
    }

    // After feeding, the resources should be depleted
    BOOST_AUTO_TEST_CASE(checkResourceIsDepleted)
    {
        Deme pop = Deme(10u, s, max, k, r, arch);
        vecDbl before = { pop.getResource(0u), pop.getResource(1u) };
        pop.consume();
        vecDbl after = { pop.getResource(0u), pop.getResource(1u) };

        for (size_t res = 0u; res < 2u; ++res)
            BOOST_CHECK(before[res] > after[res]);
    }

BOOST_AUTO_TEST_SUITE_END()
