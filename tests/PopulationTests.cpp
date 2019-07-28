#define BOOST_TEST_DYN_LINK

#include "Population.h"
#include "Individual.h"
#include "TestFixtureMakeAPopulation.h"
#include <boost/test/unit_test.hpp>


/// A suite of tests of the habitat sorting functionality of class Population
BOOST_FIXTURE_TEST_SUITE(testHabitatSorting, makeADefaultPopulation)

    // Now test the correct sorting by habitat
    // Do I really need a Population, with protected access to everything, to test that function?
    // A simpler option might be to use a vector of individuals simply

    BOOST_AUTO_TEST_CASE(testThatFourZeroesStayTheSame) {

        // If habitats are 0 0 0 0 then habitats become 0 0 0 0

        expectedInitialHabitats = {0u, 0u, 0u, 0u};
        expectedSortedHabitats = {0u, 0u, 0u, 0u};

        setPopulationWithSortedHabitats();

        // Check that they are equal to what we wanted
        BOOST_CHECK_EQUAL_COLLECTIONS(sortedHabitats.begin(), sortedHabitats.end(), expectedSortedHabitats.begin(), expectedSortedHabitats.end());

    }

    BOOST_AUTO_TEST_CASE(testThatFourOnesStayTheSame) {

        // If habitats are 1 1 1 1 then habitats become 1 1 1 1

        expectedInitialHabitats = {1u, 1u, 1u, 1u};
        expectedSortedHabitats = {1u, 1u, 1u, 1u};

        setPopulationWithSortedHabitats();

        // Check that they are equal to what we wanted
        BOOST_CHECK_EQUAL_COLLECTIONS(sortedHabitats.begin(), sortedHabitats.end(), expectedSortedHabitats.begin(), expectedSortedHabitats.end());

    }

    BOOST_AUTO_TEST_CASE(testThatZeroOneZeroOneIsSortedInIncreasingOrder) {

        // If habitats are 0 1 0 1 then habitats become 0 0 1 1

        expectedInitialHabitats = {0u, 1u, 0u, 1u};
        expectedSortedHabitats = {0u, 0u, 1u, 1u};

        setPopulationWithSortedHabitats();

        BOOST_CHECK_EQUAL(1, 2); // Test is not run. Why?

        // Check that they are equal to what we wanted
        BOOST_CHECK_EQUAL_COLLECTIONS(sortedHabitats.begin(), sortedHabitats.end(), expectedSortedHabitats.begin(), expectedSortedHabitats.end());

    }

BOOST_AUTO_TEST_SUITE_END()

