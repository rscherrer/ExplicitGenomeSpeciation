#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include "doMain.h"


/// Good usage of the function
BOOST_AUTO_TEST_SUITE(doMainUse)

    /// Test of the default behavior of the main function
    BOOST_AUTO_TEST_CASE(defaultMainReturnsZero)
    {
        BOOST_CHECK_EQUAL(doMain( { "ExplicitGenomeSpeciation" } ), 0);

    }

BOOST_AUTO_TEST_SUITE_END()


/// Bad usage of the function
BOOST_AUTO_TEST_SUITE(doMainAbuse)

    /// Test that main takes at most one argument
    BOOST_AUTO_TEST_CASE(mainTakesOneArgAtMost)
    {
        BOOST_CHECK_EQUAL(doMain( { "ExplicitGenomeSpeciation", "1", "2" } ), 1);
    }

BOOST_AUTO_TEST_SUITE_END()


