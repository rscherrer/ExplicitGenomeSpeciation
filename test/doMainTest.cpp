#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include "doMain.h"

/// Test of the default behavior of the main function
BOOST_AUTO_TEST_CASE(defaultMainReturnsZero)
{
    BOOST_CHECK_EQUAL(doMain( { "ExplicitGenomeSpeciation" } ), 0);

}

/// Test that main takes at most one argument
BOOST_AUTO_TEST_CASE(mainTakesOneArgAtMost)
{
    BOOST_CHECK_EQUAL(doMain( { "ExplicitGenomeSpeciation", "1", "2" } ), 1);
}
