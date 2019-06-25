#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>
#include "doMain.h"

BOOST_AUTO_TEST_CASE(mainReturnsZeroExitCode)
{
    BOOST_CHECK_EQUAL(doMain(), 0);
}
