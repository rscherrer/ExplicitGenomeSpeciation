#define BOOST_TEST_MAIN

#include "library/doMain.h"
#include <boost/test/unit_test.hpp>


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
