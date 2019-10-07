#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include "library/Simul.h"
#include "tests/TestUtilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

// Black box testing of the proper run of the main function

// Check that the program can run without arguments
BOOST_AUTO_TEST_CASE(testUseNoArgs)
{
    BOOST_CHECK_EQUAL(simulate({ "EGS_test" }), 0);
}

// Check that the program cannot run with more than one argument
BOOST_AUTO_TEST_CASE(testAbuseTooManyArgs)
{
    BOOST_CHECK_EQUAL(simulate({ "EGS_test", "arg1", "arg2" }), 1);
}

BOOST_AUTO_TEST_CASE(testAbuseInvalidFilename)
{
    tst::makeValidParamFile();
    BOOST_CHECK_EQUAL(simulate({ "EGS_test", "nonsense.txt" }), 1);
}


BOOST_AUTO_TEST_CASE(testUseValidFilename)
{
    BOOST_CHECK_EQUAL(simulate({ "EGS_test", "valid_paramfile_test.txt" }), 0);
}


BOOST_AUTO_TEST_CASE(testAbuseInvalidParamName)
{
    tst::makeInvalidParamName();
    BOOST_CHECK_EQUAL(simulate({ "EGS_test", "invalid_paramname_test.txt" }), 1);
}

BOOST_AUTO_TEST_CASE(testAbuseInvalidParamValue)
{
    tst::makeInvalidParamValue();
    BOOST_CHECK_EQUAL(simulate({"EGS_test", "invalid_paramvalue_test.txt"}), 1);
    tst::makeInvalidParamValue2();
    BOOST_CHECK_EQUAL(simulate({"EGS_test", "invalid_paramvalue_test2.txt"}), 1);
}
