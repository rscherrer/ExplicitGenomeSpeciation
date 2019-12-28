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
    std::clog << "Testing run without arguments...\n";
    BOOST_CHECK_EQUAL(simulate({ "EGS_test" }), 0);
}

// Check that the program cannot run with more than one argument
BOOST_AUTO_TEST_CASE(testAbuseTooManyArgs)
{
    std::clog << "Testing run with too many arguments...\n";
    BOOST_CHECK_EQUAL(simulate({ "EGS_test", "arg1", "arg2" }), 1);
}

BOOST_AUTO_TEST_CASE(testAbuseInvalidFilename)
{
    std::clog << "Testing run with invalid parameter file...\n";
    BOOST_CHECK_EQUAL(simulate({ "EGS_test", "nonsense.txt" }), 1);
}


BOOST_AUTO_TEST_CASE(testUseValidFilename)
{
    std::clog << "Testing run with valid parameter file...\n";
    tst::makeValidParamFile();
    BOOST_CHECK_EQUAL(simulate({ "EGS_test", "validparamfile.txt" }), 0);
}


BOOST_AUTO_TEST_CASE(testAbuseInvalidParamName)
{
    std::clog << "Testing run with invalid parameter names...\n";
    tst::makeInvalidParamName();
    BOOST_CHECK_EQUAL(simulate({ "EGS_test", "invalidparamname.txt" }), 1);
}

BOOST_AUTO_TEST_CASE(testAbuseInvalidParamValue)
{
    std::clog << "Testing run with invalid parameter values...\n";
    tst::makeInvalidParamValue();
    BOOST_CHECK_EQUAL(simulate({"EGS_test", "invalidparamvalue.txt"}), 1);
    tst::makeInvalidParamValue2();
    BOOST_CHECK_EQUAL(simulate({"EGS_test", "invalidparamvalue2.txt"}), 1);
}

/*
BOOST_AUTO_TEST_CASE(testUseLongSim) // takes a while
{
    std::clog << "Testing run for a hundred generations...\n";
    tst::makeValidParamFile2();
    BOOST_CHECK_EQUAL(simulate({ "EGS_test", "validparamfile2.txt" }), 0);
}
*/



