#define BOOST_TEST_MAIN

#include "library/doMain.h"
#include "library/Deme.h"
#include "library/MetaPop.h"
#include "tests/GenFixture.h"
#include "tests/testUtilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

typedef std::vector<Network> MultiNet;

// Check that the program can run without arguments
BOOST_AUTO_TEST_CASE(testUseNoArgs)
{
    std::cout << "Testing that the main runs without arguments...\n";
    BOOST_CHECK_EQUAL(doMain({ "EGS_test" }), 0);
}


// Check that the program cannot run with more than one argument
BOOST_AUTO_TEST_CASE(testAbuseTooManyArgs)
{
    std::cout << "Testing providing too many arguments to the main...\n";
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "arg1", "arg2" }), 1);
}

BOOST_AUTO_TEST_CASE(testAbuseInvalidFilename)
{
    std::cout << "Testing invalid parameter file name...\n";
    makeValidParamFile();
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "nonsense.txt" }), 1);
}


BOOST_AUTO_TEST_CASE(testUseValidFilename)
{
    std::cout << "Testing valid parameter file name...\n";
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "valid_paramfile_test.txt" }), 0);
}


BOOST_AUTO_TEST_CASE(testAbuseInvalidParamName)
{
    std::cout << "Testing invalid parameter name...\n";
    makeInvalidParamName();
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "invalid_paramname_test.txt" }), 1);
}

BOOST_AUTO_TEST_CASE(testAbuseInvalidParamValue)
{
    std::cout << "Testing invalid parameter value...\n";
    makeInvalidParamValue();
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "invalid_paramvalue_test.txt" }), 1);
    makeInvalidParamValue2();
    BOOST_CHECK_EQUAL(doMain({ "EGS_test", "invalid_paramvalue_test.txt" }), 1);
}
