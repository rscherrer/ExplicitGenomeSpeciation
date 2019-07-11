#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE testTheMainFunction

#include "doMain.h"
#include "TestFixtureCreateParameterFiles.h"
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <string>
#include <iostream>


/// Test that the main can run without argument
BOOST_AUTO_TEST_CASE(mainCanRunWithoutArgument)
{
    BOOST_TEST_MESSAGE("Testing the main function...\n");
    BOOST_CHECK_EQUAL(doMain( { "ExplicitGenomeSpeciation" } ), 0);
}


/// Test that main cannot take more than one argument
BOOST_AUTO_TEST_CASE(mainCannotTakeMoreThanOneArgument)
{
    BOOST_CHECK_EQUAL(doMain( { "ExplicitGenomeSpeciation", "1", "2" } ), 1);
}


/// Test that main errors if the argument is not a valid input file name
BOOST_AUTO_TEST_CASE(mainCannotTakeInvalidInputFileAsArgument)
{
    BOOST_CHECK_EQUAL(doMain( { "ExplicitGenomeSpeciation", "nonsense" } ), 1);
}


/// Test the behavior of main in relation to input file reading
BOOST_FIXTURE_TEST_SUITE(testMainReadingAbilities, createValidParameterFile)

    /// Test that main can run if a valid input file name is provided
    BOOST_AUTO_TEST_CASE(mainCanTakeValidInputFileAsArgument)
    {
        BOOST_CHECK_EQUAL(doMain( { "ExplicitGenomeSpeciation", "valid_parameters.txt" }), 0);
    }

BOOST_AUTO_TEST_SUITE_END()


/// Test that main errors if invalid parameters are provided in a parameter file
BOOST_FIXTURE_TEST_SUITE(testMainWithInvalidParameters, createInvalidParameterFile)

    BOOST_AUTO_TEST_CASE(mainCanDetectInvalidParametersFromAFile)
    {
        BOOST_CHECK_EQUAL(doMain( { "ExplicitGenomeSpeciation", "invalid_parameters.txt" } ), 1);
    }

BOOST_AUTO_TEST_SUITE_END()



// Characteristics of each interaction



// Create a dummy invalid genetic architecture for testing purposes


// Check that doMain errors if an the genetic architecture file does not exist
// Check that doMain errors if a genetic architecture file has parameters that do not match previously set parameters
// Check that doMain errors if the genetic architecture file does not have the right dimensions

