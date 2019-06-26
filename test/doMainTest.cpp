#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE testTheMainFunction

#include "doMain.h"
#include "TestFixtureCreateParameterFiles.h"
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <string>


/// Test that the main can run without argument
BOOST_AUTO_TEST_CASE(mainCanRunWithoutArgument)
{
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

// Wait, what is a genetic architecture file made of?

// A text file with the following
// The number of loci underlying the ecological, mating and neutral trait, respectively, each separated by a new line
// After a new line, the number of epistatic interactions underlying the ecological, mating and neutral trait, respectively
// After a new line, the number of chromosomes
// After a new line, the sizes of each chromosome separated by tabulations
// After a new line, for each locus separated by a new line, the index of the trait underlain by the locus, its location, its effect size and its dominance coefficient, separated by tabulations
// After a new line, for each pair of interacting loci separated by a new line, the index of the first locus, the index of the second locus and the interaction weight, separated by tabulations
// Note: each pair of interacting loci appears only once in the file

// Number of loci for each trait
2
2
2

// Number of interactions for each trait
1
1
1

// Number of chromosomes
1

// Size of each chromosome
1

// Characteristics of each locus
1
2
3
4
5
6

// Characteristics of each interaction

/// Test fixture for creating a dummy valid genetic architecture
struct createValidArchitectureFile
{

    createValidArchitectureFile() : validArchitectureFileName("valid_architecture.txt") {
        validArchitectureFile.open(validArchitectureFileName);
        validArchitectureFile <<
    }
    ~createValidArchitectureFile() {}

    std::ofstream validArchitectureFile;
    std::string validArchitectureFileName;

};

// Create a dummy invalid genetic architecture for testing purposes


// Check that doMain errors if an the genetic architecture file does not exist
// Check that doMain errors if a genetic architecture file has parameters that do not match previously set parameters
// Check that doMain errors if the genetic architecture file does not have the right dimensions

// Check that the genetic architecture is read properly