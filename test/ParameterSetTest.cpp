#define BOOST_TEST_DYN_LINK

#include "ParameterSet.h"
#include "TestFixtureCreateParameterFiles.h"
#include <boost/test/unit_test.hpp>


/// Test that parameters are read in correctly from a file
BOOST_FIXTURE_TEST_SUITE(testParametersAreReadInCorrectly, createValidParameterFile)

    BOOST_AUTO_TEST_CASE(readParametersCanReadParameters)
    {
        std::ifstream inputFile(validParameterFileName);
        ParameterSet parameters;
        parameters.readParameters(inputFile);
        BOOST_CHECK_EQUAL(parameters.initialPopSize, 1000u);
        BOOST_CHECK_EQUAL(parameters.birthRate, 4.0);
        const std::vector<double> expectedScaleA{ 0.0, 0.5, 1.0 };
        BOOST_CHECK_EQUAL_COLLECTIONS(parameters.scaleA.begin(), parameters.scaleA.end(), expectedScaleA.begin(), expectedScaleA.end());
    }

BOOST_AUTO_TEST_SUITE_END()
