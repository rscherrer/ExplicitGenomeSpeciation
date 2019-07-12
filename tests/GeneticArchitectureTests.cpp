#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE testGeneticArchitecture

#include "GeneticArchitecture.h"
//#include "TestFixtureCreateArchitectureFile.h"
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <cassert>


BOOST_AUTO_TEST_CASE(dummyTest)
{
    BOOST_CHECK_EQUAL(1, 1);
}


struct defaultGeneticArchitecture
{

    defaultGeneticArchitecture() : nChromosomes(3u), geneticArchitecture(GeneticArchitecture(nChromosomes)) {}
    ~defaultGeneticArchitecture() {}

    const size_t nChromosomes;
    GeneticArchitecture geneticArchitecture;

};


/// A genetic architecture created with 3 equal sized chromosomes must have chromosome sizes 0.33, 0.66 and 1
BOOST_FIXTURE_TEST_CASE(threeEqualSizedChromosomes, defaultGeneticArchitecture)
{

    std::vector<double> expectedChromosomeSizes {(0.0 + 1.0) / 3.0, (1.0 + 1.0) / 3.0, (2.0 + 1.0) / 3.0};
    std::vector<double> realizedChromosomeSizes = geneticArchitecture.getChromosomeSizes();

    assert(expectedChromosomeSizes.size() == realizedChromosomeSizes.size());

    BOOST_CHECK_EQUAL_COLLECTIONS(
            realizedChromosomeSizes.begin(),
            realizedChromosomeSizes.end(),
            expectedChromosomeSizes.begin(),
            expectedChromosomeSizes.end()
            );
}


/*
/// Test to check that the genetic architecture is read properly from a file
BOOST_FIXTURE_TEST_CASE(geneticArchitectureIsLoadedProperly, createValidArchitectureFile)
{

    GeneticArchitecture geneticArchitecture;
    parameters.architectureFileName = validArchitectureFileName;
    geneticArchitecture.loadGeneticArchitecture(parameters);

    std::clog << "Architecture tests.\n";

    // Check that the values are read in correctly
    BOOST_CHECK_EQUAL_COLLECTIONS(
            geneticArchitecture.chromosomeSizes.begin(),
            geneticArchitecture.chromosomeSizes.end(),
            expectedChromosomeSizes.begin(),
            expectedChromosomeSizes.end()
            );
    BOOST_CHECK_EQUAL(geneticArchitecture.locusConstants[0u].trait, 0u);
    BOOST_CHECK_EQUAL(geneticArchitecture.locusConstants[3u].location, 0.4);
    BOOST_CHECK_EQUAL(geneticArchitecture.locusConstants[4u].effectSize, 0.01);
    BOOST_CHECK_EQUAL(geneticArchitecture.locusConstants[5u].dominanceCoeff, 0.0);
    BOOST_CHECK_EQUAL(geneticArchitecture.locusConstants[4u].neighbors.front().first, 3u);
    BOOST_CHECK_EQUAL(geneticArchitecture.locusConstants[4u].neighbors.front().second, 0.01);

}
 */