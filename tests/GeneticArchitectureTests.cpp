#include "GeneticArchitecture.h"
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <cassert>


/// Set up a simple genetic architecture
struct defaultGeneticArchitectureParams
{

    defaultGeneticArchitectureParams() {
        pars.setNChromosomes(3u);
        pars.setNTraits(3u);
        pars.setNLociPerTrait({ 3u, 3u, 4u });
        pars.setNEdgesPerTrait({ 1u, 0u, 2u });
        pars.setSkewnesses({ 1.0, 1.0, 1.0 });
    }
    ~defaultGeneticArchitectureParams() {}

    ParameterSet pars;

};


/// Battery of tests using the same default genetic architecture
BOOST_FIXTURE_TEST_SUITE(testSuiteDefaultGeneticArchitecture, defaultGeneticArchitectureParams)

    /// A genetic architecture created with 3 equal sized chromosomes must have chromosome sizes 0.33, 0.66 and 1
    BOOST_AUTO_TEST_CASE(threeEqualSizedChromosomes)
    {
        GeneticArchitecture arch = GeneticArchitecture(pars);
        std::vector<double> expectedChromosomeSizes {(0.0 + 1.0) / 3.0, (1.0 + 1.0) / 3.0, (2.0 + 1.0) / 3.0};
        std::vector<double> realizedChromosomeSizes = arch.getChromosomeSizes();

        assert(expectedChromosomeSizes.size() == realizedChromosomeSizes.size());

        BOOST_CHECK_EQUAL_COLLECTIONS(
                realizedChromosomeSizes.begin(),
                realizedChromosomeSizes.end(),
                expectedChromosomeSizes.begin(),
                expectedChromosomeSizes.end()
                );
    }

    /// There should be one edge in trait 0, zero edges in trait 1 and two edges in trait 3
    BOOST_AUTO_TEST_CASE(checkNetworkEdges)
    {
        GeneticArchitecture arch = GeneticArchitecture(pars);
        BOOST_CHECK_EQUAL(arch.getTraitNetworkMaps()[0u].map[0u].first, 0u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworkMaps()[0u].map[0u].second, 1u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworkMaps()[2u].map[0u].first, 0u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworkMaps()[2u].map[0u].second, 1u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworkMaps()[2u].map[1u].first, 0u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworkMaps()[2u].map[1u].second, 2u);
    }

BOOST_AUTO_TEST_SUITE_END()


/*
GeneticArchitecture arch = GeneticArchitecture(nChromosomes, nTraits, nLociPerTrait, nEdgesPerTrait,
                                               skewnesses);

std::cout << "Architecture created\n";

for (size_t trait = 0u; trait < nTraits; ++trait)
{
std::cout << "trait = " << trait << '\n';
std::cout << "nEdges = " << arch.getTraitNetworkMaps()[trait].nEdges << '\n';
std::cout << "nVertices = " << arch.getTraitNetworkMaps()[trait].nVertices << '\n';
std::cout << "skewness = " << arch.getTraitNetworkMaps()[trait].skewness << '\n';
std::cout << "number of interactions realized = " << arch.getTraitNetworkMaps()[trait].map.size() << '\n';
for (size_t edge = 0u; edge < arch.getTraitNetworkMaps()[trait].map.size(); ++edge)
{
std::cout << "Edge " << edge << ": gene " << arch.getTraitNetworkMaps()[trait].map[edge].first <<
" and gene " << arch.getTraitNetworkMaps()[trait].map[edge].second << '\n';
}
}

 */



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