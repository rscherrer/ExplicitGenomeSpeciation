#include "GeneticArchitecture.h"
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <cassert>


/// Parameters for a simple genetic architecture
struct simpleParams
{

    simpleParams()
    {
        pars.setNChromosomes(3u);
        pars.setNTraits(3u);
        pars.setNLoci(10u);
        pars.setNLociPerTrait({ 3u, 3u, 4u });
        pars.setNEdgesPerTrait({ 1u, 0u, 2u });
        pars.setSkewnesses({ 1.0, 1.0, 1.0 });
        pars.setSeed(42u);
    }
    ~simpleParams() {}

    ParameterSet pars;

};


/// A simple genetic architecture
struct simpleArch
{
    simpleArch() : arch(GeneticArchitecture(simplepars.pars)) {}
    ~simpleArch() {}

    simpleParams simplepars;
    GeneticArchitecture arch;
};


/// Function to compare vectors of doubles to their expectations
void checkVectorOfDoubles(std::vector<double> exp, std::vector<double> real, const double &factor = 1000000.0)
{
    for (size_t i = 0u; i < real.size(); ++i) {
        real[i] = round(real[i] * 1000000.0);
        exp[i] = round(exp[i] * 1000000.0);
    }

    BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(), exp.end());
}


/// Battery of tests using the same default genetic architecture
BOOST_FIXTURE_TEST_SUITE(testSuiteDefaultGeneticArchitectureParams, simpleArch)

    /// A genetic architecture created with 3 equal sized chromosomes must have chromosome sizes 0.33, 0.66 and 1
    BOOST_AUTO_TEST_CASE(threeEqualSizedChromosomes)
    {
        std::vector<double> exp {(0.0 + 1.0) / 3.0, (1.0 + 1.0) / 3.0, (2.0 + 1.0) / 3.0};
        std::vector<double> real = arch.getChromosomeSizes();

        assert(exp.size() == real.size());

        BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(), exp.end());
    }

    /// There should be one edge in trait 0, zero edges in trait 1 and two edges in trait 3
    BOOST_AUTO_TEST_CASE(checkNetworkEdges)
    {
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[0u].map[0u].first, 0u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[0u].map[0u].second, 1u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[2u].map[0u].first, 0u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[2u].map[0u].second, 1u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[2u].map[1u].first, 0u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[2u].map[1u].second, 2u);
    }

    BOOST_AUTO_TEST_CASE(checkEncodedTraits)
    {
        std::vector<size_t> exp { 1u, 1u, 0u, 2u, 2u, 2u, 0u, 2u, 1u, 0u };
        std::vector<size_t> real = arch.getGenome().encodedTraits;
        assert(exp.size() == real.size());
        BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(), exp.end());
    }

    BOOST_AUTO_TEST_CASE(checkLocations)
    {
        std::vector<double> exp { 0.00717771, 0.0634439, 0.141306, 0.234415, 0.50963, 0.569106, 0.590133, 0.602823,
                                  0.756531, 0.960691 };
        std::vector<double> real = arch.getGenome().locations;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }

    BOOST_AUTO_TEST_CASE(checkEffectSizes)
    {
        std::vector<double> exp { -0.401348, 0.112783, 0.690384, -0.601905, -0.757984, 0.158833, 0.569097, -0.632165,
                                  0.114488, -0.806784 };
        std::vector<double> real = arch.getGenome().effectSizes;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }

    BOOST_AUTO_TEST_CASE(checkDominanceCoeffs)
    {
        std::vector<double> exp { 0.948098, 0.00938372, 0.62366, 0.781639, 0.367136, 0.335983, 0.310815, 0.0671098,
                                  0.432369, 0.75192 };
        std::vector<double> real = arch.getGenome().dominanceCoeffs;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }

    BOOST_AUTO_TEST_CASE(checkInteractionWeights)
    {
        
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