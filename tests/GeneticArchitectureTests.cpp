#include "library/GeneticArchitecture.h"
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <cassert>


/// Parameters for a simple genetic architecture
struct SimpleParams
{

    SimpleParams() : pars(ParameterSet())
    {
        pars.setNChromosomes(3u);
        pars.setNTraits(3u);
        pars.setNLoci(10u);
        pars.setNLociPerTrait({ 3u, 3u, 4u });
        pars.setNEdgesPerTrait({ 1u, 0u, 2u });
        pars.setSkewnesses({ 1.0, 1.0, 1.0 });
        pars.setSeed(42u);

    }
    ~SimpleParams() {}

    ParameterSet pars;

};


/// A simple genetic architecture
struct SimpleArch
{
    SimpleArch() :
        simplepars(SimpleParams()),
        arch(GeneticArchitecture(simplepars.pars))
    {}
    ~SimpleArch() {}

    SimpleParams simplepars;
    GeneticArchitecture arch;
};


/// Function to compare vectors of doubles to their expectations
void checkVectorOfDoubles(std::vector<double> exp, std::vector<double> real, const double &factor = 1000000.0)
{
    for (size_t i = 0u; i < real.size(); ++i) {
        real[i] = round(real[i] * factor);
        exp[i] = round(exp[i] * factor);
    }

    BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(), exp.end());
}


/// Battery of tests using the same default genetic architecture
BOOST_FIXTURE_TEST_SUITE(testSuiteDefaultGeneticArchitectureParams, SimpleArch)

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
        std::vector<size_t> exp { 1u, 2u, 0u, 0u, 1u, 1u, 2u, 2u, 2u, 0u };
        std::vector<size_t> real = arch.getGenome().encodedTraits;
        assert(exp.size() == real.size());
        BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(), exp.end());
    }

    BOOST_AUTO_TEST_CASE(checkLocations)
    {
        std::vector<double> exp { 0.207918, 0.234415, 0.254931, 0.256963, 0.311207, 0.358894, 0.409707, 0.521725, 0.590133, 0.802339 };
        std::vector<double> real = arch.getGenome().locations;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }

    BOOST_AUTO_TEST_CASE(checkEffectSizes)
    {

        std::vector<double> exp { 0.950368, 0.283178, -0.468101, -0.337048, -0.948958, -0.154836, 0.128879, -0.816872, 0.0910685, 0.259251 };
        std::vector<double> real = arch.getGenome().effectSizes;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }

    BOOST_AUTO_TEST_CASE(checkDominanceCoeffs)
    {
        for (size_t locus = 0u; locus < simplepars.pars.getNLoci(); ++locus)
            std::cout << arch.getGenome().dominanceCoeffs[locus] << ' ';
        std::cout << '\n';

        std::vector<double> exp { 0.829453, 0.241587, 0.686339, 0.0740643, 0.442264, 0.132761, 0.54257, 0.737295, 0.510684, 0.681974 };
        std::vector<double> real = arch.getGenome().dominanceCoeffs;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
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
