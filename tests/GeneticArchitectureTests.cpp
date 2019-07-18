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
        real[i] = round(real[i] * factor);
        exp[i] = round(exp[i] * factor);
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
        std::vector<size_t> exp { 1u, 2u, 0u, 1u, 2u, 2u, 2u, 0u, 0u, 1u };
        std::vector<size_t> real = arch.getGenome().encodedTraits;
        assert(exp.size() == real.size());
        BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(), exp.end());
    }

    BOOST_AUTO_TEST_CASE(checkLocations)
    {
        std::vector<double> exp { 0.0245672, 0.121113, 0.462783, 0.485052, 0.498453, 0.557455, 0.668395, 0.792946,
                                  0.85107, 0.986548 };
        std::vector<double> real = arch.getGenome().locations;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }

    BOOST_AUTO_TEST_CASE(checkEffectSizes)
    {
        std::vector<double> exp { -0.76702, 0.428359, 0.72599, 0.14135, 0.323973, 0.553825, 0.414186, -0.854193,
                                  0.548988, 0.258618 };
        std::vector<double> real = arch.getGenome().effectSizes;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }

    BOOST_AUTO_TEST_CASE(checkDominanceCoeffs)
    {
        for (size_t locus = 0u; locus < simplepars.pars.getNLoci(); ++locus)
            std::cout << arch.getGenome().dominanceCoeffs[locus] << ' ';
        std::cout << '\n';

        std::vector<double> exp { 0.48246, 0.269249, 0.483417, 0.147075, 0.4715, 0.951774, 0.738187, 0.618192,
                                  0.617521, 0.0530536 };
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