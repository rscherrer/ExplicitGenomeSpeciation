#include "library/GeneticArchitecture.h"
#include "library/Random.h"
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
        rnd(Random(simplepars.pars.getSeed())),
        arch(GeneticArchitecture(simplepars.pars, rnd))
    {}
    ~SimpleArch() {}

    SimpleParams simplepars;
    Random rnd;
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

    //

    /*
    BOOST_AUTO_TEST_CASE(checkEncodedTraits)
    {
        std::vector<size_t> exp { 2u, 2u, 0u, 2u, 2u, 1u, 1u, 0u, 0u, 1u };
        std::vector<size_t> real = arch.getGenome().encodedTraits;
        assert(exp.size() == real.size());
        BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(), exp.end());
    }

    BOOST_AUTO_TEST_CASE(checkLocations)
    {
        std::vector<double> exp { 0.0645945, 0.138685, 0.402559, 0.416049, 0.443436, 0.444609, 0.493076, 0.686054, 0.735338, 0.773019 };
        std::vector<double> real = arch.getGenome().locations;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }

    BOOST_AUTO_TEST_CASE(checkEffectSizes)
    {
        std::vector<double> exp { -0.151555, 0.337905, -0.280014, -0.716991, -0.590573, -0.192771, 0.973196, 0.328876, 0.901905, 0.125417 };
        std::vector<double> real = arch.getGenome().effectSizes;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }


    BOOST_AUTO_TEST_CASE(checkDominanceCoeffs)
    {
        std::vector<double> exp { 0.997809, 0.0533615, 0.269254, 0.0368082, 0.0132147, 0.555865, 0.822442, 0.958672, 0.0919276, 0.120843 };
        std::vector<double> real = arch.getGenome().dominanceCoeffs;
        assert(exp.size() == real.size());
        checkVectorOfDoubles(exp, real);
    }
    */

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
