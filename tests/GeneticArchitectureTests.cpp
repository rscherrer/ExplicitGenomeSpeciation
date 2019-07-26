#include "library/GeneticArchitecture.h"
#include "library/Random.h"
#include <boost/test/included/unit_test.hpp>
#include <iostream>
#include <cassert>


/// Parameters for a simple genetic architecture
struct SimpleParams
{

    SimpleParams() : pars(ParameterSet()) {
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

    SimpleArch() : simplepars(SimpleParams()),
     rnd(Random(simplepars.pars.getSeed())),
      arch(GeneticArchitecture(simplepars.pars, rnd)) {}
    ~SimpleArch() {}

    SimpleParams simplepars;
    Random rnd;
    GeneticArchitecture arch;
};


/// A parameter set with a thousand genes
struct BigParams {

    BigParams() : pars(ParameterSet()) {
        pars.setNChromosomes(3u);
        pars.setNTraits(3u);
        pars.setNLoci(1000u);
        pars.setNLociPerTrait({ 300u, 300u, 400u });
        pars.setNEdgesPerTrait({ 0u, 0u, 0u });
        pars.setSkewnesses({ 1.0, 1.0, 1.0 });
        pars.setSeed(42u);
    }
    ~BigParams() {}

    ParameterSet pars;
};


/// A genetic architecture with a thousand genes
struct BigArch {

    BigArch() : bigpars(BigParams()), rnd(Random(bigpars.pars.getSeed())),
     arch(GeneticArchitecture(bigpars.pars, rnd)) {}
    ~BigArch() {}

    BigParams bigpars;
    Random rnd;
    GeneticArchitecture arch;
};


/// Function to compare vectors of doubles to their expectations
void checkVectorOfDoubles(std::vector<double> exp, std::vector<double> real,
 const double &factor = 1.0e6)
{
    for (size_t i = 0u; i < real.size(); ++i) {
        real[i] = round(real[i] * factor);
        exp[i] = round(exp[i] * factor);
    }

    BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(),
     exp.end());
}


/// Battery of tests using the same default genetic architecture
BOOST_FIXTURE_TEST_SUITE(testSuiteDefaultGeneticArchitectureParams, SimpleArch)

    /// A genetic architecture created with 3 equal sized chromosomes must have
    /// chromosome sizes
    /// 0.33, 0.66 and 1
    BOOST_AUTO_TEST_CASE(threeEqualSizedChromosomes)
    {
        std::vector<double> exp {(0.0 + 1.0) / 3.0, (1.0 + 1.0) / 3.0,
         (2.0 + 1.0) / 3.0};
        std::vector<double> real = arch.getChromosomeSizes();

        assert(exp.size() == real.size());

        BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(),
         exp.end());
    }

    /// There should be one edge in trait 0, zero edges in trait 1 and two edges
    /// in trait 3
    BOOST_AUTO_TEST_CASE(checkNetworkEdges)
    {
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[0u].map[0u].first, 0u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[0u].map[0u].second, 1u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[2u].map[0u].first, 0u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[2u].map[0u].second, 1u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[2u].map[1u].first, 0u);
        BOOST_CHECK_EQUAL(arch.getTraitNetworks()[2u].map[1u].second, 2u);
    }

BOOST_AUTO_TEST_SUITE_END()
