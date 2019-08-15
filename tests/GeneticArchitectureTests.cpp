#include "library/GeneticArchitecture.h"
#include "library/Random.h"
#include "library/utils.h"
#include "tests/testUtilities.h"
#include <boost/test/unit_test.hpp>
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
        rnd::rng.seed(pars.getSeed());
    }
    ~SimpleParams() {}

    ParameterSet pars;

};


/// A simple genetic architecture
struct SimpleArch
{

    SimpleArch() : simplepars(SimpleParams()),
      arch(GeneticArchitecture(simplepars.pars)) {}
    ~SimpleArch() {}

    SimpleParams simplepars;
    GeneticArchitecture arch;
};


// Tests on a small genome
BOOST_FIXTURE_TEST_SUITE(smallGenomeTests, SimpleArch)

    // A genetic architecture created with 3 equal sized chromosomes must have
    // chromosome sizes 0.33, 0.66 and 1
    BOOST_AUTO_TEST_CASE(threeEqualSizedChromosomes) {

        std::vector<double> exp {(0.0 + 1.0) / 3.0, (1.0 + 1.0) / 3.0,
         (2.0 + 1.0) / 3.0};
        std::vector<double> real = arch.getChromosomeSizes();

        assert(exp.size() == real.size());

        BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(),
         exp.end());
    }


BOOST_AUTO_TEST_SUITE_END()
