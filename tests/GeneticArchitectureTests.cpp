#include "library/GeneticArchitecture.h"
#include "library/Random.h"
#include "library/utils.h"
#include "tests/testUtilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cassert>


// A genetic architecture created with 3 equal sized chromosomes
BOOST_AUTO_TEST_CASE(threeEqualSizedChromosomes) {

    std::vector<double> exp {(0.0 + 1.0) / 3.0, (1.0 + 1.0) / 3.0,
     (2.0 + 1.0) / 3.0};

    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    std::vector<double> real = arch.getChromosomeSizes();

    assert(exp.size() == real.size());

    BOOST_CHECK_EQUAL_COLLECTIONS(real.begin(), real.end(), exp.begin(),
     exp.end());
}

