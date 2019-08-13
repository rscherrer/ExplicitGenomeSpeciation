#include "library/Individual.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

// Check that a fully homogamous female will always accept identical mate
BOOST_AUTO_TEST_CASE(checkAssortative)
{
    std::cout << "Testing homogamous female...\n";
    Individual ind = Individual();
    ind.setEcoTrait(0.0);
    ind.setMatePref(1.0);
    BOOST_CHECK(ind.acceptMate(0.0, 1.0));
}

// Check that a fully heterogamous female will always reject identical mate
BOOST_AUTO_TEST_CASE(checkDisassortative)
{
    std::cout << "Testing heterogamous female...\n";
    Individual ind = Individual();
    ind.setEcoTrait(0.0);
    ind.setMatePref(-1.0);
    BOOST_CHECK(!ind.acceptMate(0.0, 1.0));
}

// Check that a random mating female will accept everybody
BOOST_AUTO_TEST_CASE(checkRandomMating)
{
    std::cout << "Testing random mating female...\n";
    Individual ind = Individual();
    ind.setEcoTrait(0.0);
    ind.setMatePref(0.0);
    BOOST_CHECK(ind.acceptMate(0.0, 1.0));
    BOOST_CHECK(ind.acceptMate(1.0, 1.0));
    BOOST_CHECK(ind.acceptMate(-1.0, 1.0));
}
