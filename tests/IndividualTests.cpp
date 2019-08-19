#include "library/Individual.h"
#include "library/utils.h"
#include "tests/GenFixture.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_FIXTURE_TEST_SUITE(indTestSuite, GenFixture)

    // Check the genome generation function
    BOOST_AUTO_TEST_CASE(checkGenomeGenerator)
    {
        std::cout << "Testing generating a random genome sequence...\n";
        Individual ind = Individual(genome, networks, 0.0);
        Diplotype seq = ind.getSequence();
        BOOST_CHECK_EQUAL(sumbool(seq[0u]), 0u);
        BOOST_CHECK_EQUAL(sumbool(seq[1u]), 0u);
    }

    // Check that a fully homogamous female will always accept identical mate
    BOOST_AUTO_TEST_CASE(checkAssortative)
    {
        std::cout << "Testing homogamous female...\n";
        Individual ind = Individual(genome, networks);
        ind.setEcoTrait(0.0, 1.0);
        ind.setMatePref(1.0);
        BOOST_CHECK(ind.acceptMate(0.0, 1.0));
    }

    // Check that a fully heterogamous female will always reject identical mate
    BOOST_AUTO_TEST_CASE(checkDisassortative)
    {
        std::cout << "Testing heterogamous female...\n";
        Individual ind = Individual(genome, networks);
        ind.setEcoTrait(0.0, 1.0);
        ind.setMatePref(-1.0);
        BOOST_CHECK(!ind.acceptMate(0.0, 1.0));
    }

    // Check that a random mating female will accept everybody
    BOOST_AUTO_TEST_CASE(checkRandomMating)
    {
        std::cout << "Testing random mating female...\n";
        Individual ind = Individual(genome, networks);
        ind.setEcoTrait(0.0, 1.0);
        ind.setMatePref(0.0);
        BOOST_CHECK(ind.acceptMate(0.0, 1.0));
        BOOST_CHECK(ind.acceptMate(1.0, 1.0));
        BOOST_CHECK(ind.acceptMate(-1.0, 1.0));
    }

    // We know exactly the fitness that a default individual should get
    BOOST_AUTO_TEST_CASE(checkFeeding)
    {
        std::cout << "Testing fitness obtained by feeding...\n";
        Individual ind = Individual(genome, networks);
        ind.setEcoTrait(-1.0, 1.0);
        ind.feed({ 100.0, 100.0 });
        BOOST_CHECK(ind.getFitness() == 0.04 + 0.0004 * exp(-4.0) * 100);
    }

BOOST_AUTO_TEST_SUITE_END()
