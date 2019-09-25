#include "library/Individual.h"
#include "library/Utilities.h"
#include "tests/GenFixture.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_FIXTURE_TEST_SUITE(indTestSuite, GenFixture)

    // Check the genome generation function
    BOOST_AUTO_TEST_CASE(checkOnlyAllelesZero)
    {
        std::clog << "Testing generating a sequence of only zeros...\n";
        Individual ind = Individual(arch, 1.0, 4.0E-4, 0.0);
        BOOST_CHECK_EQUAL(ind.getAlleleSum(0u), 0u);
        BOOST_CHECK_EQUAL(ind.getAlleleSum(1u), 0u);
    }

    // Check the genome generation function
    BOOST_AUTO_TEST_CASE(checkOnlyAllelesOne)
    {
        std::clog << "Testing generating a sequence of only ones...\n";
        Individual ind = Individual(arch, 1.0, 4.0E-4, 1.0);
        BOOST_CHECK_EQUAL(ind.getAlleleSum(0u), arch.nLoci);
        BOOST_CHECK_EQUAL(ind.getAlleleSum(1u), arch.nLoci);
    }

    // Check that a fully homogamous female will always accept identical mate
    BOOST_AUTO_TEST_CASE(checkAssortative)
    {
        std::clog << "Testing homogamous female...\n";
        Individual ind = Individual(arch, 1.0, 4.0E-4, 0.5);
        ind.setEcoTrait(0.0, 1.0, 4.0E-4);
        ind.setMatePref(1.0);
        BOOST_CHECK(ind.acceptMate(0.0, 1.0));
    }

    // Check that a fully heterogamous female will always reject identical mate
    BOOST_AUTO_TEST_CASE(checkDisassortative)
    {
        std::clog << "Testing heterogamous female...\n";
        Individual ind = Individual(arch, 1.0, 4.0E-4, 0.5);
        ind.setEcoTrait(0.0, 1.0, 4.0E-4);
        ind.setMatePref(-1.0);
        BOOST_CHECK(!ind.acceptMate(0.0, 1.0));
    }

    // Check that a random mating female will accept everybody
    BOOST_AUTO_TEST_CASE(checkRandomMating)
    {
        std::clog << "Testing random mating female...\n";
        Individual ind = Individual(arch, 1.0, 4.0E-4, 0.5);
        ind.setEcoTrait(0.0, 1.0, 4.0E-4);
        ind.setMatePref(0.0);
        BOOST_CHECK(ind.acceptMate(0.0, 1.0));
        BOOST_CHECK(ind.acceptMate(1.0, 1.0));
        BOOST_CHECK(ind.acceptMate(-1.0, 1.0));
    }

    // We know exactly the fitness that a default individual should get
    BOOST_AUTO_TEST_CASE(checkFeeding)
    {
        std::clog << "Testing fitness obtained by feeding...\n";
        Individual ind = Individual(arch, 1.0, 4.0E-4, 0.5);
        ind.setEcoTrait(-1.0, 1.0, 4.0E-4);
        ind.feed({ 100.0, 100.0 });
        BOOST_CHECK(ind.getFitness() == 0.04 + 4.0E-4 * exp(-4.0) * 100);
    }

    // Check that fecundation fuses the two gametes
    BOOST_AUTO_TEST_CASE(checkFecundation)
    {
        std::clog << "Testing that fecundation fuses two gametes...\n";
        Individual mom = Individual(arch,1.0, 4.0E-4, 0.0);
        Individual dad = Individual(arch, 1.0, 4.0E-4, 1.0);
        Haplotype egg = mom.recombine(arch);
        Haplotype sperm = dad.recombine(arch);
        Individual baby = Individual(arch, egg, sperm, 1.0, 4.0E-4);
        BOOST_CHECK_EQUAL(baby.getAlleleSum(0u), 0u);
        BOOST_CHECK_EQUAL(baby.getAlleleSum(1u), arch.nLoci);
    }

    BOOST_AUTO_TEST_CASE(checkNoMutation)
    {
        std::clog << "Testing absence of mutations...\n";
        Individual ind = Individual(arch, 1.0, 4.0E-4, 0.0);
        Haplotype gamete = ind.recombine(arch);
        ind.mutate(gamete, 0.0);
        BOOST_CHECK_EQUAL(gamete.count(), 0u);
    }

    BOOST_AUTO_TEST_CASE(checkHighMutation)
    {
        std::clog << "Testing high mutation rate...\n";
        Individual ind = Individual(arch, 1.0, 4.0E-4, 0.0);
        Haplotype gamete = ind.recombine(arch);
        ind.mutate(gamete, 100.0);
        BOOST_CHECK(gamete.count() != 0u);
    }

    BOOST_AUTO_TEST_CASE(checkHighRecombination)
    {
        std::clog << "Testing high recombination rate...\n";
        Individual mom = Individual(arch, 1.0, 4.0E-4, 0.0);
        Individual dad = Individual(arch, 1.0, 4.0E-4, 1.0);
        Haplotype egg = mom.recombine(arch);
        Haplotype sperm = dad.recombine(arch);
        Individual baby = Individual(arch, egg, sperm, 1.0, 4.0E-4);
        arch.recombinationRate = 10.0;
        Haplotype gamete = baby.recombine(arch);
        BOOST_CHECK(gamete.count() != 0u);
        BOOST_CHECK(gamete.count() != arch.nLoci);
    }

    BOOST_AUTO_TEST_CASE(checkExpression)
    {
        std::clog << "Testing gene expression...\n";
        Individual ind1 = Individual(arch, 1.0, 4.0E-4, 0.0);
        vecDbl expression1 = ind1.getExpression();
        BOOST_CHECK_EQUAL(utl::sum(expression1), -1.0 * arch.nLoci);
        Individual ind2 = Individual(arch, 1.0, 4.0E-4, 1.0);
        vecDbl expression2 = ind2.getExpression();
        BOOST_CHECK_EQUAL(utl::sum(expression2), arch.nLoci);
    }

BOOST_AUTO_TEST_SUITE_END()

// Check the absence of recombination
BOOST_AUTO_TEST_CASE(checkNoRecombination)
{
    std::clog << "Testing meiosis without recombination...\n";
    Param pars;
    pars.setNChromosomes(1u); // to avoid free recombination
    GenArch arch = GenArch(pars);
    Individual mom = Individual(arch, 1.0, 4.0E-4, 0.0);
    Individual dad = Individual(arch, 1.0, 4.0E-4, 1.0);
    Haplotype egg = mom.recombine(arch);
    Haplotype sperm = dad.recombine(arch);
    Individual baby = Individual(arch, egg, sperm, 1.0, 4.0E-4);
    Haplotype gam = baby.recombine(arch);
    BOOST_CHECK_EQUAL(gam.count() / arch.nLoci, gam[0u]);
}

BOOST_AUTO_TEST_CASE(checkDevelopment)
{
    std::clog << "Testing developing individual...\n";
    Param pars;
    pars.setDominanceVariance(0.0);
    GenArch arch = GenArch(pars);
    Individual mom = Individual(arch, 1.0, 4.0E-4, 0.0);
    Individual dad = Individual(arch, 1.0, 4.0E-4, 1.0);
    Haplotype egg = mom.recombine(arch);
    Haplotype sperm = dad.recombine(arch);
    Individual baby = Individual(arch, egg, sperm, 1.0, 4.0E-4);
    BOOST_CHECK_EQUAL(baby.getEcoTrait(), 0.0);
    BOOST_CHECK_EQUAL(baby.getMatePref(), 0.0);
    BOOST_CHECK_EQUAL(baby.getNeutral(), 0.0);
}






