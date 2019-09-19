#include "library/Individual.h"
#include "library/utils.h"
#include "tests/GenFixture.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_FIXTURE_TEST_SUITE(indTestSuite, GenFixture)

    // Check the genome generation function
    BOOST_AUTO_TEST_CASE(checkOnlyAllelesZero)
    {
        std::cout << "Testing generating a sequence of only zeros...\n";
        Individual ind = Individual(arch, 0.0);
        BOOST_CHECK_EQUAL(sumbool(ind.getSequence(0u)), 0u);
        BOOST_CHECK_EQUAL(sumbool(ind.getSequence(1u)), 0u);
    }

    // Check the genome generation function
    BOOST_AUTO_TEST_CASE(checkOnlyAllelesOne)
    {
        std::cout << "Testing generating a sequence of only ones...\n";
        Individual ind = Individual(arch, 1.0);
        BOOST_CHECK_EQUAL(sumbool(ind.getSequence(0u)), arch.nLoci);
        BOOST_CHECK_EQUAL(sumbool(ind.getSequence(1u)), arch.nLoci);
    }

    // Check that a fully homogamous female will always accept identical mate
    BOOST_AUTO_TEST_CASE(checkAssortative)
    {
        std::cout << "Testing homogamous female...\n";
        Individual ind = Individual(arch);
        ind.setEcoTrait(0.0, 1.0);
        ind.setMatePref(1.0);
        BOOST_CHECK(ind.acceptMate(0.0, 1.0));
    }

    // Check that a fully heterogamous female will always reject identical mate
    BOOST_AUTO_TEST_CASE(checkDisassortative)
    {
        std::cout << "Testing heterogamous female...\n";
        Individual ind = Individual(arch);
        ind.setEcoTrait(0.0, 1.0);
        ind.setMatePref(-1.0);
        BOOST_CHECK(!ind.acceptMate(0.0, 1.0));
    }

    // Check that a random mating female will accept everybody
    BOOST_AUTO_TEST_CASE(checkRandomMating)
    {
        std::cout << "Testing random mating female...\n";
        Individual ind = Individual(arch);
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
        Individual ind = Individual(arch);
        ind.setEcoTrait(-1.0, 1.0);
        ind.feed({ 100.0, 100.0 });
        BOOST_CHECK(ind.getFitness() == 0.04 + 0.0004 * exp(-4.0) * 100);
    }

    // Check that fecundation fuses the two gametes
    BOOST_AUTO_TEST_CASE(checkFecundation)
    {
        std::cout << "Testing that fecundation fuses two gametes...\n";
        Individual mom = Individual(arch, 0.0);
        Individual dad = Individual(arch, 1.0);
        Haplotype egg = mom.recombine(arch);
        Haplotype sperm = dad.recombine(arch);
        Individual baby = Individual(arch, egg, sperm);
        BOOST_CHECK_EQUAL(sumbool(baby.getSequence(0u)), 0u);
        BOOST_CHECK_EQUAL(sumbool(baby.getSequence(1u)), arch.nLoci);
    }

    BOOST_AUTO_TEST_CASE(checkNoMutation)
    {
        std::cout << "Testing absence of mutations...\n";
        Individual ind = Individual(arch, 0.0);
        Haplotype gamete = ind.recombine(arch);
        ind.mutate(gamete, 0.0);
        BOOST_CHECK_EQUAL(sumbool(gamete), 0u);
    }

    BOOST_AUTO_TEST_CASE(checkHighMutation)
    {
        std::cout << "Testing high mutation rate...\n";
        Individual ind = Individual(arch, 0.0);
        Haplotype gamete = ind.recombine(arch);
        ind.mutate(gamete, 100.0);
        BOOST_CHECK(sumbool(gamete) != 0u);
    }

    BOOST_AUTO_TEST_CASE(checkHighRecombination)
    {
        std::cout << "Testing high recombination rate...\n";
        Individual mom = Individual(arch, 0.0);
        Individual dad = Individual(arch, 1.0);
        Haplotype egg = mom.recombine(arch);
        Haplotype sperm = dad.recombine(arch);
        Individual baby = Individual(arch, egg, sperm);
        arch.recombinationRate = 10.0;
        Haplotype gamete = baby.recombine(arch);
        BOOST_CHECK(sumbool(gamete) != 0u);
        BOOST_CHECK(sumbool(gamete) != arch.nLoci);
    }

    BOOST_AUTO_TEST_CASE(checkExpression)
    {
        std::cout << "Testing gene expression...\n";
        Individual ind1 = Individual(arch, 0.0);
        std::vector<double> expression1 = ind1.getExpression();
        BOOST_CHECK_EQUAL(sum(expression1), -1.0 * arch.nLoci);
        Individual ind2 = Individual(arch, 1.0);
        std::vector<double> expression2 = ind2.getExpression();
        BOOST_CHECK_EQUAL(sum(expression2), arch.nLoci);
    }

BOOST_AUTO_TEST_SUITE_END()

// Check the absence of recombination
BOOST_AUTO_TEST_CASE(checkNoRecombination)
{
    std::cout << "Testing meiosis without recombination...\n";
    ParameterSet pars;
    pars.setNChromosomes(1u); // to avoid free recombination
    GeneticArchitecture arch = GeneticArchitecture(pars);
    MultiNet networks = arch.getNetworks();
    Individual mom = Individual(arch, 0.0);
    Individual dad = Individual(arch, 1.0);
    Haplotype egg = mom.recombine(arch);
    Haplotype sperm = dad.recombine(arch);
    Individual baby = Individual(arch, egg, sperm);
    Haplotype gam = baby.recombine(arch);
    BOOST_CHECK_EQUAL(sumbool(gam) / arch.nLoci, gam[0u]);
}

BOOST_AUTO_TEST_CASE(checkDevelopment)
{
    std::cout << "Testing developing individual...\n";
    ParameterSet pars;
    pars.setDominanceVariance(0.0);
    GeneticArchitecture arch = GeneticArchitecture(pars);
    MultiNet networks = arch.getNetworks();
    Individual mom = Individual(arch, 0.0);
    Individual dad = Individual(arch, 1.0);
    Haplotype egg = mom.recombine(arch);
    Haplotype sperm = dad.recombine(arch);
    Individual baby = Individual(arch, egg, sperm);
    BOOST_CHECK_EQUAL(baby.getEcoTrait(), 0.0);
    BOOST_CHECK_EQUAL(baby.getMatePref(), 0.0);
    BOOST_CHECK_EQUAL(baby.getNeutral(), 0.0);
}






