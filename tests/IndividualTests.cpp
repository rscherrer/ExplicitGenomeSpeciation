#include "library/Individual.h"
#include "library/Utilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

// Test functions associated with individual

// Test genome generation

BOOST_AUTO_TEST_CASE(GenerateOnlyZeroAlleles)
{
    Param pars;
    pars.allfreq = 0.0;
    GenArch arch = GenArch(pars);
    Individual ind = Individual(pars, arch);
    BOOST_CHECK_EQUAL(ind.getAlleleSum(), 0u);
}

BOOST_AUTO_TEST_CASE(GenerateOnlyOneAlleles)
{
    Param pars;
    pars.allfreq = 1.0;
    GenArch arch = GenArch(pars);
    Individual ind = Individual(pars, arch);
    BOOST_CHECK_EQUAL(ind.getAlleleSum(), 2u * pars.nloci);
}

// Test mate choice

BOOST_AUTO_TEST_CASE(HomogamousFemaleAcceptsIdenticalMale)
{
    Param pars;
    pars.allfreq = 0.5;
    pars.sexsel = 1.0;
    GenArch arch = GenArch(pars);
    Individual ind = Individual(pars, arch);
    ind.resetEcoTrait(0.0, pars);
    ind.resetMatePref(1.0); // full assortative mating
    BOOST_CHECK_EQUAL(ind.mate(0.0, pars), 1.0);
}

BOOST_AUTO_TEST_CASE(HeterogamousFemaleRejectsIdenticalMale)
{
    Param pars;
    pars.allfreq = 0.5;
    pars.sexsel = 1.0;
    GenArch arch = GenArch(pars);
    Individual ind = Individual(pars, arch);
    ind.resetEcoTrait(0.0, pars);
    ind.resetMatePref(-1.0); // full disassortative mating
    BOOST_CHECK_EQUAL(ind.mate(0.0, pars), 0.0);
}

BOOST_AUTO_TEST_CASE(RandomMatingFemaleAcceptsAnyone)
{
    Param pars;
    pars.allfreq = 0.5;
    pars.sexsel = 1.0;
    GenArch arch = GenArch(pars);
    Individual ind = Individual(pars, arch);
    ind.resetEcoTrait(0.0, pars);
    ind.resetMatePref(0.0); // random mating
    BOOST_CHECK_EQUAL(ind.mate(0.0, pars), 1.0);
    BOOST_CHECK_EQUAL(ind.mate(1.0, pars), 1.0);
    BOOST_CHECK_EQUAL(ind.mate(-1.0, pars), 1.0);
}

// Test fecundation
// (this is also a test of no mutation)
BOOST_AUTO_TEST_CASE(MatingBetweenAlternativeHomozygotesGivesHeterozygoteNoMut)
{
    Param pars;
    pars.mutation = 0.0; // no mutation
    GenArch arch = GenArch(pars);
    pars.allfreq = 0.0;
    Individual mom = Individual(pars, arch);
    pars.allfreq = 1.0;
    Individual dad = Individual(pars, arch);
    Individual baby = Individual(pars, arch, mom, dad);
    BOOST_CHECK_EQUAL(baby.getAlleleSum(), pars.nloci);
}

// Test mutation
BOOST_AUTO_TEST_CASE(MutationAltersTheGenome)
{
    Param pars;
    pars.mutation = 100.0; // high mutation
    GenArch arch = GenArch(pars);
    pars.allfreq = 0.0;
    Individual mom = Individual(pars, arch);
    Individual dad = Individual(pars, arch);
    Individual baby = Individual(pars, arch, mom, dad);
    BOOST_CHECK(baby.getAlleleSum() > 0u);
}

// Test recombination

BOOST_AUTO_TEST_CASE(BackcrossPredictableIfNoRecombination)
{
    Param pars;
    pars.nchrom = 1u; // one chromosome to avoid free recombination
    pars.mutation = 0.0; // no mutation
    pars.recombination = 0.0; // no recombination
    GenArch arch = GenArch(pars);
    pars.allfreq = 0.0;
    Individual mom = Individual(pars, arch); // aa
    pars.allfreq = 1.0;
    Individual dad = Individual(pars, arch); // AA
    Individual f1 = Individual(pars, arch, mom, dad); // Aa
    Individual backcross = Individual(pars, arch, mom, f1);

    // Backcross should be either full aa or full Aa if no recombination
    const size_t allsum = backcross.getAlleleSum();
    BOOST_CHECK(allsum == 0u || allsum == pars.nloci);
}

BOOST_AUTO_TEST_CASE(BackcrossUnpredictableIfRecombination)
{
    Param pars;
    pars.mutation = 0.0; // no mutation
    pars.recombination = 10.0; // high recombination
    GenArch arch = GenArch(pars);
    pars.allfreq = 0.0;
    Individual mom = Individual(pars, arch); // aa
    pars.allfreq = 1.0;
    Individual dad = Individual(pars, arch); // AA
    Individual f1 = Individual(pars, arch, mom, dad); // Aa
    Individual backcross = Individual(pars, arch, mom, f1);

    // Backcross should be either full aa or full Aa if no recombination
    const size_t allsum = backcross.getAlleleSum();
    BOOST_CHECK(allsum > 0u && allsum != pars.nloci);
}

// Test gene expression

BOOST_AUTO_TEST_CASE(AlleleZeroInhibitsGeneExpression)
{
    Param pars;
    pars.allfreq = 0.0;
    GenArch arch = GenArch(pars);
    Individual ind = Individual(pars, arch);
    BOOST_CHECK_EQUAL(ind.getExpression(), -1.0 * pars.nloci);
}

BOOST_AUTO_TEST_CASE(AlleleOneEnhancesGeneExpression)
{
    Param pars;
    pars.allfreq = 1.0;
    GenArch arch = GenArch(pars);
    Individual ind = Individual(pars, arch);
    BOOST_CHECK_EQUAL(ind.getExpression(), pars.nloci);
}

// Test development
BOOST_AUTO_TEST_CASE(HybridDevelopsWithZeroTraitValuesIfCodominance)
{
    Param pars;
    pars.mutation = 0.0;
    pars.dominancevar = 0.0; // codominance
    GenArch arch = GenArch(pars);
    pars.allfreq = 0.0;
    Individual mom = Individual(pars, arch);
    pars.allfreq = 1.0;
    Individual dad = Individual(pars, arch);
    Individual baby = Individual(pars, arch, mom, dad);

    // Without dominance the expression of each trait should be zero
    BOOST_CHECK_EQUAL(baby.getEcoTrait(), 0.0);
    BOOST_CHECK_EQUAL(baby.getMatePref(), 0.0);
    BOOST_CHECK_EQUAL(baby.getNeutral(), 0.0);
}







