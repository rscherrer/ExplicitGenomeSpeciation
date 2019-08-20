#include "library/Genome.h"
#include "library/utils.h"
#include "tests/GenFixture.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_CASE(checkChromosomes)
{
    std::cout << "Testing chromosome lengths...\n";
    ParameterSet pars;
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    BOOST_CHECK_EQUAL(genome.chromosomes[0u], 1.0 / 3.0);
    BOOST_CHECK_EQUAL(genome.chromosomes[1u], 2.0 / 3.0);
    BOOST_CHECK_EQUAL(genome.chromosomes[2u], 1.0);
}

BOOST_AUTO_TEST_CASE(checkEncodedTraits)
{
    std::cout << "Testing a genome encoding just one trait...\n";
    ParameterSet pars;
    pars.setNLociPerTrait({10u, 2u, 2u});
    GeneticArchitecture arch = GeneticArchitecture(pars);
    Genome genome = arch.getGenome();
    BOOST_CHECK(sumu(genome.traits) == 6u);
}
