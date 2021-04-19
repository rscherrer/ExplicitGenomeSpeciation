#include "library/MetaPop.h"
#include "library/Utilities.h"
#include <boost/test/unit_test.hpp>

  // Tests of the good behavior of a metapopulation object

  // Simulation should reach tmax in the absence of mortality
BOOST_AUTO_TEST_CASE(ConstantPopSizeWhenNoBirthNoDeath)
{
  // std::clog << "Testing constant population size...\n";
  Param pars;
  pars.tend = 10;
  pars.tburnin = 0;
  pars.demesizes = { 10u, 10u };
  pars.dispersal = 0.0;
  pars.survival = 1.0; // no death
  pars.birth = 0.0; // no birth
  pars.sexsel = 0.0;

  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  for (int t = 0; t < pars.tend; ++t)
    metapop.cycle(pars, arch);

  BOOST_CHECK_EQUAL(metapop.getSize(), 20u);
}

BOOST_AUTO_TEST_CASE(InstantExtinctionIfSurvivalIsZero)
{
  // std::clog << "Testing instantaneous extinction...\n";
  Param pars;
  pars.tend = 100;
  pars.tburnin = 0;
  pars.demesizes = { 10u, 10u };
  pars.dispersal = 0.0;
  pars.survival = 0.0; // 0% chance of survival
  pars.birth = 0.0; // no birth
  pars.sexsel = 0.0;

  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  metapop.cycle(pars, arch);
  BOOST_CHECK(metapop.isextinct());
}

  // Simulation should end prematurely with high mortality
BOOST_AUTO_TEST_CASE(ProgressiveExtinctionWhenLowSurvival)
{
  // std::clog << "Testing progressive extinction...\n";
  Param pars;
  pars.tend = 100;
  pars.tburnin = 0;
  pars.demesizes = { 10u, 10u };
  pars.dispersal = 0.0;
  pars.survival = 0.1; // 10% chance of survival
  pars.birth = 0.0; // no birth
  pars.sexsel = 0.0;

  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);

  bool extinction = false;

  for (int t = 0; t < pars.tend; ++t) {
    metapop.cycle(pars, arch);
    if (metapop.isextinct()) {
      extinction = true;
      break;
    }
  }

  BOOST_CHECK(extinction);
}

  // Test that habitats are initialized with only one resource if
  // habitat symmetry is zero
BOOST_AUTO_TEST_CASE(HabitatsHaveOneResourceIfCompleteAsymmetry)
{
  // std::clog << "Testing complete habitat asymmetry...\n";
  Param pars;
  GenArch arch = GenArch(pars);
  pars.rdynamics = 1u;
  pars.inflow = 400.0;
  pars.outflow = 100.0;
  pars.hsymmetry = 0.0; // full habitat asymmetry
  pars.demesizes = {10u, 10u};
  pars.ecosel = 0.0;
  pars.tburnin = 0;
  pars.dispersal = 0.0;
  MetaPop metapop = MetaPop(pars, arch);
  metapop.cycle(pars, arch);
  BOOST_CHECK_EQUAL(metapop.getResource(0u, 0u), 400.0/110.0);
  BOOST_CHECK_EQUAL(metapop.getResource(1u, 1u), 400.0/110.0);
  BOOST_CHECK_EQUAL(metapop.getResource(0u, 1u), 0.0);
  BOOST_CHECK_EQUAL(metapop.getResource(1u, 0u), 0.0);

}

BOOST_AUTO_TEST_CASE(NoDispersalLeavesHabitatsWithSameNumberOfIndividuals)
{
  // std::clog << "Testing no dispersal...\n";
  Param pars;
  pars.demesizes = { 15u, 10u };
  pars.dispersal = 0.0; // no dispersal
  pars.survival = 1.0; // no death
  pars.birth = 0.0; // no birth
  pars.tburnin = 0;
  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  metapop.cycle(pars, arch);
  BOOST_CHECK_EQUAL(metapop.getDemeSize(0u), 15u);
  BOOST_CHECK_EQUAL(metapop.getDemeSize(1u), 10u);
}

BOOST_AUTO_TEST_CASE(AllIndividualsMigrateIfDispersalIsMax)
{
  // std::clog << "Testing mass exodus...\n";
  Param pars;
  pars.demesizes = { 10u, 0u };
  pars.dispersal = 1.0; // 100% chance dispersal
  pars.survival = 1.0; // no death
  pars.birth = 0.0; // no birth
  pars.tburnin = 0;
  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  metapop.cycle(pars, arch);
  BOOST_CHECK_EQUAL(metapop.getDemeSize(0u), 0u);
  BOOST_CHECK_EQUAL(metapop.getDemeSize(1u), 10u);
}

  // Check that a population has grown after reproduction
BOOST_AUTO_TEST_CASE(ReproductionHasProducedNewIndividuals)
{
  // std::clog << "Testing reproduction...\n";
  Param pars;
  pars.capacity = 100.0;
  pars.dispersal = 0.0;
  pars.birth = 1000.0; // high birth rate
  pars.demesizes = { 10u, 0u };
  pars.survival = 1.0; // 100% chance survival
  pars.tburnin = 0u;
  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  metapop.cycle(pars, arch);
  BOOST_CHECK(metapop.getSize() >= 10u);
  BOOST_CHECK(metapop.getDemeSize(0u) >= 10u);
  BOOST_CHECK_EQUAL(metapop.getDemeSize(1u), 0u);
}


  // Newborns should not die (this one is annoying)

BOOST_AUTO_TEST_CASE(PopulationWipeOutLeavesOnlyNewborns)
{
  // std::clog << "Testing that newborns do not die...\n";
  Param pars;
  pars.birth = 100.0; // relatively high birth rate
  pars.demesizes = { 10u, 0u };
  pars.survival = 0.0; // all adults should die
  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  metapop.cycle(pars, arch);
  BOOST_CHECK(metapop.getSize() > 0u);
  BOOST_CHECK(metapop.getDemeSize(0u) > 0u);
  BOOST_CHECK_EQUAL(metapop.getDemeSize(1u), 0u);
}


  // After feeding, the resources should be depleted
BOOST_AUTO_TEST_CASE(ResourceIsDepletedAfterConsumption)
{
  // std::clog << "Testing resource depletion...\n";
  Param pars;
  pars.hsymmetry = 1.0;
  pars.capacity = 10.0;
  pars.demesizes = { 5u, 5u };
  pars.dispersal = 0.0; // everybody feeds in one habitat only
  pars.tburnin = 0;
  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  metapop.cycle(pars, arch);
  BOOST_CHECK(metapop.getResource(0u, 0u) < 10.0);
  BOOST_CHECK(metapop.getResource(0u, 1u) < 10.0);
  BOOST_CHECK(metapop.getResource(1u, 0u) < 10.0);
  BOOST_CHECK(metapop.getResource(1u, 1u) < 10.0);
}


  // Test logistic resource dynamics
BOOST_AUTO_TEST_CASE(KnownLogisticResourceEquilibrium)
{
  // std::clog << "Testing logistic resource equilibrium...\n";
  Param pars;
  pars.rdynamics = 0u;
  pars.dispersal = 0.0;
  pars.birth = 0.0;
  pars.survival = 1.0;
  pars.capacity = 1.0;
  pars.hsymmetry = 0.0;
  pars.replenish = 2000.0;
  pars.demesizes = { 10u, 0u };
  pars.ecosel = 1.0;
  pars.tburnin = 0;
  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  metapop.resetTraits(0u, -1.0, pars); // optimally adapted individuals
  metapop.cycle(pars, arch);

  // Predict resource equilibrium after consumption
  const double R0 = utl::round(199.0/200.0, 4u);
  const double R1 = 0.0;

  // Fitness should sum up to the amount of food consumed
  const double sumw = utl::round(10.0 * R0, 2u);

  BOOST_CHECK_EQUAL(utl::round(metapop.getResource(0u, 0u), 4u), R0);
  BOOST_CHECK_EQUAL(utl::round(metapop.getResource(0u, 1u), 4u), R1);
  BOOST_CHECK_EQUAL(utl::round(metapop.getSumFitness(), 2u), sumw);
  BOOST_CHECK_EQUAL(utl::round(metapop.getVarFitness(), 4u), 0.0);
}

  // Test logistic resource dynamics
BOOST_AUTO_TEST_CASE(KnownChemostatResourceEquilibrium)
{
  // std::clog << "Testing chemostat resource equilibrium...\n";
  Param pars;
  pars.rdynamics = 1u;
  pars.inflow = 400.0;
  pars.outflow = 100.0;
  pars.dispersal = 0.0;
  pars.birth = 0.0;
  pars.survival = 1.0;
  pars.hsymmetry = 0.0;
  pars.demesizes = { 10u, 0u };
  pars.ecosel = 1.0;
  pars.tburnin = 0;
  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  metapop.resetTraits(0u, -1.0, pars); // optimally adapted individuals
  metapop.cycle(pars, arch);

    // Predict resource equilibrium after consumption
  const double R0 = utl::round(400.0/110.0, 4u);
  const double R1 = 0.0;

    // Fitness should sum up to the amount of food consumed
  const double sumw = utl::round(10.0 * R0, 2u);

  BOOST_CHECK_EQUAL(utl::round(metapop.getResource(0u, 0u), 4u), R0);
  BOOST_CHECK_EQUAL(utl::round(metapop.getResource(0u, 1u), 4u), R1);
  BOOST_CHECK_EQUAL(utl::round(metapop.getSumFitness(), 2u), sumw);
  BOOST_CHECK_EQUAL(utl::round(metapop.getVarFitness(), 4u), 0.0);
}

BOOST_AUTO_TEST_CASE(EcotypeClassification)
{
  // std::clog << "Testing ecotype classification...\n";

  Param pars;
  pars.rdynamics = 1u;
  pars.dispersal = 0.0;
  pars.birth = 0.0;
  pars.survival = 1.0;
  pars.hsymmetry = 1.0;
  pars.demesizes = { 100u, 0u };
  pars.allfreq = 0.5;
  pars.ecosel = 1.0;
  pars.tburnin = 0u;

  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  metapop.cycle(pars, arch);

    // Test that ecological trait values in ecotype 1 are always higher than
    // in ecotype 0

  double xmax0 = -1.0;
  double xmin1 = 1.0;

    // For each individual
  for (size_t i = 0u; i < metapop.getSize(); ++i) {

      // Record trait value and ecotype
    const double x = metapop.getTrait(i, 0u);
    const size_t e = metapop.getEcotype(i);

    if (e == 0u && x > xmax0) xmax0 = x;
    if (e == 1u && x < xmin1) xmin1 = x;

  }

  BOOST_CHECK(xmax0 < xmin1);

}

BOOST_AUTO_TEST_CASE(KnownInitialGeneticDiversityFromAlleleFreq)
{
    // std::clog << "Testing known initial genetic diversity...\n";
    Param pars = Param();
    const GenArch arch = GenArch(pars);
    pars.allfreq = 0.0;
    MetaPop pop0 = MetaPop(pars, arch);
    pars.allfreq = 1.0;
    MetaPop pop1 = MetaPop(pars, arch);
    BOOST_CHECK_EQUAL(pop0.getAlleleSum(), 0u);
    BOOST_CHECK_EQUAL(pop1.getAlleleSum(), pop1.getSize() * 2u * pars.nloci);
    // pick random individual and random locus
    BOOST_CHECK_EQUAL(pop0.getZygosity(10u, 1u, pars.nloci), 0u);
    BOOST_CHECK_EQUAL(pop1.getZygosity(10u, 1u, pars.nloci), 2u);

}
