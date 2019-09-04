#include "doMain.h"
#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "Random.h"
#include "Population.h"
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cassert>


/// Function to run a simulation
size_t runSimulation(MetaPop &pops, const size_t &tmax, const size_t &tsave,
 const double &dispersal, const double &survival, const double &birth,
  const double &mating, const Genome &genome, const MultiNet &networks,
   const bool &record)
{
    size_t t = 0u;

    // Loop through time...
    for (; t < tmax; ++t) {

        // Time to save?
        if (record && t % tsave == 0u) {
            std::cout << "t = " << t << '\n';
        }

        // Dispersal
        Crowd migrants1 = pops[0u].emigrate(dispersal);
        Crowd migrants2 = pops[1u].emigrate(dispersal);
        pops[0u].immigrate(migrants2);
        pops[1u].immigrate(migrants1);


        // Resource acquisition
        pops[0u].consume();
        pops[1u].consume();

        // Reproduction
        pops[0u].reproduce(birth, mating, genome, networks);
        pops[1u].reproduce(birth, mating, genome, networks);

        // Survival
        if (!pops[0u].survive(survival) && !pops[1u].survive(survival)) {
            std::cout << "The population went extinct at t = " << t << '\n';
            break;
        }

    }

    return t;

}


/// Program to run the main function
int doMain(const vecStr &args)
{

    try
    {
        // Return an error if there are more than one argument
        if (args.size() > 2u)
            throw std::runtime_error("More than one argument was supplied");

        // Create a default parameter set
        ParameterSet pars;

        // Create and seed a random number generator
        rnd::rng.seed(pars.getSeed());

        // Create a genetic architecture
        GeneticArchitecture arch = GeneticArchitecture(pars);

        // Create a population of individuals
        Genome genome = arch.getGenome();
        MultiNet networks = arch.getNetworks();
        const size_t n0 = pars.getInitialPopSize();

        const double foodmax = pars.getMaxResourceCapacity();
        const double foodgrowth = pars.getMaxResourceGrowth();

        // Symmetry in resource partitioning
        const double symmetry = pars.getHabitatSymmetry();
        const vecDbl foodmax1 = {foodmax, symmetry * foodmax};
        const vecDbl foodmax2 = {symmetry * foodmax, foodmax};
        const vecDbl foodgrows = {foodgrowth, foodgrowth};

        Population pop1 = Population(n0, genome, networks, foodmax1, foodgrows);
        Population pop2 = Population(n0, genome, networks, foodmax2, foodgrows);
        MetaPop metapop = {pop1, pop2};

        std::cout << "Simulation started\n";

        // Run the simulation
        size_t tmax = pars.getTEndSim();
        size_t tsave = pars.getTSave();
        const double dispersal = pars.getDispersalRate();
        const double survival = pars.getSurvivalProb();
        const double birth = pars.getBirthRate();
        const double mating = pars.getMatePreferenceStrength();
        const bool record = pars.getRecord();

        runSimulation(metapop, tmax, tsave, dispersal, survival, birth, mating,
         genome, networks, record);

        std::cout << "Simulation ended\n";

    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
