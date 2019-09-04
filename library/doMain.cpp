#include "doMain.h"
#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "Random.h"
#include "Population.h"
#include "MetaPop.h"
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cassert>


/// Function to run a simulation
size_t runSimulation(vecPop &pops, const size_t &tmax, const size_t &tsave,
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

        // Create populations
        Population pop1 = Population(n0, genome, networks, foodmax1, foodgrows);
        Population pop2 = Population(n0, genome, networks, foodmax2, foodgrows);
        vecPop metapop = {pop1, pop2};

        // Create a metapopulation
        MetaPop meta = MetaPop(metapop, pars);

        // Launch simulation
        std::cout << "Simulation started\n";
        meta.evolve(genome, networks);
        std::cout << "Simulation ended\n";

        // Open a data file
        std::ofstream out;
        out.open("output.dat");
        if (!out.is_open())
            throw std::runtime_error("Unable to open output file");

        out.close();

    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
