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

        // Open a data file
        std::ofstream out;
        out.open("output.dat");
        if (!out.is_open())
            throw std::runtime_error("Unable to open output file");

        // Launch simulation
        std::cout << "Simulation started\n";
        meta.evolve(genome, networks);
        std::cout << "Simulation ended\n";

        out.close();

    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
