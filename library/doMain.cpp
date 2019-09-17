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
            throw std::runtime_error("More than one argument were supplied");

        // Create a default parameter set
        ParameterSet pars;

        // Now is the time to update parameters if some are provided
        if (args.size() == 2u) {

            std::string filename = args[1u];
            std::ifstream inputfile;
            inputfile.open(filename);
            if (!inputfile.is_open()) {
                std::string msg = "Unable to open parameter file ";
                throw std::runtime_error(msg + filename);
            }

            pars.readParams(inputfile);
            inputfile.close();

        }

        // Create and seed a random number generator
        rnd::rng.seed(pars.getSeed());

        // Create a genetic architecture
        GeneticArchitecture arch = GeneticArchitecture(pars);

        // Create a population of individuals
        Genome genome = arch.getGenome();
        MultiNet networks = arch.getNetworks();

        // Create populations
        const size_t n0 = pars.getInitialPopSize();
        const double foodmax = pars.getMaxResourceCapacity();
        const double foodgrowth = pars.getMaxResourceGrowth();
        const double symmetry = pars.getHabitatSymmetry();
        const vecDbl foodmax1 = {foodmax, symmetry * foodmax};
        const vecDbl foodmax2 = {symmetry * foodmax, foodmax};
        const vecDbl foodgrows = {foodgrowth, foodgrowth};

        Population pop1 = Population(n0, genome, networks, foodmax1, foodgrows);
        Population pop2 = Population(0u, genome, networks, foodmax2, foodgrows);
        vecPop pops = {pop1, pop2};

        // Create a metapopulation
        MetaPop metapop = MetaPop(pops, pars);

        // Open a data file
        std::ofstream out;
        out.open("output.dat");
        if (!out.is_open())
            throw std::runtime_error("Unable to open output file");

        // Launch simulation
        std::cout << "Simulation started\n";
        metapop.evolve(genome, networks);
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
