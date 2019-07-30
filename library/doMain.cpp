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
void runSimulation(size_t &t, Population &pop, const size_t &tmax)

{
    // Loop through time...
    for (; t < tmax; ++t) {

        // Dispersion

        // Resource acquisition

        // Reproduction

        // Survival
        pop.survive(1.0);

    }

    pop.getPopSize();


}


/// Program to run the main function
int doMain(const std::vector<std::string> &args)
{

    try
    {
        // Return an error if there are more than one argument
        if (args.size() > 2u)
            throw std::runtime_error("More than one argument was supplied");

        // Create a default parameter set
        ParameterSet pars;

        // Create and seed a random number generator
        rnd::rng.seed(42u);

        // Create a genetic architecture
        GeneticArchitecture arch = GeneticArchitecture(pars);

        // Create a population of individuals
        Population pop = Population(pars.getInitialPopSize());

        std::cout << "Simulation startes\n";

        // Run the simulation
        size_t t = 0u;
        size_t tmax = pars.getTEndSim();
        runSimulation(t, pop, tmax);

        std::cout << "Simulation ended\n";

    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
