#include "doMain.h"
#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "Random.h"
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cassert>


/// Function to run a simulation
void runSimulation(size_t &t, Population &pop, const size_t &tmax)

{
    // Loop through time
    // At every generation perform the steps of the life cycle
    // and update the population accordingly
    for (; t < tmax; ++t);

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
        Random rnd = Random(42u);

        // Create a genetic architecture
        GeneticArchitecture arch = GeneticArchitecture(pars, rnd);

        std::cout << "Architecture created\n";

        // Create a population of individuals
        Population pop = Population(pars.getInitialPopSize());

        // Run the simulation
        size_t t = 0u;
        size_t tmax = pars.getTEndSim();
        runSimulation(t, pop, tmax);

    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
