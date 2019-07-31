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
void runSimulation(size_t &t, Population &pop, const size_t &tmax,
 const double &survival, const double &birth, const double &strength)
{
    // Loop through time...
    for (; t < tmax; ++t) {

        // Dispersion

        // Resource acquisition

        // Reproduction
        pop.reproduce(birth, strength);

        // Survival
        if (!pop.survive(survival)) {
            std::cout << "The population went extinct at t = " << t << '\n';
            break;
        }

    }

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

        std::cout << "Simulation started\n";

        // Run the simulation
        size_t t = 0u;
        size_t tmax = pars.getTEndSim();
        const double survival = pars.getSurvivalProb();
        const double birth = pars.getBirthRate();
        const double strength = pars.getMatePreferenceStrength();

        runSimulation(t, pop, tmax, survival, birth, strength);

        std::cout << "Simulation ended\n";

    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
