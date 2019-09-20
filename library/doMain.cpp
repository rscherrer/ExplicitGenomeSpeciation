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
int doMain(const vecStrings &args)
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

        // Create a metapopulation
        const vecUns popsizes = { pars.getInitialPopSize(), 0u };
        MetaPop metapop = MetaPop(popsizes, pars, arch);

        // Launch simulation
        std::cout << "Simulation started\n";
        metapop.evolve(arch);
        std::cout << "Simulation ended\n";

    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
