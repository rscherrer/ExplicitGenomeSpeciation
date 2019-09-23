#include "doMain.h"
#include "Param.h"
#include "GenArch.h"
#include "Random.h"
#include "Deme.h"
#include "MetaPop.h"
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <cassert>

int doMain(const vecStrings &args)
{

    try
    {

        if (args.size() > 2u)
            throw std::runtime_error("More than one argument were supplied");

        Param pars;

        // Read parameters from a file if supplied
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

        // Random number generator
        rnd::rng.seed(pars.getSeed());
        rnd::rng.seed(42u);

        GenArch arch = GenArch(pars);

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
