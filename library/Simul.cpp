#include "Simul.h"

int simulate(const vecStrings &args)
{

    try
    {

        if (args.size() > 2u)
            throw std::runtime_error("More than one argument were supplied");

        // Create a default parameter set
        Param pars;

        // Read parameters from a file if supplied
        if (args.size() == 2) pars.read(args[1u]);

        // Random number generator
        rnd::rng.seed(pars.getSeed());

        // Create a genetic architecture
        GenArch arch = GenArch(pars);

        // Create a metapopulation with two demes, one of which is empty
        MetaPop metapop = MetaPop(pars, arch, true);

        // Evolve the metapopulation
        std::clog << "Simulation started\n";
        metapop.evolve(arch);
        std::clog << "Simulation ended\n";

    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
