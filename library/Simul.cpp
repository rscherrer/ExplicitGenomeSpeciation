#include "Simul.h"

bool timetosave(const int &t,const Param &p)
{
    return (t > 0 && p.record && t % p.tsave == 0u);
}

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
        rnd::rng.seed(pars.seed);

        // Create a genetic architecture
        GenArch arch = GenArch(pars);

        // Create a metapopulation with two demes
        MetaPop metapop = MetaPop(pars, arch);

        // Create an analytical module
        Collector collector = Collector(arch);

        // Loop through time
        for (int t = -pars.tburnin; t < pars.tend; ++t) {

            if (t == 0) metapop.exitburnin();

            // Life cycle of the metapopulation
            metapop.cycle(pars, arch);

            // Is the population still there?
            if (metapop.isextinct()) break;

            // Analyze the metapopulation if needed
            if (timetosave(t, pars)) collector.analyze(metapop, pars);

        }
    }
    catch (const std::runtime_error &err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
        return 1;
    }

    return 0;
}
