#include "Simul.h"

bool timetosave(const int &t,const Param &p)
{
    return p.record && t >= 0 && t % p.tsave == 0;
}

bool timetofreeze(const int &t,const Param &p)
{
    return p.record && p.gensave && t >= 0 && t % p.tfreeze == 0;
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

        // Create a genetic architecture
        GenArch arch = GenArch(pars);

        // Load the genetic architecture if necessary
        if (pars.archload) arch.load(pars);

        // Save parameters if necessary
        if (pars.parsave) pars.save();

        // Create a metapopulation with two demes
        MetaPop metapop = MetaPop(pars, arch);

        // Create an analytical module
        Collector collector = Collector(arch);

        std::cout << "Simulation started.\n";

        // Loop through time
        for (int t = -pars.tburnin; t < pars.tend; ++t) {

            if (t == 0) metapop.exitburnin();

            if (pars.talkative) std::clog << t << '\n';

            // Life cycle of the metapopulation
            metapop.disperse(pars);
            metapop.consume(pars);

            // Analyze the metapopulation if needed
            if (timetosave(t, pars)) {
                collector.analyze(metapop, pars, arch); // collect stats
                if (pars.datsave) collector.print(t, metapop); // save them to files
            }

            // Save whole genomes if needed (space-consuming)
            if (timetofreeze(t, pars)) {
                collector.freeze(metapop, pars);
            }

            metapop.reproduce(pars, arch);
            metapop.survive(pars);

            // Is the population still there?
            if (metapop.isextinct()) {
                std::cout << "The population went extinct at t = " << t << '\n';
                break;
            }
        }

        std::cout << "Simulation ended.\n";
        return 0;
    }
    catch (const std::exception& err)
    {
        std::cerr << "Exception: " << err.what() << '\n';
    }
    catch (const char* err)
    {
        std::cerr << "Exception: " << err << '\n';
    }
    catch (...)
    {
        std::cerr << "Unknown Exception\n";
    }
    return 1;
}
