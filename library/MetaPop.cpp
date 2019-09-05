#include "MetaPop.h"
#include <iostream>

size_t MetaPop::evolve(const Genome &genome, const MultiNet &networks)
{
    size_t t = 0u;

    std::ofstream output;
    output.open("record.dat");
    if (!output.is_open())
        throw std::runtime_error("Unable to open file record.dat\n");

    for (; t < tmax; ++t) {

        output << "t = " << t << '\n';

        // Dispersal
        Crowd migrants1 = pops[0u].emigrate(dispersal);
        Crowd migrants2 = pops[1u].emigrate(dispersal);
        pops[0u].immigrate(migrants2);
        pops[1u].immigrate(migrants1);

        // Resource acquisition
        pops[0u].consume();
        pops[1u].consume();

        // Reproduction
        pops[0u].reproduce(birth, mating, genome, networks);
        pops[1u].reproduce(birth, mating, genome, networks);

        // Survival
        if (!pops[0u].survive(survival) && !pops[1u].survive(survival)) {
            std::cout << "The population went extinct at t = " << t << '\n';
            break;
        }
    }

    output.close();

    return t;
}
