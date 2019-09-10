#include "MetaPop.h"
#include <iostream>
#include <fstream>


size_t MetaPop::evolve(const Genome &genome, const MultiNet &networks)
{
    size_t t = 0u;

    std::ofstream output;

    if (record) {
        output.open("record.dat", std::ios::binary);
        if (!output.is_open())
            throw std::runtime_error("Unable to open file record.dat\n");
    }

    for (; t < tmax; ++t) {

        if (record && t % tsave == 0u) {

            loadBuffer(t);
            buffer.write(output);

        }

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

    if (record) output.close();

    return t;
}


void MetaPop::loadBuffer(const size_t &t)
{
    buffer.time = static_cast<double>(t);
    buffer.popsize0 = static_cast<double>(pops[0u].getPopSize());
    buffer.popsize1 = static_cast<double>(pops[1u].getPopSize());
}

void Buffer::write(std::ofstream &out)
{
    out.write((char *) &time, sizeof(time));
    out.write((char *) &popsize0, sizeof(popsize0));
    out.write((char *) &popsize1, sizeof(popsize1));
}
