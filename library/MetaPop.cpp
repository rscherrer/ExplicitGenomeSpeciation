#include "MetaPop.h"
#include <iostream>
#include <fstream>


size_t MetaPop::evolve(const Genome &genome, const MultiNet &networks)
{
    size_t t = 0u;

    std::ofstream outTime;
    std::ofstream outPopSize0;
    std::ofstream outPopSize1;

    // There should not be one but many records, one per variable
    if (record) {

        outTime.open("time.dat", std::ios::binary);
        outPopSize0.open("popsize0.dat", std::ios::binary);
        outPopSize1.open("popsize1.dat", std::ios::binary);

        if (!outPopSize0.is_open())
            throw std::runtime_error("Unable to open file popsize0.dat\n");

        if (!outPopSize1.is_open())
            throw std::runtime_error("Unable to open file popsize1.dat\n");

        if (!outTime.is_open())
            throw std::runtime_error("Unable to open file time.dat\n");
    }

    for (; t < tmax; ++t) {

        if (record && t % tsave == 0u) {

            loadBuffer(t);
            buffer.write(outTime, buffer.time);
            buffer.write(outPopSize0, buffer.popsize0);
            buffer.write(outPopSize1, buffer.popsize1);

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

    if (record) {
        outTime.close();
        outPopSize0.close();
        outPopSize1.close();
    }

    return t;
}


void MetaPop::loadBuffer(const size_t &t)
{
    buffer.time = static_cast<double>(t);
    buffer.popsize0 = static_cast<double>(pops[0u].getPopSize());
    buffer.popsize1 = static_cast<double>(pops[1u].getPopSize());
}

void Buffer::write(std::ofstream &out, const double &value)
{
    out.write((char *) &value, sizeof(value));
}
