#include "MetaPop.h"
#include <iostream>
#include <fstream>


size_t MetaPop::evolve(const Genome &genome, const MultiNet &networks)
{
    size_t t = 0u;

    StreamBag outfiles;

    if (record) {

        outfiles.openAll();
        outfiles.checkAll();
    }

    for (; t < tmax; ++t) {

        if (record && t % tsave == 0u) {

            loadBuffer(t);
            buffer.write(outfiles.outTime, buffer.time);
            buffer.write(outfiles.outPopSize0, buffer.popsize0);
            buffer.write(outfiles.outPopSize1, buffer.popsize1);

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
        outfiles.closeAll();
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

void StreamBag::openAll()
{
    outTime.open("time.dat", std::ios::binary);
    outPopSize0.open("popsize0.dat", std::ios::binary);
    outPopSize1.open("popsize1.dat", std::ios::binary);
}

void checkOpen(std::ofstream &out, std::string var)
{
    if (!out.is_open()) {
        std::string msg = "Unable to open file " + var + ".dat";
        throw std::runtime_error(msg);
    }
}

void StreamBag::checkAll()
{
    checkOpen(outTime, "time");
    checkOpen(outPopSize0, "popsize0");
    checkOpen(outPopSize1, "popsize1");
}

void StreamBag::closeAll()
{
    outTime.close();
    outPopSize0.close();
    outPopSize1.close();
}
