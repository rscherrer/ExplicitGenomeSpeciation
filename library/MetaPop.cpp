#include "MetaPop.h"
#include <iostream>
#include <fstream>


size_t MetaPop::evolve(const Genome &genome, const MultiNet &networks)
{
    size_t t = 0u;

    StreamBag out;

    if (record) {

        out.openAll();
    }

    for (; t < tmax; ++t) {

        if (record && t % tsave == 0u) {

            loadBuffer(t);

            for (size_t f = 0u; f < out.names.size(); ++f)
                buffer.write(out.files[f], buffer.fields[f]);

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
        out.closeAll();
    }

    return t;
}


void MetaPop::loadBuffer(const size_t &t)
{
    buffer.fields.clear();
    buffer.fields.push_back(static_cast<double>(t));
    buffer.fields.push_back(static_cast<double>(pops[0u].getPopSize()));
    buffer.fields.push_back(static_cast<double>(pops[1u].getPopSize()));
    buffer.fields.push_back(static_cast<double>(pops[0u].getResources()[0u]));
    buffer.fields.push_back(static_cast<double>(pops[0u].getResources()[1u]));
    buffer.fields.push_back(static_cast<double>(pops[1u].getResources()[0u]));
    buffer.fields.push_back(static_cast<double>(pops[1u].getResources()[1u]));
}

void Buffer::write(std::ofstream * &out, const double &value)
{
    out->write((char *) &value, sizeof(value));
}

void StreamBag::openAll()
{
    for (size_t f = 0u; f < names.size(); ++f) {
        std::string filename = names[f] + ".dat";
        files.push_back(new std::ofstream());
        files.back()->open(filename, std::ios::binary);
        if (!files.back()->is_open()) {
            std::string msg = "Unable to open output file " + filename;
            throw std::runtime_error(msg);
        }
    }
}

void StreamBag::closeAll()
{
    for (size_t f = 0u; f < names.size(); ++f) {
        files[f]->close();
    }
}
