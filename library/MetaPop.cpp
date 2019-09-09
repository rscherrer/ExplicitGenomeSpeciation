#include "MetaPop.h"
#include <iostream>
#include <fstream>


size_t MetaPop::evolve(const Genome &genome, const MultiNet &networks)
{
    size_t t = 0u;

    std::ofstream output;

    if (record) {
        output.open("record.txt", std::ios::binary);
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

    /* // Turn the following into a test
    std::ifstream input;
    int size = 0;
    char elt0[25];
    char elt1[25];
    char elt2[25];
    input.open("record.txt", std::ios::binary);
    input.seekg(0, std::ios::end);
    size = (int) input.tellg();
    input.seekg(0, std::ios::beg);
    while(input.tellg() < size)
    {
        input.read((char *) elt0, sizeof(elt0));
        input.read((char *) elt1, sizeof(elt1));
        input.read((char *) elt2, sizeof(elt2));
    }
    input.close();
    std::cout << elt0 << '\t' << elt1 << '\t' << elt2 << '\n';
    */

    return t;
}


void MetaPop::loadBuffer(const size_t &t)
{
    sprintf(buffer.time, "%lu", t);
    sprintf(buffer.popsize0, "%lu", pops[0u].getPopSize());
    sprintf(buffer.popsize1, "%lu", pops[1u].getPopSize());
}

void Buffer::write(std::ofstream &out)
{
    out.write((char *) &time, sizeof(time));
    out.write((char *) &popsize0, sizeof(popsize0));
    out.write((char *) &popsize1, sizeof(popsize1));
}
