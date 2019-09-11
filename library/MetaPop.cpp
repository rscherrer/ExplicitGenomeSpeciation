#include "MetaPop.h"
#include "Individual.h"
#include "utils.h"
#include <iostream>
#include <fstream>


size_t MetaPop::evolve(const Genome &genome, const MultiNet &networks)
{
    size_t t = 0u;

    StreamBag out;

    if (record)
        out.openAll();

    for (; t < tmax; ++t) {

        // Sort out the sexes
        pops[0u].sortSexes();
        pops[1u].sortSexes();

        if (record && t % tsave == 0u) {

            for (size_t p = 0u; p < 2u; ++p) {
                // The means are used multiple times
                pops[p].calcMeanEcoTrait();
                pops[p].calcMeanMatePref();
                pops[p].calcMeanNtrTrait();
            }

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

    if (record)
        out.closeAll();

    return t;
}

double MetaPop::getEcoIsolation()
{
    // Ecological isolation is the standard deviation in ecological trait value

    double ei = 0.0;

    for (size_t p = 0u; p < 2u; ++p) {
        auto pop = pops[p];
        for (size_t i = 0u; i < pop.getPopSize(); ++i) {
            auto ind = pop.individuals[i];
            ei += sqr(ind->getEcoTrait());
        }
    }

    const size_t n0 = pops[0u].getPopSize();
    const size_t n1 = pops[1u].getPopSize();

    ei /= (n0 + n1);

    const double x0 = pops[0u].getMeanEcoTrait();
    const double x1 = pops[1u].getMeanEcoTrait();

    ei -= sqr((x0 * n0 + x1 * n1)); // minus the square of the mean

    return sqrt(ei); // return standard deviation
}


double size2dbl(const size_t &x)
{
    return static_cast<double>(x);
}

void Buffer::flush()
{
    fields.clear();
}

void Buffer::add(const double &x)
{
    fields.push_back(x);
}

void MetaPop::loadBuffer(const size_t &t)
{
    buffer.flush();
    buffer.add(size2dbl(t));
    buffer.add(size2dbl(pops[0u].getPopSize()));
    buffer.add(size2dbl(pops[1u].getPopSize()));
    buffer.add(size2dbl(pops[0u].getNFemales()));
    buffer.add(size2dbl(pops[1u].getNFemales()));
    buffer.add(pops[0u].getResources()[0u]);
    buffer.add(pops[0u].getResources()[1u]);
    buffer.add(pops[1u].getResources()[0u]);
    buffer.add(pops[1u].getResources()[1u]);
    buffer.add(pops[0u].getMeanEcoTrait());
    buffer.add(pops[1u].getMeanEcoTrait());
    buffer.add(pops[0u].getMeanMatePref());
    buffer.add(pops[1u].getMeanMatePref());
    buffer.add(pops[0u].getMeanNtrTrait());
    buffer.add(pops[1u].getMeanNtrTrait());
    buffer.add(getEcoIsolation());
    // buffer.add(getSpatialIsolation());
    // buffer.add(getMatingIsolation());
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
