#include "MetaPop.h"



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

            // Prepare analysis for output

            ecomean = 0.0; // global mean trait value

            for (size_t p = 0u; p < 2u; ++p) {

                // The means are used multiple times
                pops[p].calcMeanEcoTrait();
                pops[p].calcMeanMatePref();
                pops[p].calcMeanNtrTrait();

                ecomean += pops[p].getPopSize() * pops[p].getMeanEcoTrait();

            }

            // Assign ecotypes
            for (size_t p = 0u; p < 2u; ++p)
                pops[p].assignEcotypes(ecomean);

            // Load output to buffer
            loadBuffer(t);

            // Write to files
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
        pops[0u].reproduce(birth, sexsel, genome, networks);
        pops[1u].reproduce(birth, sexsel, genome, networks);

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
    ei -= sqr(ecomean); // minus the square of the mean
    return sqrt(ei); // return standard deviation
}

double MetaPop::getSpatialIsolation()
{

    // Calculated as some non-overlap between ecotypes

    // Ecotype-by-habitat table
    std::vector<vecUns> n = { uzeros(2u), uzeros(2u) };

    for (size_t p = 0u; p < 2u; ++p) {
        auto pop = pops[p];
        for (size_t i = 0u; i < pop.getPopSize(); ++i) {
            auto ind = pop.individuals[i];
            size_t ecotype = ind->getEcotype();
            ++n[p][ecotype];
        }
    }

    double si = n[0u][0u] * n[1u][1u] - n[0u][1u] * n[1u][0u];
    si /= sqrt(n[0u][0u] + n[0u][1u]);
    si /= sqrt(n[1u][0u] + n[1u][1u]);
    si /= sqrt(n[0u][0u] + n[1u][0u]);
    si /= sqrt(n[0u][1u] + n[1u][1u]);

    return si;
}


double MetaPop::getMatingIsolation()
{

    // Count homogamic and heterogamic crossings

    // For each of them sample a number of males to encounter, from the metapop
    // Evaluate each male by a yes or no
    // Update the table of matings accordingly by looking at ecotypes
    // RI = 0 is mating tests are not possible

    size_t nfemales = 0u;
    size_t nmales = 0u;

    for (size_t p = 0u; p < 2u; ++p) {
        nfemales += pops[p].getNFemales();
        nmales += pops[p].getNMales();
    }

    if (nfemales == 0u || nmales == 0u)
        return 0.0;

    // Table of crossings
    std::vector<vecUns> m = { uzeros(2u), uzeros(2u) };

    // Make a vector with all males of the metapop
    Crowd allMales;
    for (size_t p = 0u; p < 2u; ++p)
        for (size_t i = 0u; i < pops[p].getNMales(); ++i)
            allMales.push_back(pops[p].individuals[i]);

    // Loop through the females of the metapop and test their preference
    for (size_t p = 0u; p < 2u; ++p) {
        for (size_t i = 0u; i < pops[p].getNFemales(); ++i) {
            auto fem = pops[p].females[i];
            size_t nencounters = rnd::poisson(1.0 / matingcost);
            while (nencounters) {
                auto candidate = allMales[rnd::random(allMales.size())];
                if (fem->acceptMate(candidate->getEcoTrait(), sexsel))
                    ++m[fem->getEcotype()][candidate->getEcotype()];
                --nencounters;
            }
        }
    }

    double ri = m[0u][0u] * m[1u][1u] - m[0u][1u] * m[1u][0u];
    ri /= sqrt(m[0u][0u] * m[1u][1u] * m[0u][1u] * m[1u][0u]);
    return ri;

}


double size2dbl(const size_t &x)
{
    return static_cast<double>(x);
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
    buffer.add(getSpatialIsolation());
    buffer.add(getMatingIsolation());
}




