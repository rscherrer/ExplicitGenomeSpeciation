#include "MetaPop.h"

void MetaPop::resetEcoTraits(const size_t &p, const double &x)
{
    pops[p]->resetEcoTraits(x, ecosel, maxfeed);
}

void MetaPop::resetMatePrefs(const size_t &p, const double &y)
{
    pops[p]->resetMatePrefs(y);
}

void MetaPop::resetEcotypes(const size_t &p, const size_t &e)
{
    pops[p]->resetEcotypes(e);
}

void MetaPop::resetGenders(const size_t &p, const bool &sex)
{
    pops[p]->resetGenders(sex);
    pops[p]->sortSexes();
}

void MetaPop::analyze(const GenArch& arch)
{
    stats.reset(t, arch);
    stats.analyze(pops, arch);
    stats.setEcoIsolation();
    stats.setSpatialIsolation(pops);
    stats.setMatingIsolation(pops, matingcost, sexsel);
}

void MetaPop::consume()
{
    for (size_t p = 0u; p < 2u; ++p)
        pops[p]->consume();
}

int MetaPop::evolve(const GenArch &arch)
{
    bool isBurnin = true;
    t = - tburnin;

    Output out;
    if (record) out.openAll();

    for (; t < tmax; ++t) {

        std::clog << "t = " << t << '\n';

        // Sort out the sexes
        pops[0u]->sortSexes();
        pops[1u]->sortSexes();

        if (t > 0 && isBurnin) {
            isBurnin = false;
            for (auto pop : pops) pop->exitBurnIn();
        }

        // Analyze and record
        if (record && t % tsave == 0u && t > 0) {
            analyze(arch);
            stats.save(out);
        }

        // Dispersal (only if not burnin)
        if (t > 0) {
            Crowd migrants1 = pops[0u]->emigrate(dispersal);
            Crowd migrants2 = pops[1u]->emigrate(dispersal);
            pops[0u]->immigrate(migrants2);
            pops[1u]->immigrate(migrants1);
        }

        size_t isExtant = 0u;

        for (auto pop : pops) {
            pop->consume();
            pop->reproduce(birth, sexsel, matingcost, ecosel, maxfeed, arch);
            isExtant += pop->survive(survival);
        }

        if (isExtant == 0u) {
            std::clog << "The population went extinct at t = " << t << '\n';
            break;
        }
    }

    if (record) out.closeAll();

    return t;
}

size_t MetaPop::getNOffspring(const size_t &p) const
{
    return pops[p]->getNOffspring();
}

double MetaPop::getResource(const size_t &p, const size_t &r) const
{
    return pops[p]->getResource(r);
}
