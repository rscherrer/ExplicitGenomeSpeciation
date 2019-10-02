#include "MetaPop.h"


// Initialization

vecPop MetaPop::makeDemes(const GenArch &arch, const bool &isburnin)
{
    // Create the demes
    vecPop demes;
    for (size_t p = 0u; p < 2u; ++p) {
        for (size_t res = 0u; res < 2u; ++res) {
            resources[p][res] = maxresources;
            if (p != res) resources[p][res] *= symmetry;
            replenish[p][res] = maxreplenish;
        }
        const size_t n = popsizes[p];
        const double max = maxfeed;
        const vecDbl k = resources[p];
        const vecDbl r = replenish[p];
        demes.push_back(Deme(n, ecosel, max, k, r, arch, isburnin));
    }
    return demes;
}


// Main function

int MetaPop::evolve(const GenArch &arch)
{
    bool isBurnin = true;
    t = - tburnin;

    Output out;
    if (record) out.openAll();

    for (; t < tmax; ++t) {

        std::clog << "t = " << t << '\n';

        // Sort out the sexes
        pops[0u].sortSexes();
        pops[1u].sortSexes();

        if (t > 0 && isBurnin) {
            isBurnin = false;
            for (auto &pop : pops) pop.exitBurnIn();
        }        

        // Dispersal (only if not burnin)
        if (t > 0) {
            Crowd migrants1 = pops[0u].emigrate(dispersal);
            Crowd migrants2 = pops[1u].emigrate(dispersal);
            pops[0u].immigrate(migrants2);
            pops[1u].immigrate(migrants1);
        }

        size_t isExtant = 0u;

        for (auto &pop : pops) {
            pop.consume();
            pop.reproduce(birth, sexsel, matingcost, ecosel, maxfeed, arch);
            isExtant += pop.survive(survival);
        }

        if (isExtant == 0u) {
            std::clog << "The population went extinct at t = " << t << '\n';
            break;
        }

        // Analyze and record
        if (record && t % tsave == 0u && t > 0) {
            analyze(arch);
            stats.save(out);
        }
    }

    if (record) out.closeAll();

    return t;
}

void MetaPop::analyze(const GenArch& arch)
{
    stats.reset(t, arch);
    stats.analyze(pops, arch);
    stats.setEcoIsolation();
    stats.setSpatialIsolation(pops);
    stats.setMatingIsolation(pops, matingcost, sexsel);
}


// Getters

size_t MetaPop::getNOffspring(const size_t &p) const
{
    return pops[p].getNOffspring();
}

size_t MetaPop::getSumEcotypes(const size_t &p) const
{
    return pops[p].getSumEcotypes();
}

double MetaPop::getResource(const size_t &p, const size_t &r) const
{
    return pops[p].getResource(r);
}

double MetaPop::getVarP(const size_t &trait, const size_t &eco) const
{
    return stats.getVarP(trait, eco);
}

double MetaPop::getSsqPhe(const size_t &trait, const size_t &eco) const
{
    return stats.getSsqPhe(trait, eco);
}

double MetaPop::getSumPhe(const size_t &trait, const size_t &eco) const
{
    return stats.getSumPhe(trait, eco);
}

double MetaPop::getSumTrait(const size_t &trait, const size_t &pop) const
{
    return pops[pop].getSumTrait(trait);
}


// Resetters used in tests

void MetaPop::consume() // metapopulation-level fitness/ecotype resetter
{
    for (size_t p = 0u; p < 2u; ++p)
        pops[p].consume();
}

void MetaPop::resetEcoTraits(const size_t &p, const double &x)
{
    pops[p].resetEcoTraits(x, ecosel, maxfeed);
}

void MetaPop::resetMatePrefs(const size_t &p, const double &y)
{
    pops[p].resetMatePrefs(y);
}

void MetaPop::resetEcotypes(const size_t &p, const size_t &e)
{
    pops[p].resetEcotypes(e);
}

void MetaPop::resetGenders(const size_t &p, const bool &sex)
{
    pops[p].resetGenders(sex);
    pops[p].sortSexes();
}
