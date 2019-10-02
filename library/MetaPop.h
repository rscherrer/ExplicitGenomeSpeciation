#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Deme.h"
#include "Output.h"
//#include "Individual.h"
#include "Utilities.h"
#include "Types.h"
#include "Stats.h"
#include <cassert>

typedef std::vector<Deme> vecPop;
typedef std::vector<std::ofstream *> vecStreams;

class MetaPop
{

public:

    MetaPop(const vecUns& popSizes, const Param &pars,
     const GenArch &arch, const bool &isburnin) :
        popsizes(popSizes),
        dispersal(pars.getDispersalRate()),
        survival(pars.getSurvivalProb()),
        birth(pars.getBirthRate()),
        matingcost(pars.getMateEvaluationCost()),
        sexsel(pars.getMatePreferenceStrength()),
        ecosel(pars.getEcoSelCoeff()),
        symmetry(pars.getHabitatSymmetry()),
        maxfeed(pars.getMaxFeedingRate()),
        maxresources(pars.getMaxResourceCapacity()),
        maxreplenish(pars.getMaxResourceGrowth()),
        resources(utl::matzeros(2u, 2u)),
        replenish(utl::matzeros(2u, 2u)),
        t(0),
        tmax(pars.getTEndSim()),
        tsave(pars.getTSave()),
        tburnin(pars.getTBurnIn()),
        record(pars.getRecord()),
        pops(makeDemes(arch, isburnin)),
        stats(Stats(arch))
    {}

    ~MetaPop() {}

    size_t getPopSize(const size_t &p) const { return pops[p].getPopSize(); }
    size_t getNFemales(const size_t &p) const { return pops[p].getNFemales(); }
    size_t getNOffspring(const size_t&) const;
    size_t getSumEcotypes(const size_t&) const;
    double getResource(const size_t&, const size_t&) const;
    double getEcoIsolation() const { return stats.getEcoIsolation(); }
    double getSpatialIsolation() const { return stats.getSpatialIsolation(); }
    double getMatingIsolation() const { return stats.getMatingIsolation(); }
    double getPst(const size_t &trait) const { return stats.getPst(trait); }
    double getVarP(const size_t&, const size_t&) const;

    int evolve(const GenArch&);
    void analyze(const GenArch&);
    void consume();

    void resetEcoTraits(const size_t&, const double&);
    void resetMatePrefs(const size_t&, const double&);
    void resetEcotypes(const size_t&, const size_t&);
    void resetGenders(const size_t&, const bool&);

private:

    vecPop makeDemes(const GenArch&, const bool&);

    vecUns popsizes;
    double dispersal;
    double survival;
    double birth;
    double matingcost;
    double sexsel;
    double ecosel;
    double symmetry;
    double maxfeed;
    double maxresources;
    double maxreplenish;
    Matrix resources;
    Matrix replenish;
    int t;
    int tmax;
    int tsave;
    int tburnin;
    bool record;

    vecPop pops;
    Stats stats;

    static constexpr size_t x = 0u;
    static constexpr size_t y = 1u;
    static constexpr size_t z = 2u;

};

#endif
