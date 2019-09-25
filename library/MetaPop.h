#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Deme.h"
#include "Output.h"
//#include "Individual.h"
#include "Utilities.h"
#include "Types.h"
#include <cassert>

typedef std::vector<Deme *> vecPop;
typedef std::vector<std::ofstream *> vecStreams;

class MetaPop
{

public:

    MetaPop(const vecUns& popSizes, const Param &pars,
     const GenArch &arch) :
        pops({ }),
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
        ecotypes({ { }, { } }),
        meanPhenotypes(utl::matzeros(3u, 3u)),
        pheVariances(utl::matzeros(3u, 3u)),
        genVariances(utl::matzeros(3u, 3u)),
        addVariances(utl::matzeros(3u, 3u)),
        nadVariances(utl::matzeros(3u, 3u)),
        domVariances(utl::zeros(3u)),
        intVariances(utl::zeros(3u)),
        Pst(utl::zeros(3u)),
        Gst(utl::zeros(3u)),
        Qst(utl::zeros(3u)),
        Cst(utl::zeros(3u)),
        Fst(utl::zeros(3u)),
        varPScan(utl::zeros(pars.getNLoci())),
        varGScan(utl::zeros(pars.getNLoci())),
        varAScan(utl::zeros(pars.getNLoci())),
        varNScan(utl::zeros(pars.getNLoci())),
        PstScan(utl::zeros(pars.getNLoci())),
        GstScan(utl::zeros(pars.getNLoci())),
        QstScan(utl::zeros(pars.getNLoci())),
        CstScan(utl::zeros(pars.getNLoci())),
        FstScan(utl::zeros(pars.getNLoci()))
    {
        // Create the demes
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
            pops.push_back(new Deme(n, ecosel, max, k, r, arch, true));
        }
    }

    ~MetaPop() {}

    size_t getPopSize(const size_t &p) const { return pops[p]->getPopSize(); }
    size_t getNFemales(const size_t &p) const { return pops[p]->getNFemales(); }
    size_t getNOffspring(const size_t&) const;
    size_t getEcotypeSize(const size_t &e) const { return ecotypes[e].size(); }
    double getResource(const size_t&, const size_t&) const;
    double getEcoIsolation();
    double getSpatialIsolation();
    double getMatingIsolation();

    int evolve(const GenArch&);
    void analyze(const GenArch&);
    void save(Output&);
    void write(const double&, std::ofstream *&);
    void write(const vecDbl&, std::ofstream *&);

    void resetEcoTraits(const size_t&, const double&);
    void resetMatePrefs(const size_t&, const double&);
    void resetEcotypes(const size_t&, const size_t&);
    void resetGenders(const size_t&, const bool&);

private:

    vecPop pops;
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

    // Variables for analysis
    std::vector<Crowd> ecotypes;
    Matrix meanPhenotypes;
    Matrix pheVariances;
    Matrix genVariances;
    Matrix addVariances;
    Matrix nadVariances;
    vecDbl domVariances;
    vecDbl intVariances;
    vecDbl Pst;
    vecDbl Gst;
    vecDbl Qst;
    vecDbl Cst;
    vecDbl Fst;
    vecDbl varPScan;
    vecDbl varGScan;
    vecDbl varAScan;
    vecDbl varNScan;
    vecDbl PstScan;
    vecDbl GstScan;
    vecDbl QstScan;
    vecDbl CstScan;
    vecDbl FstScan;

};

#endif
