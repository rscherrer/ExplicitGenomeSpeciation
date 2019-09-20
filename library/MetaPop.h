#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Deme.h"
#include "Buffer.h"
#include "Output.h"
#include "Individual.h"
#include "Utilities.h"
#include "types.h"

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
        tmax(pars.getTEndSim()),
        tsave(pars.getTSave()),
        tburnin(pars.getTBurnIn()),
        record(pars.getRecord()),
        buffer(Buffer()),
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
        for (int p = 0; p < 2; ++p) {
            for (int res = 0; res < 2; ++res) {
                resources[p][res] = maxresources * fabs(p - res) * symmetry;
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
    double getEcoIsolation();
    double getSpatialIsolation();
    double getMatingIsolation();

    int evolve(const GenArch&);
    void analyze(const GenArch&);
    void loadBuffer(const size_t &t);
    void save(Output&);

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
    int tmax;
    int tsave;
    int tburnin;
    bool record;
    Buffer buffer;

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
