#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Population.h"
#include "Buffer.h"
#include "StreamBag.h"
#include "Individual.h"
#include "utils.h"
#include "types.h"

typedef std::vector<Population> vecPop;
typedef std::vector<std::ofstream *> vecStreams;

class MetaPop
{

public:

    MetaPop(const vecUns& popSizes, const ParameterSet &pars,
     const GeneticArchitecture &arch) :
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
        resources(matzeros(2u, 2u)),
        replenish(matzeros(2u, 2u)),
        tmax(pars.getTEndSim()),
        tsave(pars.getTSave()),
        tburnin(pars.getTBurnIn()),
        record(pars.getRecord()),
        buffer(Buffer()),
        ecotypes({ { }, { } }),
        meanPhenotypes(matzeros(3u, 3u)),
        pheVariances(matzeros(3u, 3u)),
        genVariances(matzeros(3u, 3u)),
        addVariances(matzeros(3u, 3u)),
        nadVariances(matzeros(3u, 3u)),
        domVariances(zeros(3u)),
        intVariances(zeros(3u)),
        Pst(zeros(3u)),
        Gst(zeros(3u)),
        Qst(zeros(3u)),
        Cst(zeros(3u)),
        Fst(zeros(3u)),
        varPScan(zeros(pars.getNLoci())),
        varGScan(zeros(pars.getNLoci())),
        varAScan(zeros(pars.getNLoci())),
        varNScan(zeros(pars.getNLoci())),
        PstScan(zeros(pars.getNLoci())),
        GstScan(zeros(pars.getNLoci())),
        QstScan(zeros(pars.getNLoci())),
        CstScan(zeros(pars.getNLoci())),
        FstScan(zeros(pars.getNLoci()))
    {
        // Create the demes
        for (size_t p = 0u; p < 2u; ++p) {
            for (size_t res = 0u; res < 2u; ++res) {
                resources[p][res] = maxresources * fabs(p - res) * symmetry;
                replenish[p][res] = maxreplenish;
            }
            const size_t n = popsizes[p];
            const double max = maxfeed;
            const vecDbl k = resources[p];
            const vecDbl r = replenish[p];
            pops.push_back(Population(n, ecosel, max, k, r, arch));
        }
    }
    ~MetaPop() {}

    vecPop getPops() const { return pops; }
    double getEcoIsolation();
    double getSpatialIsolation();
    double getMatingIsolation();

    int evolve(const GeneticArchitecture&);
    void analyze(const GeneticArchitecture&);
    void loadBuffer(const size_t &t);
    void save(StreamBag&);

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
