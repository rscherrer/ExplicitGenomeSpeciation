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

    MetaPop(const vecPop &populations, const ParameterSet &pars) :
        pops(populations),
        dispersal(pars.getDispersalRate()),
        survival(pars.getSurvivalProb()),
        birth(pars.getBirthRate()),
        matingcost(pars.getMateEvaluationCost()),
        sexsel(pars.getMatePreferenceStrength()),
        ecosel(pars.getEcoSelCoeff()),
        tmax(pars.getTEndSim()),
        tsave(pars.getTSave()),
        tburnin(pars.getTBurnIn()),
        record(pars.getRecord()),
        buffer(Buffer()),
        ecotypes({ { }, { } }),
        meanPhenotypes({ zeros(3u), zeros(3u), zeros(3u) }),
        pheVariances({ zeros(3u), zeros(3u), zeros(3u) }),
        genVariances({ zeros(3u), zeros(3u), zeros(3u) }),
        addVariances({ zeros(3u), zeros(3u), zeros(3u) }),
        nadVariances({ zeros(3u), zeros(3u), zeros(3u) }),
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
    {}
    ~MetaPop() {}

    // Getters
    vecPop getPops() const { return pops; }
    double getEcoIsolation();
    double getSpatialIsolation();
    double getMatingIsolation();

    // Setters
    int evolve(const Genome&, const MultiNet&, const vecDbl& = zeros(3u));
    void analyze(const size_t&, const vecUns&, const vecDbl&);
    void loadBuffer(const size_t &t);
    void save(StreamBag&);

private:

    vecPop pops;
    double dispersal;
    double survival;
    double birth;
    double matingcost;
    double sexsel;
    double ecosel;
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
