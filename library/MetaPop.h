#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Population.h"
#include "Buffer.h"
#include "StreamBag.h"
#include "Individual.h"
#include "utils.h"

typedef std::vector<Population> vecPop;
typedef std::vector<std::ofstream *> vecStreams;
typedef std::vector<std::string> vecStrings;

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
        tmax(pars.getTEndSim()),
        tsave(pars.getTSave()),
        buffer(Buffer()),
        record(pars.getRecord()),
        ecomean(0.0)
    {}
    ~MetaPop() {}

    vecPop getPops() const { return pops; }

    size_t evolve(const Genome&, const MultiNet&);
    void loadBuffer(const size_t &t);

    double getEcoIsolation();
    double getSpatialIsolation();
    double getMatingIsolation();

private:

    vecPop pops;
    double dispersal;
    double survival;
    double birth;
    double matingcost;
    double sexsel;
    size_t tmax;
    size_t tsave;
    Buffer buffer;
    bool record;
    double ecomean;

};

#endif
