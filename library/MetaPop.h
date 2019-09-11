#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Population.h"

typedef std::vector<Population> vecPop;
typedef std::vector<std::ofstream *> vecStreams;
typedef std::vector<std::string> vecStrings;

class StreamBag
{

    friend class MetaPop;
    friend class Buffer;

private:

    vecStreams files;
    vecStrings names;

public:

    StreamBag() :
        files({ }),
        names({
              "time",
              "popsize0",
              "popsize1",
              "nfemales0",
              "nfemales1",
              "resource00",
              "resource01",
              "resource10",
              "resource11",
              "meanx0",
              "meanx1",
              "meany0",
              "meany1",
              "meanz0",
              "meanz1",
              "ecological_iso",
              "spatial_iso",
              "mating_iso"
        })
    {}
    ~StreamBag() {}

    void openAll();
    void closeAll();


};

class Buffer
{

    friend class MetaPop;

private:

    vecDbl fields;

public:

    Buffer() : fields({ }) {}
    ~Buffer() {}

    void flush();
    void add(const double&);
    void write(std::ofstream *&, const double&);
};

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
