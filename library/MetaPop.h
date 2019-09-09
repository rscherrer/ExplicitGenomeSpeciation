#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Population.h"

typedef std::vector<Population> vecPop;

class Buffer
{

    friend class MetaPop;

private:
    char time[25];
    char popsize0[25];
    char popsize1[25];

public:
    Buffer() {}
    ~Buffer() {}
    void load(const size_t&);
    void write(std::ofstream&);
};

class MetaPop
{

public:

    MetaPop(const vecPop &populations, const ParameterSet &pars) :
        pops(populations),
        dispersal(pars.getDispersalRate()),
        survival(pars.getSurvivalProb()),
        birth(pars.getBirthRate()),
        mating(pars.getMatePreferenceStrength()),
        tmax(pars.getTEndSim()),
        tsave(pars.getTSave()),
        buffer(Buffer())
    {}
    ~MetaPop() {}

    vecPop getPops() const { return pops; }

    size_t evolve(const Genome&, const MultiNet&);
    void loadBuffer(const size_t &t);

private:

    vecPop pops;
    double dispersal;
    double survival;
    double birth;
    double mating;
    size_t tmax;
    size_t tsave;
    Buffer buffer;

};

#endif
