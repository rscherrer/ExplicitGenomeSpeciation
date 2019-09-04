#ifndef EXPLICITGENOMESPECIATION_METAPOP_H
#define EXPLICITGENOMESPECIATION_METAPOP_H

#include "Population.h"

typedef std::vector<Population> vecPop;

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
        tsave(pars.getTSave())
    {}
    ~MetaPop() {}

private:

    vecPop pops;
    double dispersal;
    double survival;
    double birth;
    double mating;
    size_t tmax;
    size_t tsave;

};

#endif
