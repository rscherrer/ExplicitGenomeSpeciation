#ifndef EXPLICITGENOMESPECIATION_POPFIXTURE_H
#define EXPLICITGENOMESPECIATION_POPFIXTURE_H

#include "library/utils.h"

struct PopFixture
{

    PopFixture() :
        pars(ParameterSet()),
        arch(GeneticArchitecture(pars)),
        n0(pars.getInitialPopSize()),
        s(pars.getEcoSelCoeff()),
        max(pars.getMaxFeedingRate()),
        k(rep(pars.getMaxResourceCapacity(), 2u)),
        r(rep(pars.getMaxResourceGrowth(), 2u))
    {}

    ~PopFixture() {}

    ParameterSet pars;
    GeneticArchitecture arch;
    size_t n0;  // initial population size
    double s;   // ecological selection coefficient
    double max; // maximum feeding rate
    vecDbl k;   // maximum resource capacities
    vecDbl r;   // maximum resource replenishment rate

};


#endif

