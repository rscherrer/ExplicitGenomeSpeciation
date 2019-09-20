#ifndef EXPLICITGENOMESPECIATION_POPFIXTURE_H
#define EXPLICITGENOMESPECIATION_POPFIXTURE_H

#include "library/Param.h"
#include "library/GenArch.h"
#include "library/Utilities.h"

struct PopFixture
{

    PopFixture() :
        pars(Param()),
        arch(GenArch(pars)),
        n0(pars.getInitialPopSize()),
        s(pars.getEcoSelCoeff()),
        max(pars.getMaxFeedingRate()),
        k(utl::rep(pars.getMaxResourceCapacity(), 2u)),
        r(utl::rep(pars.getMaxResourceGrowth(), 2u))
    {}

    ~PopFixture() {}

    Param pars;
    GenArch arch;
    size_t n0;  // initial population size
    double s;   // ecological selection coefficient
    double max; // maximum feeding rate
    vecDbl k;   // maximum resource capacities
    vecDbl r;   // maximum resource replenishment rate

};


#endif

