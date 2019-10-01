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
        n0(100u),
        s(1.0),
        max(0.004),
        k(utl::rep(100.0, 2u)),
        r(utl::rep(1.0, 2u))
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

