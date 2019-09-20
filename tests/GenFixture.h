#ifndef EXPLICITGENOMESPECIATION_GENFIXTURE_H
#define EXPLICITGENOMESPECIATION_GENFIXTURE_H

#include "library/GenArch.h"

typedef std::vector<Network> MultiNet;

struct GenFixture {

    GenFixture() :
        pars(Param()),
        arch(GenArch(pars))
    {}
    ~GenFixture() {}

    Param pars;
    GenArch arch;

};

#endif
