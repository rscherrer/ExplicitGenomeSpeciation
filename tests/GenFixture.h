#ifndef EXPLICITGENOMESPECIATION_GENFIXTURE_H
#define EXPLICITGENOMESPECIATION_GENFIXTURE_H

#include "library/GeneticArchitecture.h"

typedef std::vector<Network> MultiNet;

struct GenFixture {

    GenFixture() :
        pars(ParameterSet()),
        arch(GenArch(pars))
    {}
    ~GenFixture() {}

    ParameterSet pars;
    GenArch arch;

};

#endif
