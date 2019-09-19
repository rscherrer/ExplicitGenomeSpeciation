#ifndef EXPLICITGENOMESPECIATION_GENFIXTURE_H
#define EXPLICITGENOMESPECIATION_GENFIXTURE_H

#include "library/GeneticArchitecture.h"

typedef std::vector<Network> MultiNet;

struct GenFixture {

    GenFixture() :
        pars(ParameterSet()),
        arch(GeneticArchitecture(pars))
    {}
    ~GenFixture() {}

    ParameterSet pars;
    GeneticArchitecture arch;

};

#endif
