#ifndef EXPLICITGENOMESPECIATION_GENFIXTURE_H
#define EXPLICITGENOMESPECIATION_GENFIXTURE_H

#include "library/GeneticArchitecture.h"

struct GenFixture {

    GenFixture() :
        pars(ParameterSet()),
        arch(GeneticArchitecture(pars)),
        genome(arch.getGenome())
    {}
    ~GenFixture() {}

    ParameterSet pars;
    GeneticArchitecture arch;
    Genome genome;

};

#endif
