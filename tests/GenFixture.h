#ifndef EXPLICITGENOMESPECIATION_GENFIXTURE_H
#define EXPLICITGENOMESPECIATION_GENFIXTURE_H

#include "library/GeneticArchitecture.h"

typedef std::vector<Network> MultiNet;

struct GenFixture {

    GenFixture() :
        pars(ParameterSet()),
        arch(GeneticArchitecture(pars)),
        genome(arch.getGenome()),
        networks(arch.getNetworks())
    {}
    ~GenFixture() {}

    ParameterSet pars;
    GeneticArchitecture arch;
    Genome genome;
    MultiNet networks;

};

#endif
