#ifndef EXPLICITGENOMESPECIATION_TESTFIXTUREDEFAULTPARAMETERANDARCHITECTUREFILES_H
#define EXPLICITGENOMESPECIATION_TESTFIXTUREDEFAULTPARAMETERANDARCHITECTUREFILES_H

#include "GeneticArchitecture.h"
#include "ParameterSet.h"



/// Test fixture for default parameters and architecture file
struct createPopulationStartingKit
{

    createPopulationStartingKit()
    {
        geneticArchitecture.generateGeneticArchitecture(parameters);
    }

    ~createPopulationStartingKit() {}

    ParameterSet parameters;
    GeneticArchitecture geneticArchitecture;

};


#endif //EXPLICITGENOMESPECIATION_TESTFIXTUREDEFAULTPARAMETERANDARCHITECTUREFILES_H
