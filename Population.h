//
// Created by p278834 on 7-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_POPULATION_H
#define EXPLICITGENOMESPECIATION_POPULATION_H

#include <vector>
#include <set>
#include "individual.h"

class Population {

public:



    // Static variables -- these should belong to population

    std::vector<std::vector<double> > // 3 by 3
            avgG, varP, varG, varA, varI;
    std::vector<double> varD, F_st, P_st, G_st, Q_st, C_st; // size nCharacter



    // How about a class Population for statistics measured across individuals, and a class Genome for the statistics measured on a per-locus basis?
    // That would make more sense


private:

};


#endif //EXPLICITGENOMESPECIATION_POPULATION_H
