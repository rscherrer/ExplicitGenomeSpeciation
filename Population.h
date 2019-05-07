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

    // A Character object is locus-specific and population-wide
    struct Character
    {
        size_t character, linkageGroup;
        double location, effectSize, dominanceCoeff, avgEffectOfSubstitution,
                varD, F_it, F_is, F_st, P_st, G_st, Q_st, C_st;
        std::array<double, 3u> alleleFrequency, meanEffect, varP, varG, varA, varI;
        std::list<std::pair<size_t, double> > edges;
    };

    // Static variables -- these should belong to population

    std::vector<std::vector<double> > // 3 by 3
            avgG, varP, varG, varA, varI;
    std::vector<double> varD, F_st, P_st, G_st, Q_st, C_st; // size nCharacter
    std::vector<double> chromosomeSize; // size nchromosomes - 1
    std::vector<std::set<size_t> > vertices; // size nCharacter
    std::vector<Character> characterLocus; // size nLoci


private:

};


#endif //EXPLICITGENOMESPECIATION_POPULATION_H
