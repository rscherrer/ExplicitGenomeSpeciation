//
// Created by p278834 on 7-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_POPULATION_H
#define EXPLICITGENOMESPECIATION_POPULATION_H

#include <vector>
#include <set>
#include "Individual.h"
#include "ParameterSet.h"

class Individual;
typedef Individual const * PInd;

class Population {

public:

    std::list<PInd> individuals;

    std::vector<std::vector<double> > // 3 by 3
            avgG, varP, varG, varA, varI;
    std::vector<double> varD, F_st, P_st, G_st, Q_st, C_st; // size nCharacter

    void dispersal(const ParameterSet&);

    void competitionAndReproduction(const size_t,
                                    const ParameterSet&,
                                    std::vector<std::pair<double, double> >&,
                                    Individual::TradeOffPt&,
                                    std::vector<std::pair<double, double> >&,
                                    std::vector<std::pair<size_t, size_t> >&,
                                    const Genome&,
                                    const size_t)




private:

};


#endif //EXPLICITGENOMESPECIATION_POPULATION_H
