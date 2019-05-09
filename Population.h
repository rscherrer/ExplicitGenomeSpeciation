//
// Created by p278834 on 7-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_POPULATION_H
#define EXPLICITGENOMESPECIATION_POPULATION_H

#include <vector>
#include <set>
#include "Individual.h"
#include "ParameterSet.h"
#include "Genome.h"

// Forward declaration
class Individual;
typedef Individual const * PInd;
typedef std::pair<double, double> TradeOffPt;

class Population {

public:

    // The population
    std::list<PInd> individuals;

    // Ecological state
    std::vector<std::pair<double, double> > resourceConsumption;
    std::vector<std::pair<double, double> > resourceEql;
    std::vector<std::pair<size_t, size_t> > genderCounts;
    TradeOffPt breakEvenPoint;
    size_t nAccessibleResource;

    // Genome-wide genetics
    std::vector<std::vector<double> > avgG;
    std::vector<std::vector<double> > varP;
    std::vector<std::vector<double> > varG;
    std::vector<std::vector<double> > varA;
    std::vector<std::vector<double> > varI;
    std::vector<double> varD;
    std::vector<double> F_st;
    std::vector<double> P_st;
    std::vector<double> G_st;
    std::vector<double> Q_st;
    std::vector<double> C_st;

    // Locus-specific genetics
    struct LocusVariables {

        double avgEffectOfSubstitution;
        double varD;
        double F_it;
        double F_is;
        double F_st;
        double P_st;
        double G_st;
        double Q_st;
        double C_st;
        std::array<double, 3u> alleleFrequency;
        std::array<double, 3u> meanEffect;
        std::array<double, 3u> varP;
        std::array<double, 3u> varG;
        std::array<double, 3u> varA;
        std::array<double, 3u> varI;

    };

    std::vector<LocusVariables> variablesPerLocus;

    // Member functions
    void dispersal(const ParameterSet&);
    void competitionAndReproduction(const size_t, const ParameterSet&, const Genome&);

private:

};


#endif //EXPLICITGENOMESPECIATION_POPULATION_H
