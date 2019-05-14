//
// Created by p278834 on 7-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_POPULATION_H
#define EXPLICITGENOMESPECIATION_POPULATION_H

#include <vector>
#include <set>
#include "Individual.h"
#include "ParameterSet.h"

// Forward declaration
class Individual;
typedef Individual const * PInd;

class Population {

public:

    Population(const ParameterSet&, const GeneticArchitecture&);

    // High-level member functions
    void dispersal(const ParameterSet&);
    void sortByHabitat();
    void resourceDynamics(const size_t&, const double&);
    void reproduction(const size_t&, const ParameterSet&, const GeneticArchitecture&);
    void survival(const double&);

    // Low-level member functions
    void setResourceCapacities(const double&, const double&);
    void setReplenishRates(const double&);
    void setResourceConsumption(const size_t&);
    void setResourceEquilibrium(const size_t&);
    void assignFitnesses(const size_t&, const double&);
    void classifyGenders(const bool&);
    void setMaleFitnesses(const size_t&, const double&);
    void birth(const PInd&, const ParameterSet&, const GeneticArchitecture&);
    void emptyPopulation();

private:

    // The population
    std::list<PInd> individuals;
    std::vector<PInd> females;
    std::vector<PInd> males;
    std::vector<PInd> offspring;
    std::vector<std::pair<size_t, size_t> > genderCounts;
    std::vector<std::pair<size_t, size_t> > idHabitatBoundaries;

    // Mating features
    std::vector<double> maleSuccesses;

    // Ecological state
    size_t nAccessibleResources;
    std::vector<std::pair<double, double> > resourceCapacities;
    std::vector<std::pair<double, double> > replenishRates;
    std::vector<std::pair<double, double> > resourceConsumption;
    std::vector<std::pair<double, double> > resourceEql;
    std::pair<double, double> breakEvenPoint;

    // Genome-wide genetic variables
    std::vector<std::vector<double> > avgG;
    std::vector<std::vector<double> > varP;
    std::vector<std::vector<double> > varG;
    std::vector<std::vector<double> > varA;
    std::vector<std::vector<double> > varI;
    std::vector<double> varD;
    std::vector<double> Fst;
    std::vector<double> Pst;
    std::vector<double> Gst;
    std::vector<double> Qst;
    std::vector<double> Cst;

    // Locus-specific genetic variables
    struct LocusVariables {

        double avgEffectOfSubstitution;
        double locusvarD;
        double locusFit;
        double locusFis;
        double locusFst;
        double locusPst;
        double locusGst;
        double locusQst;
        double locusCst;
        std::array<double, 3u> alleleFrequency;
        std::array<double, 3u> meanEffect;
        std::array<double, 3u> locusvarP;
        std::array<double, 3u> locusvarG;
        std::array<double, 3u> locusvarA;
        std::array<double, 3u> locusvarI;

    };

    std::vector<LocusVariables> locusVariables;

};


#endif //EXPLICITGENOMESPECIATION_POPULATION_H
