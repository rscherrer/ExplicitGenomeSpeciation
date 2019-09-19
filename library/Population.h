#ifndef EXPLICITGENOMESPECIATION_POPULATION_H
#define EXPLICITGENOMESPECIATION_POPULATION_H

#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "Random.h"
#include "utils.h"
#include "types.h"
#include <stddef.h>

class Individual;

typedef Individual * PInd;
typedef std::vector<PInd> Crowd;

class Population {

    friend class MetaPop;

public:

    Population(const size_t &popsize, const GeneticArchitecture &arch,
     const vecDbl &foodmax = rep(100.0, 2u),
      const vecDbl &foodgrowth = ones(2u)) :
        individuals(populate(popsize, arch)),
        females({ }),
        males({ }),
        offspring({ }),
        survivors({ }),
        capacity(foodmax),
        replenish(foodgrowth),
        resources(capacity)
    {
        // Sexes should be sorted out upon creation of the population
        sortSexes();
    }

    ~Population() {};

    // Getters
    size_t getPopSize() const { return individuals.size(); }
    size_t getNOffspring() const { return offspring.size(); }
    size_t getNFemales() const { return females.size(); }
    size_t getNMales() const { return males.size(); }
    vecDbl getResources() const { return resources; }
    PInd getInd(const size_t &i) const { return individuals[i]; }

    // Life cycle
    void sortSexes();
    Crowd emigrate(const double& = 0.01);
    void immigrate(const Crowd&);
    void consume();
    void burninConsume();
    void reproduce(const double&, const double&, const double&,
     const GeneticArchitecture&);
    void burninReproduce(const double&, const double&, const double&,
     const double&, const GeneticArchitecture&);
    bool survive(const double&);

    void resetEcoTraits(const double&, const double&);
    void resetMatePrefs(const double&);
    void resetGenders(const bool&);
    void resetEcotypes(const size_t&);

private:

    // Makers
    Crowd populate(const size_t&, const GeneticArchitecture&);

    // The population
    Crowd individuals;
    Crowd females;
    Crowd males;
    Crowd offspring;
    Crowd survivors;

    vecDbl capacity;
    vecDbl replenish;
    vecDbl resources;

};


#endif
