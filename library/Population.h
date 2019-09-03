#ifndef EXPLICITGENOMESPECIATION_POPULATION_H
#define EXPLICITGENOMESPECIATION_POPULATION_H

#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "Random.h"
#include <stddef.h>


class Individual;
class LocusVariables;
typedef Individual * PInd;
typedef std::vector<PInd> Crowd;
typedef std::vector<double> dVector;

class Population {

public:

    Population(const size_t&, const Genome&, const MultiNet&,
     const dVector& = {100.0, 100.0}, const dVector& = {1.0, 1.0});
    ~Population() {};

    // Getters
    size_t getPopSize() const { return individuals.size(); }
    size_t getNOffspring() const { return offspring.size(); }
    dVector getResources() const { return resources; }

    // Life cycle
    Crowd emigrate(const double& = 0.01);
    void immigrate(const Crowd&);
    void consume();
    void reproduce(const double&, const double&, const Genome&,
     const MultiNet&, const double& = 0.01);
    void reproduceAsexual(const double&, const Genome&,
     const MultiNet&);
    bool survive(const double&);



private:

    // Makers
    Crowd populate(const size_t&, const Genome&,
     const MultiNet&);

    // The population
    Crowd individuals;
    Crowd females;
    Crowd males;
    Crowd offspring;
    Crowd survivors;

    dVector capacity;
    dVector replenish;
    dVector resources;



};


#endif
