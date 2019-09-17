#ifndef EXPLICITGENOMESPECIATION_POPULATION_H
#define EXPLICITGENOMESPECIATION_POPULATION_H

#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "Random.h"
#include "utils.h"
#include <stddef.h>


class Individual;
class LocusVariables;
typedef Individual * PInd;
typedef std::vector<PInd> Crowd;
typedef std::vector<double> vecDbl;

class Population {

    friend class MetaPop;

public:

    Population(const size_t&, const Genome&, const MultiNet&,
     const vecDbl& = {100.0, 100.0}, const vecDbl& = {1.0, 1.0});
    ~Population() {};

    // Getters
    size_t getPopSize() const { return individuals.size(); }
    size_t getNOffspring() const { return offspring.size(); }
    size_t getNFemales() const { return females.size(); }
    size_t getNMales() const { return males.size(); }
    vecDbl getResources() const { return resources; }
    double getMeanEcoTrait() const { return meanEcoTrait; }
    double getMeanMatePref() const { return meanMatePref; }
    double getMeanNtrTrait() const { return meanNtrTrait; }

    // Life cycle
    void sortSexes();
    Crowd emigrate(const double& = 0.01);
    void immigrate(const Crowd&);
    void consume();
    void burninConsume();
    void reproduce(const double&, const double&, const Genome&,
     const MultiNet&, const double& = 0.01);
    void burninReproduce(const double&, const double&, const Genome&,
     const MultiNet&, const double& = 0.01, const double& = 0.0);
    void reproduceAsexual(const double&, const Genome&,
     const MultiNet&);
    bool survive(const double&);

    // Setters
    void calcMeanEcoTrait();
    void calcMeanMatePref();
    void calcMeanNtrTrait();


private:

    // Makers
    Crowd populate(const size_t&, const Genome&, const MultiNet&,
     const double& = 0.5, const vecDbl& = { 1.0, 1.0, 1.0 },
      const vecDbl& = zeros(3u), const vecDbl& = zeros(3u),
       const vecDbl& = zeros(3u));

    // The population
    Crowd individuals;
    Crowd females;
    Crowd males;
    Crowd offspring;
    Crowd survivors;

    vecDbl capacity;
    vecDbl replenish;
    vecDbl resources;

    double meanEcoTrait = 0.0;
    double meanMatePref = 0.0;
    double meanNtrTrait = 0.0;


};


#endif
