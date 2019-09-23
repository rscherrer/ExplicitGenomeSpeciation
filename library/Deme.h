#ifndef EXPLICITGENOMESPECIATION_DEME_H
#define EXPLICITGENOMESPECIATION_DEME_H

#include "Param.h"
#include "GenArch.h"
#include "Random.h"
#include "Utilities.h"
#include "Types.h"
#include <stddef.h>
#include <iostream>

class Individual;

typedef Individual * PInd;
typedef std::vector<PInd> Crowd;

class Deme {

    friend class MetaPop;

public:

    Deme(const size_t &popsize, const double &ecosel,
     const double &maxfeeding, const vecDbl &foodmax, const vecDbl &foodgrowth,
      const GenArch &arch, const bool &isBurnin = false) :
        individuals(populate(popsize, ecosel, maxfeeding, arch)),
        females({ }),
        males({ }),
        offspring({ }),
        survivors({ }),
        capacity(foodmax),
        replenish(foodgrowth),
        resources(capacity),
        burnin(isBurnin)
    {
        // Sexes should be sorted out upon creation of the population
        sortSexes();
    }

    ~Deme() {};

    // Getters
    size_t getPopSize() const { return individuals.size(); }
    size_t getNOffspring() const { return offspring.size(); }
    size_t getNFemales() const { return females.size(); }
    size_t getNMales() const { return males.size(); }
    double getResource(const size_t &r) const { return resources[r]; }
    PInd getInd(const size_t &i) const { return individuals[i]; }

    // Life cycle
    void sortSexes();
    Crowd emigrate(const double& = 0.01);
    void immigrate(const Crowd&);
    void consume();
    void reproduce(const double&, const double&, const double&, const double&,
     const double&, const GenArch&);
    bool survive(const double&);
    void exitBurnIn();

    void resetEcoTraits(const double&, const double&, const double&);
    void resetMatePrefs(const double&);
    void resetGenders(const bool&);
    void resetEcotypes(const size_t&);

private:

    // Makers
    Crowd populate(const size_t&, const double&, const double&,
     const GenArch&);

    // The population
    Crowd individuals;
    Crowd females;
    Crowd males;
    Crowd offspring;
    Crowd survivors;

    vecDbl capacity;
    vecDbl replenish;
    vecDbl resources;

    bool burnin;

};


#endif
