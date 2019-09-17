#ifndef EXPLICITGENOMESPECIATION_INDIVIDUAL_H
#define EXPLICITGENOMESPECIATION_INDIVIDUAL_H

#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "utils.h"
#include <vector>
#include <random>
#include <stddef.h>


typedef std::vector<bool> Haplotype;
typedef std::vector<Haplotype > Diplotype;
typedef std::vector<double> vecDbl;

/// Function to calculate feeding rates
vecDbl calcFeedingRates(const double&, const double&,
 const double& = 0.0004);

class Individual {

public:

    typedef Individual const * PInd;

    Individual(const Genome&, const MultiNet&, const double& = 0.5,
     const vecDbl& = { 1.0, 1.0, 1.0 }, const vecDbl& = zeros(3u),
      const vecDbl& = zeros(3u), const vecDbl& = zeros(3u));
    Individual(const Genome&, const MultiNet&, const Haplotype&,
     const Haplotype&, const vecDbl& = { 1.0, 1.0, 1.0 },
      const vecDbl& = zeros(3u), const vecDbl& = zeros(3u),
       const vecDbl& = zeros(3u));
    ~Individual() {}

    // Getters
    bool getGender() const { return isFemale; }
    double getEcoTrait() const { return ecoTrait; }
    double getMatePref() const { return matePref; }
    double getNeutral() const { return neutral; }
    double getFitness() const { return fitness; }
    vecDbl getTraits() const { return traitvalues; }
    vecDbl getFeedingRates() const { return feedingRates; }
    vecDbl getGenValues() const { return genvalues; }
    Diplotype getSequence() const { return sequence; }
    vecDbl getExpression() const { return genexp; }
    size_t getEcotype() const { return ecotype; }
    size_t getZygosity(const size_t&);
    double getLocusValue(const size_t&);


    // Actions
    void feed(const vecDbl&);
    bool acceptMate(const double&, const double&) const;
    Haplotype recombine(const vecDbl&, const vecDbl&,
     const double& = 3.0);
    void mutate(Haplotype&, const double& = 1.0e-5);

    // Setters
    void setEcoTrait(const double &value, const double &sel) {
        ecoTrait = value;
        feedingRates = calcFeedingRates(sel, value);
    }
    void setMatePref(const double &value) { matePref = value; }
    void setEcotype(const double&);



private:

    friend class Population;

    // Makers
    Diplotype makeSequence(const size_t&, const double& = 0.5);
    Diplotype fecundate(const Haplotype&, const Haplotype&);
    bool determineSex(const bool&);

    // Setters
    void develop(const Genome&, const MultiNet&, const vecDbl&, const vecDbl&,
     const vecDbl&, const vecDbl&);

    // Fields
    Diplotype sequence;
    vecDbl genexp;
    vecDbl locivalues;
    bool isFemale;
    vecDbl genvalues;
    vecDbl traitvalues;
    double ecoTrait;
    double matePref;
    double neutral;    
    double fitness;
    vecDbl feedingRates;
    size_t ecotype;


};

#endif


