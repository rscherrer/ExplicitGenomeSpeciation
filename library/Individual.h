#ifndef EXPLICITGENOMESPECIATION_INDIVIDUAL_H
#define EXPLICITGENOMESPECIATION_INDIVIDUAL_H

//#include "Param.h"
#include "GenArch.h"
#include "Utilities.h"
#include "Types.h"
//#include "Deme.h"
#include "Random.h"
//#include <vector>
//#include <random>
//#include <cmath> // exponential?
#include <cassert>
//#include <algorithm>
//#include <stddef.h>
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<> Haplotype;
typedef std::vector<Haplotype> Genome;

class Individual {

    friend class Population;
    friend class MetaPop;

    typedef Individual const * PInd;

public:

    /// Spontaneous creation
    Individual(const GenArch &arch, const double &ecosel,
     const double &maxfeeding, const double &snpFreq) :
        genome(generateGenome(arch, snpFreq)),
        transcriptome(utl::zeros(arch.nLoci)),
        locivalues(utl::zeros(arch.nLoci)),
        gender(rnd::bernoulli(0.5)),
        genvalues(utl::zeros(3u)),
        traitvalues(utl::zeros(3u)),
        ecotrait(0.0),
        matepref(0.0),
        neutrait(0.0),
        fitness(1.0),
        feedingRates(calcFeedingRates(ecosel, ecotrait, maxfeeding)),
        ecotype(0u)
    {
        develop(arch);
        checkIndividual(arch.nLoci);
    }


    /// Newborn
    Individual(const GenArch &arch, const Haplotype &egg,
     const Haplotype &sperm, const double &ecosel, const double &maxfeeding) :
        genome(fecundate(egg, sperm)),
        transcriptome(utl::zeros(arch.nLoci)),
        locivalues(utl::zeros(arch.nLoci)),
        gender(rnd::bernoulli(0.5)),
        genvalues(utl::zeros(3u)),
        traitvalues(utl::zeros(3u)),
        ecotrait(0.0),
        matepref(0.0),
        neutrait(0.0),
        fitness(1.0),
        feedingRates(calcFeedingRates(ecosel, ecotrait, maxfeeding)),
        ecotype(0u)
    {
        develop(arch);
        checkIndividual(arch.nLoci);
    }

    ~Individual() {}

    bool getGender() const { return gender; }
    double getEcoTrait() const { return ecotrait; }
    double getMatePref() const { return matepref; }
    double getNeutral() const { return neutrait; }
    double getFitness() const { return fitness; }
    double getTraitValue(const size_t &t) const { return traitvalues[t]; }
    double getGenValue(const size_t &t) const { return genvalues[t]; }
    vecDbl getTraits() const { return traitvalues; }
    vecDbl getFeedingRates() const { return feedingRates; }
    vecDbl getGenValues() const { return genvalues; }
    Haplotype getSequence(const size_t &i) const { return genome[i]; }
    vecDbl getExpression() const { return transcriptome; }
    size_t getEcotype() const { return ecotype; }
    size_t getZygosity(const size_t&);
    double getLocusValue(const size_t&);
    size_t getAlleleSum(const size_t&);
    bool acceptMate(const double&, const double&) const;
    Haplotype recombine(const GenArch&);
    vecDbl calcFeedingRates(const double&, const double&, const double&);

    void setEcoTrait(const double&, const double&, const double&);
    void setMatePref(const double&);
    void resetEcotype(const size_t &e) { ecotype = e; }
    void setGender(const bool&);
    void feed(const vecDbl&);
    void mutate(Haplotype&, const double& = 1.0e-5);

private:

    Genome generateGenome(const GenArch&, double = -1.0);
    Genome fecundate(const Haplotype&, const Haplotype&);

    void develop(const GenArch&);
    bool checkIndividual(const size_t&);

    Genome genome;
    vecDbl transcriptome;
    vecDbl locivalues;
    bool gender;
    vecDbl genvalues;
    vecDbl traitvalues;
    double ecotrait;
    double matepref;
    double neutrait;
    double fitness;
    vecDbl feedingRates;
    size_t ecotype;

};

#endif


