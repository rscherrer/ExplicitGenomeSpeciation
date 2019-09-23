#ifndef EXPLICITGENOMESPECIATION_INDIVIDUAL_H
#define EXPLICITGENOMESPECIATION_INDIVIDUAL_H

#include "Param.h"
#include "GenArch.h"
#include "Utilities.h"
#include "Types.h"
#include <vector>
#include <random>
#include <stddef.h>

typedef std::vector<bool> Haplotype;
typedef std::vector<Haplotype > Diplotype;

class Individual {

    friend class Population;
    friend class MetaPop;

    typedef Individual const * PInd;

public:

    /// Spontaneous creation
    Individual(const GenArch &arch, const double &ecosel,
     const double &maxfeeding, const double &snpFreq) :
        genome(makeSequence(arch, snpFreq)),
        transcriptome(utl::zeros(arch.nLoci)),
        locivalues(utl::zeros(arch.nLoci)),
        isFemale(determineSex(arch.femHeterogamy)),
        genvalues(utl::zeros(3u)),
        traitvalues(utl::zeros(3u)),
        ecotrait(traitvalues[0u]),
        matepref(traitvalues[1u]),
        neutral(traitvalues[2u]),
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
        isFemale(determineSex(arch.femHeterogamy)),
        genvalues(utl::zeros(3u)),
        traitvalues(utl::zeros(3u)),
        ecotrait(traitvalues[0u]),
        matepref(traitvalues[1u]),
        neutral(traitvalues[2u]),
        fitness(1.0),
        feedingRates(calcFeedingRates(ecosel, ecotrait, maxfeeding)),
        ecotype(0u)
    {
        develop(arch);
        checkIndividual(arch.nLoci);
    }

    ~Individual() {}

    bool getGender() const { return isFemale; }
    double getEcoTrait() const { return ecotrait; }
    double getMatePref() const { return matepref; }
    double getNeutral() const { return neutral; }
    double getFitness() const { return fitness; }
    vecDbl getTraits() const { return traitvalues; }
    vecDbl getFeedingRates() const { return feedingRates; }
    vecDbl getGenValues() const { return genvalues; }
    Haplotype getSequence(const size_t &i) const { return genome[i]; }
    vecDbl getExpression() const { return transcriptome; }
    size_t getEcotype() const { return ecotype; }
    size_t getZygosity(const size_t&);
    double getLocusValue(const size_t&);
    bool acceptMate(const double&, const double&) const;
    Haplotype recombine(const GenArch&);
    vecDbl calcFeedingRates(const double&, const double&, const double&);

    void setEcoTrait(const double&, const double&, const double&);
    void setMatePref(const double&);
    void setEcotype(const double&);
    void setGender(const bool&);
    void feed(const vecDbl&);
    void mutate(Haplotype&, const double& = 1.0e-5);

private:

    Diplotype makeSequence(const GenArch&, double = -1.0);
    Diplotype fecundate(const Haplotype&, const Haplotype&);
    bool determineSex(const bool&);

    void develop(const GenArch&);
    bool checkIndividual(const size_t&);

    Diplotype genome;
    vecDbl transcriptome;
    vecDbl locivalues;
    bool isFemale;
    vecDbl genvalues;
    vecDbl traitvalues;
    double ecotrait;
    double matepref;
    double neutral;    
    double fitness;
    vecDbl feedingRates;
    size_t ecotype;

};

#endif


