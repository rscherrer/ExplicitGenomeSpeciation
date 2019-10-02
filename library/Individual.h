#ifndef EXPLICITGENOMESPECIATION_INDIVIDUAL_H
#define EXPLICITGENOMESPECIATION_INDIVIDUAL_H

//#include "Param.h"
#include "GenArch.h"
#include "Utilities.h"
#include "Types.h"
#include "Gamete.h"
//#include "Deme.h"
#include "Random.h"
//#include <vector>
//#include <random>
//#include <cmath> // exponential?
#include <cassert>
//#include <algorithm>
#include <stddef.h>

typedef std::vector<Haplotype> Genome;


class Individual {

public:

    /// Spontaneous creation
    Individual(const GenArch &arch, const double &ecosel,
     const double &maxfeed, const double &snpFreq) :
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
        feeding(utl::zeros(2u)),
        ecotype(0u)
    {
        develop(arch, ecosel, maxfeed);
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
        feeding(utl::zeros(2u)),
        ecotype(0u)
    {
        develop(arch, ecosel, maxfeeding);
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
    double getFeeding(const size_t &r) const { return feeding[r]; }
    double getExpression(const size_t &l) const { return transcriptome[l]; }
    double getExpression() const;
    double getLocusValue(const size_t&) const;
    size_t getEcotype() const { return ecotype; }
    size_t getZygosity(const size_t&) const;
    size_t getAlleleSum(const size_t&) const;

    bool acceptMate(const double&, const double&) const;
    Haplotype recombine(const GenArch&) const;

    void feed(const vecDbl&);

    void setEcoTrait(const double&, const double&, const double&);
    void setMatePref(const double&);
    void setEcotype(const size_t &e) { ecotype = e; }
    void setGender(const bool&);

private:

    Genome generateGenome(const GenArch&, double = -1.0);
    Genome fecundate(const Haplotype&, const Haplotype&);

    void develop(const GenArch&, const double&, const double&);
    void setFeeding(const size_t&, const double&, const double&);

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
    vecDbl feeding;
    size_t ecotype;

};

#endif


