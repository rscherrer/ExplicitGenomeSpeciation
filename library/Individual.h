#ifndef EXPLICITGENOMESPECIATION_INDIVIDUAL_H
#define EXPLICITGENOMESPECIATION_INDIVIDUAL_H

#include "GenArch.h"
#include "Utilities.h"
#include "Types.h"
#include "Gamete.h"
#include "Random.h"
#include <cassert>
#include <stddef.h>
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<> Genome;

class Individual {

public:

    // To generate an initial population
    Individual(const Param &pars, const GenArch &arch) :
        genome(genomize(pars)),
        transcriptome(utl::zeros(pars.nloci)),
        locivalues(utl::zeros(pars.nloci)),
        genvalues(utl::zeros(3u)),
        traitvalues(utl::zeros(3u)),
        ecotrait(0.0),
        matepref(0.0),
        neutrait(0.0),
        fitness(1.0),
        feeding(utl::zeros(2u)),
        ecotype(0u),
        habitat(0u),
        gender(rnd::bernoulli(0.5)),
        alive(true),
        adult(false)
    {
        develop(pars, arch);

        assert(genome.size() == 2u * pars.nloci);
        assert(transcriptome.size() == pars.nloci);
        assert(traitvalues.size() == 3u);
        assert(fitness > 0.0);
        assert(feeding[0u] > 0.0);
        assert(feeding[1u] > 0.0);
    }


    // Newborn
    Individual(const Param &pars, const GenArch &arch, const Individual &mom,
     const Individual &dad) :
        genome(fecundate(mom, dad, pars, arch)),
        transcriptome(utl::zeros(pars.nloci)),
        locivalues(utl::zeros(pars.nloci)),
        genvalues(utl::zeros(3u)),
        traitvalues(utl::zeros(3u)),
        ecotrait(0.0),
        matepref(0.0),
        neutrait(0.0),
        fitness(1.0),
        feeding(utl::zeros(2u)),
        ecotype(0u),
        habitat(0u),
        gender(rnd::bernoulli(0.5)),
        alive(true),
        adult(false)
    {
        develop(pars, arch);

        assert(genome.size() == 2u * pars.nloci);
        assert(transcriptome.size() == pars.nloci);
        assert(traitvalues.size() == 3u);
        assert(fitness > 0.0);
        assert(feeding[0u] > 0.0);
        assert(feeding[1u] > 0.0);
    }

    ~Individual() {}

    // Life history
    bool isalive() const;
    void disperse();
    void feed(const vecDbl&);
    bool accept(const double&, const Param&) const;
    void survive(const bool&);

    // Getters called from outside
    bool getGender() const
    {
        return gender;
    }
    size_t getEcotype() const
    {
        return ecotype;
    }
    size_t getHabitat() const
    {
        return habitat;
    }
    double getEcoTrait() const
    {
        return ecotrait;
    }
    double getMatePref() const
    {
        return matepref;
    }
    double getNeutral() const
    {
        return neutrait;
    }
    double getFitness() const
    {
        return fitness;
    }
    double getTraitValue(const size_t &trait) const
    {
        return traitvalues[trait];
    }
    double getGenValue(const size_t &trait) const
    {
        return genvalues[trait];
    }
    double getFeeding(const size_t &r) const
    {
        return feeding[r];
    }
    double getLocusValue(const size_t &locus) const
    {
        return locivalues[locus];
    }
    size_t getZygosity(const size_t &locus) const
    {
        const size_t nloci = genome.size() / 2u;
        const size_t zyg = genome.test(locus) + genome.test(locus + nloci);
        assert(zyg == 0u || zyg == 1u || zyg == 2u);
        return zyg;
    }
    size_t getAlleleSum() const
    {
        return genome.count();
    }
    double getExpression() const
    {
        double sum = 0.0;
        for (size_t locus = 0u; locus < transcriptome.size(); ++locus) {
            sum += transcriptome[locus];
        }
        return sum;
    }

    // Force resetters
    void resetEcoTrait(const double &x, const Param &p)
    {
        ecotrait = x;
        traitvalues[0u] = x;
        feeding[0u] = p.maxfeed * exp(-p.ecosel * utl::sqr(ecotrait + 1.0));
        feeding[1u] = p.maxfeed * exp(-p.ecosel * utl::sqr(ecotrait - 1.0));
        assert(feeding[0u] >= 0.0);
        assert(feeding[1u] >= 0.0);
        assert(feeding[0u] <= p.maxfeed);
        assert(feeding[1u] <= p.maxfeed);
    }
    void resetMatePref(const double &y)
    {
        matepref = y;
        traitvalues[1u] = y;
    }
    void resetEcotype(const size_t &e)
    {
        ecotype = e;
    }
    void resetGender(const bool &sex)
    {
        gender = sex;
    }

private:

    // Two ways to generate the genome
    Genome genomize(const Param&) const;
    Genome fecundate(const Individual&, const Individual&, const Param&,
     const GenArch&) const;

    void recombine(Genome&, const Param&, const GenArch&) const;
    void mutate(Genome&, const Param&) const;
    void develop(const Param&, const GenArch&);
    void setFeeding(const size_t&, const double&, const double&);

    Genome genome;
    vecDbl transcriptome;
    vecDbl locivalues;
    vecDbl genvalues;
    vecDbl traitvalues;
    double ecotrait;
    double matepref;
    double neutrait;
    double fitness;
    vecDbl feeding;
    size_t ecotype;
    size_t habitat;
    bool gender;
    bool alive;
    bool adult;

};

#endif


