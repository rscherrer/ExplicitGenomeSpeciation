#ifndef EXPLICITGENOMESPECIATION_INDIVIDUAL_H
#define EXPLICITGENOMESPECIATION_INDIVIDUAL_H

#include "GenArch.h"
#include "Utilities.h"
#include "Types.h"
#include "Random.h"
#include <cassert>
#include <stddef.h>
#include <bitset>

typedef std::bitset<10000> Genome;

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
        ecomidparent(0.0),
        matmidparent(0.0),
        neumidparent(0.0),
        fitness(1.0),
        feeding(utl::zeros(2u)),
        ecotype(0u),
        habitat(0u),
        gender(determinesex()),
        alive(true),
        adult(true)
    {
        develop(pars, arch);

        assert(transcriptome.size() == pars.nloci);
        assert(traitvalues.size() == 3u);
        assert(fitness >= 0.0);
        assert(feeding[0u] >= 0.0);
        assert(feeding[1u] >= 0.0);
        assert(feeding[0u] <= 1.0);
        assert(feeding[1u] <= 1.0);
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
        ecomidparent((mom.getEcoTrait() + dad.getEcoTrait()) / 2.0),
        matmidparent((mom.getMatePref() + dad.getMatePref()) / 2.0),
        neumidparent((mom.getNeutral() + dad.getNeutral()) / 2.0),
        fitness(1.0),
        feeding(utl::zeros(2u)),
        ecotype(0u),
        habitat(mom.getHabitat()),
        gender(determinesex()),
        alive(true),
        adult(false)
    {
        develop(pars, arch);

        assert(transcriptome.size() == pars.nloci);
        assert(traitvalues.size() == 3u);
        assert(fitness >= 0.0);
        assert(feeding[0u] >= 0.0);
        assert(feeding[1u] >= 0.0);
        assert(feeding[0u] <= 1.0);
        assert(feeding[1u] <= 1.0);
    }

    ~Individual() {}

    // Life history
    bool isalive() const;
    void disperse();
    void feed(const vecDbl&);
    double mate(const double&, const Param&) const;
    void survive(const bool&);
    void classify(const double&);

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
    double getEcoMidparent() const
    {
        return ecomidparent;
    }
    double getMatMidparent() const
    {
        return matmidparent;
    }
    double getNeuMidparent() const
    {
        return neumidparent;
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
        feeding[0u] = exp(-p.ecosel * utl::sqr(ecotrait + 1.0));
        feeding[1u] = exp(-p.ecosel * utl::sqr(ecotrait - 1.0));
        assert(feeding[0u] >= 0.0);
        assert(feeding[1u] >= 0.0);
        assert(feeding[0u] <= 1.0);
        assert(feeding[1u] <= 1.0);
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

    bool determinesex() const;

    Genome genome;
    vecDbl transcriptome;
    vecDbl locivalues;
    vecDbl genvalues;
    vecDbl traitvalues;
    double ecotrait;
    double matepref;
    double neutrait;
    double ecomidparent;
    double matmidparent;
    double neumidparent;
    double fitness;
    vecDbl feeding;
    size_t ecotype;
    size_t habitat;
    bool gender;
    bool alive;
    bool adult;

};

#endif


