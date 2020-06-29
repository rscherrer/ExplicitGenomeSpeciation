#ifndef EXPLICITGENOMESPECIATION_INDIVIDUAL_H
#define EXPLICITGENOMESPECIATION_INDIVIDUAL_H

#include "GenArch.h"
#include "Utilities.h"

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
        transcriptome(std::vector<double>(pars.nloci, 0.0)),
        locivalues(std::vector<double>(pars.nloci, 0.0)),
        genvalues(std::vector<double>(3u, 0.0)),
        traitvalues(std::vector<double>(3u, 0.0)),
        midparents(std::vector<double>(3u, 0.0)),
        fitness(1.0),
        feeding(std::vector<double>(3u, 0.0)),
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
        transcriptome(std::vector<double>(pars.nloci, 0.0)),
        locivalues(std::vector<double>(pars.nloci, 0.0)),
        genvalues(std::vector<double>(3u, 0.0)),
        traitvalues(std::vector<double>(3u, 0.0)),
        midparents(calcmidparent(mom, dad)),
        fitness(1.0),
        feeding(std::vector<double>(2u, 0.0)),
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
    void feed(const std::vector<double>&);
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
    double getFitness() const
    {
        return fitness;
    }
    double getTraitValue(const size_t &trait) const
    {
        return traitvalues[trait];
    }
    double getMidparent(const size_t &trait) const
    {
        return midparents[trait];
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
    size_t getZygosity(const size_t &locus, const size_t &nloci) const
    {
        const size_t zyg = genome.test(locus) + genome.test(locus + nloci);
        assert(zyg == 0u || zyg == 1u || zyg == 2u);
        return zyg;
    }
    unsigned long long getByte(const size_t &B) const
    {
        std::bitset<64u> byte;
        const size_t start = B * 64u;
        size_t end = (B + 1u) * 64u;
        if (end > genome.size()) end = genome.size();
        for (size_t l = start, b = 0u; l < end; ++l, ++b)
            if (genome.test(l)) byte.set(b);
        return byte.to_ullong();
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
    void resetTrait(const size_t &trait, const double &newvalue, const Param &p)
    {
        traitvalues[trait] = newvalue;
        if (trait == 0u) {
            feeding[0u] = exp(-p.ecosel * utl::sqr(traitvalues[trait] + 1.0));
            feeding[1u] = exp(-p.ecosel * utl::sqr(traitvalues[trait] - 1.0));
            assert(feeding[0u] >= 0.0);
            assert(feeding[1u] >= 0.0);
            assert(feeding[0u] <= 1.0);
            assert(feeding[1u] <= 1.0);
        }
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
    std::vector<double> calcmidparent(const Individual&, const Individual&) const;

    Genome genome;
    std::vector<double> transcriptome;
    std::vector<double> locivalues;
    std::vector<double> genvalues;
    std::vector<double> traitvalues;
    std::vector<double> midparents;
    double fitness;
    std::vector<double> feeding;
    size_t ecotype;
    size_t habitat;
    bool gender;
    bool alive;
    bool adult;

};

#endif


