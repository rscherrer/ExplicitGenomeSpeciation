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

    // Various getters
    bool getGender() const;
    size_t getEcotype() const;
    size_t getHabitat() const;
    double getFitness() const;
    double getTraitValue(const size_t&) const;
    double getMidparent(const size_t&) const;
    double getGenValue(const size_t&) const;
    double getFeeding(const size_t&) const;
    double getLocusValue(const size_t&) const;
    size_t getZygosity(const size_t&, const size_t&) const;
    unsigned long long getByte(const size_t&) const;
    size_t getAlleleSum() const;
    double getExpression() const;

    // Force resetters (for testing purposes)
    void resetTrait(const size_t&, const double&, const Param&);
    void resetEcotype(const size_t&);
    void resetGender(const bool&);

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


