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

    Individual(const Param&, const GenArch&);
    Individual(const Param&, const GenArch&, const Individual&,
     const Individual&);
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
    size_t getAlleleSum() const;
    double getExpression() const;
    Genome getFullGenome() const;
    std::bitset<64u> getGenomeChunk(const size_t&) const;

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


