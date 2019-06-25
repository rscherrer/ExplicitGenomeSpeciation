#ifndef EXPLICITGENOMESPECIATION_INDIVIDUAL_H
#define EXPLICITGENOMESPECIATION_INDIVIDUAL_H

#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include <vector>
#include <random>

class Individual {

public:

    typedef Individual const * PInd;

    struct Locus
    {
        size_t alleleCount;
        double expression;
        double locusGeneticValue;
    };

    // Constructors
    Individual(const ParameterSet&, const GeneticArchitecture&);
    Individual(const std::vector<bool>&, const ParameterSet&, const GeneticArchitecture&);
    Individual(Individual const * const, Individual const * const, const ParameterSet&, const GeneticArchitecture&);

    // Getters
    bool isFemale(const bool&) const;
    double getFitness() const { return fitness;}
    size_t getHabitat() const { return habitat; }
    size_t getEcotype() const { return ecotype; }
    std::pair<double, double> getAttackRates() const { return attackRates; }
    std::vector<size_t> getMates() const { return mates; }
    std::vector<bool> getGenomeSequence() const { return genomeSequence; }
    std::vector<double> getPhenotypes() const { return phenotypes; }
    std::vector<double> getGeneticValues() const { return geneticValues; }
    std::vector<double> getEnvirValues() const { return envirValues; }
    Locus getLocus(const size_t&) const;

private:

    friend class Population;

    // Ecological attributes
    mutable size_t habitat;
    mutable size_t ecotype;
    mutable double fitness;
    mutable double nOffspring;
    mutable std::vector<size_t> mates;
    double matePreference;
    std::pair<double, double> attackRates;

    // Genetic attributes
    bool isHeterogamous;
    std::vector<bool> genomeSequence;
    std::vector<Locus> genotypes;
    std::vector<double> phenotypes;
    std::vector<double> geneticValues;
    std::vector<double> envirValues;

    // Ecology
    void disperse(const size_t& nHabitat) const;
    void setFitness(const std::pair<double, double>&) const;
    void setBurninFitness(const std::pair<double, double>&, const double&) const;
    void setAttackRates(const double&);
    void setMatePreference(const double&);
    void chooseMates(const double&, std::discrete_distribution<size_t>&, const std::vector<PInd>&, const ParameterSet&) const;
    double assessMatingProb(const double&, const double&) const;
    bool acceptMate(Individual const * const, const ParameterSet&) const;
    size_t sampleClutchSize(const double&) const;
    bool survive(const double&) const;
    void setEcotype(const std::pair<double, double>&) const;

    // Genetics
    void mutate(const ParameterSet&);
    void develop(const ParameterSet&, const GeneticArchitecture&);
    void expressGene(const size_t&, const size_t&, const double&, const double&);
    void setAdditiveValue(const size_t&, const double&, const double&);
    void setEpistaticValue(const size_t&, const double&, const std::list<std::pair<size_t, double> >&);
    void setPhenotype(const size_t&);
    void setEnvirValue(const size_t&, const double&);
    void setGeneticValue(const size_t&, const GeneticArchitecture&);
    void setLocusGeneticValue(const size_t&, const GeneticArchitecture&, const ParameterSet&);
    void setGenomeSequence(const size_t&, const double&);
    void recombineFreely(size_t&, size_t&, const size_t&, const double&, double&) const;
    void crossOver(size_t&, const double&, const double&, double&) const;
    void inheritLocus(Individual const * const, const bool&, const size_t&, const size_t&);
    void determineSex(const bool&, const bool&, const size_t&);
    void inheritGamete(Individual const * const, const ParameterSet&, const GeneticArchitecture&);

};

#endif
