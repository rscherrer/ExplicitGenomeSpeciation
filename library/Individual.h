#ifndef EXPLICITGENOMESPECIATION_INDIVIDUAL_H
#define EXPLICITGENOMESPECIATION_INDIVIDUAL_H

#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include <vector>
#include <random>


/// Function to calculate feeding rates
std::vector<double> calcFeedingRates(const double &sel, const double &trait,
 const double &maxi = 0.0004);

class Individual {

public:

    typedef Individual const * PInd;

    Individual(const Genome&, const MultiNet&);
    ~Individual() {}

    // Getters
    bool getGender() const { return isFemale; }
    double getEcoTrait() const { return ecoTrait; }
    double getMatePref() const { return matePref; }
    double getFitness() const { return fitness; }
    std::vector<double> getFeedingRates() const { return feedingRates; }

    // Actions
    void feed(const std::vector<double>&);
    bool acceptMate(const double&, const double&) const;

    // Setters
    void setEcoTrait(const double &value, const double &sel) {
        ecoTrait = value;
        feedingRates = calcFeedingRates(sel, value);
    }
    void setMatePref(const double &value) { matePref = value; }


    /*
    struct Locus
    {
        size_t alleleCount;
        double expression;
        double locusGeneticValue;
    };
    */

    // Constructors
    // Individual(const ParameterSet&, const GeneticArchitecture&);
    // Individual(const std::vector<bool>&, const ParameterSet&,
    //  const GeneticArchitecture&);
    // Individual(Individual const * const, Individual const * const,
    //  const ParameterSet&, const GeneticArchitecture&);

    // Getters
    // bool isFemale(const bool&) const;
    // double getFitness() const { return fitness;}
    // size_t getHabitat() const { return habitat; }
    // size_t getEcotype() const { return ecotype; }
    // std::pair<double, double> getAttackRates() const { return attackRates; }
    // std::vector<size_t> getMates() const { return mates; }
    // std::vector<bool> getGenomeSequence() const { return genomeSequence; }
    // std::vector<double> getPhenotypes() const { return phenotypes; }
    // std::vector<double> getGeneticValues() const { return geneticValues; }
    // std::vector<double> getEnvirValues() const { return envirValues; }
    // Locus getLocus(const size_t&) const;

    // Function to reset one's habitat
    // void resetHabitat(const size_t&) const;

private:

    friend class Population;

    // Makers
    std::vector<std::vector<bool> > makeSequence(const size_t&);
    std::vector<double> develop(const Genome&, const MultiNet&);

    // Fields
    std::vector<std::vector<bool> > sequence;
    std::vector<double> genexp;
    bool isFemale;
    std::vector<double> traits;
    double ecoTrait;
    double matePref;
    double neutral;
    double fitness;
    std::vector<double> feedingRates;

    // Ecological attributes
    // mutable size_t habitat;
    // mutable size_t ecotype;
    // mutable double fitness;
    // mutable double nOffspring;
    // mutable std::vector<size_t> mates;
    // double matePreference;
    // std::pair<double, double> attackRates;

    // Genetic attributes

    // bool isHeterogamous;
    // std::vector<bool> genomeSequence;
    // std::vector<Locus> genotypes;
    // std::vector<double> phenotypes;
    // std::vector<double> geneticValues;
    // std::vector<double> envirValues;

    // Initialize container sizes
    // void initializeSizeGenotypeVector(const size_t&);
    // void initializeSizeIndivTraitSpecificMetrics(const size_t&);

    // Ecology
    // void disperse(const size_t& nHabitat) const;
    // void setFitness(const std::pair<double, double>&) const;
    // void setBurninFitness(const std::pair<double, double>&, const double&)
    //  const;
    // void setAttackRates(const double&);
    // void setMatePreference(const double&);
    // void chooseMates(const double&, std::discrete_distribution<size_t>&,
    //  const std::vector<PInd>&, const ParameterSet&) const;
    // double assessMatingProb(const double&, const double&) const;
    // bool acceptMate(Individual const * const, const ParameterSet&) const;
    // size_t sampleClutchSize(const double&) const;
    // bool survive(const double&) const;
    // void setEcotype(const std::pair<double, double>&) const;

    // Genetics
    // void mutate(const ParameterSet&);
    // void develop(const ParameterSet&, const GeneticArchitecture&);
    // void expressGene(const size_t&, const size_t&, const double&,
    //  const double&);
    // void setAdditiveValue(const size_t&, const double&, const double&);
    // void setEpistaticValue(const size_t&, const double&,
    //  const std::list<std::pair<size_t, double> >&);
    // void setPhenotype(const size_t&);
    // void setEnvirValue(const size_t&, const double&);
    // void setGeneticValue(const size_t&, const GeneticArchitecture&);
    // void setLocusGeneticValue(const size_t&, const GeneticArchitecture&,
    //  const ParameterSet&);
    // void setGenomeSequence(const size_t&, const double&);
    // void recombineFreely(size_t&, size_t&, const size_t&, const double&,
    //  double&) const;
    // void crossOver(size_t&, const double&, const double&, double&) const;
    // void inheritLocus(Individual const * const, const bool&, const size_t&,
    //  const size_t&);
    // void determineSex(const bool&, const bool&, const size_t&);
    // void inheritGamete(Individual const * const, const ParameterSet&,
    //  const GeneticArchitecture&);

};

#endif
