//
// Created by p278834 on 7-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_POPULATION_H
#define EXPLICITGENOMESPECIATION_POPULATION_H

#include <vector>
#include <set>
#include <random>
#include "Individual.h"
#include "ParameterSet.h"
#include "GeneticArchitecture.h"

// Forward declaration
class Individual;
typedef Individual const * PInd;

class Population {

public:

    Population(const ParameterSet&, const GeneticArchitecture&);

    // Getters
    size_t getPopSize() const;
    size_t getNResources() const;
    size_t getEcotypeSize(const size_t&) const;

    // High-level member functions
    void dispersal(const ParameterSet&);
    void sortByHabitat();
    void resourceDynamics(const size_t&, const double&);
    void reproduction(const size_t&, const ParameterSet&, const GeneticArchitecture&);
    void survival(const double&);
    void setBurnin();
    void endBurnin();

    void assignEcotypes();
    void getLocalAttackRates(std::list<std::pair<double, double> >&, const size_t&);
    void setEcotypeBoundary(const std::list<std::pair<double, double> >&, const size_t&);
    void findEcotypeBoundary(const size_t&);

    // Variance decomposition
    void decomposeVariance(const double&);
    void initializeVarianceComponents();
    void accumulateMoments(const size_t&);
    void completeMoments(const size_t&);
    void accumulateSingleLocusContributions();
    void calcEcotypeDifferentations(const size_t&, const double&);
    
    void decomposeVarianceAlongGenome();


protected:

    // The population
    std::list<PInd> individuals;
    std::vector<PInd> females;
    std::vector<PInd> males;
    std::vector<PInd> offspring;
    std::vector<std::pair<size_t, size_t> > genderCounts;
    std::vector<std::pair<std::list<PInd>::iterator, std::list<PInd>::iterator> > idHabitatBoundaries;
    std::vector<size_t> ecotypeSizes;
    size_t popSize;

    // Mating features
    std::vector<double> maleSuccesses;
    std::discrete_distribution<size_t> maleMarket;

    // Ecological state
    size_t nAccessibleResources;
    std::vector<std::pair<double, double> > resourceCapacities;
    std::vector<std::pair<double, double> > replenishRates;
    std::vector<std::pair<double, double> > resourceConsumption;
    std::vector<std::pair<double, double> > resourceEql;
    std::vector<std::pair<double, double> > ecotypeBoundaries;

    // Low-level member functions
    void setResourceCapacities(const double&, const double&);
    void setReplenishRates(const double&);
    void setResourceConsumption(const size_t&);
    void setResourceEquilibrium(const size_t&);
    void assignFitnesses(const size_t&, const double&);
    void classifyGenders(const bool&);
    void setMaleFitnesses(const size_t&, const double&);
    void birth(const PInd&, const ParameterSet&, const GeneticArchitecture&);
    void emptyPopulation();

    // Variance components
    std::vector<double> meanPhenotypes;
    std::vector<double> meanGeneticValues;
    std::vector<double> phenotypicVariances;
    std::vector<double> geneticVariances;
    std::vector<double> additiveVariances;
    std::vector<double> dominanceVariances;
    std::vector<double> interactionVariances;
    std::vector<double> nonAdditiveVariances;

    std::vector<std::vector<double> > ecotypeMeanGeneticValues;
    std::vector<std::vector<double> > ecotypePhenotypicVariances;
    std::vector<std::vector<double> > ecotypeGeneticVariances;
    std::vector<std::vector<double> > ecotypeAdditiveVariances;
    std::vector<std::vector<double> > ecotypeNonAdditiveVariances;

    std::vector<double> Pst;
    std::vector<double> Gst;
    std::vector<double> Qst;
    std::vector<double> Cst;

    std::vector<LocusVariables> locusVariables;

};

// Locus-specific genetic variables
class LocusVariables: public Population {

public:

    // Single-locus variance decomposition
    void decomposeLocusVariance();

    void initializeLocusVariables();
    void accumulateLocusGeneticMoments();
    void completeLocusGeneticMoments();
    void calcLocusPhenotypicVariances();
    void regressLocusPhenotypeAgainstGenotype();
    void calcLocusAdditiveVariance();
    void calcLocusDominanceVariance();
    void calcLocusEcotypeAdditiveVariances();
    void completeLocusInteractionVariance();
    void completeLocusNonAdditiveVariances();
    void calcLocusHeterozygosities();
    void calcLocusEcotypeDifferentiations();

    size_t locus;
    size_t trait;

    double locusMeanGeneticValue;
    double locusMeanAlleleCount;
    double locusVarAlleleCount;
    double locusCovGeneticValueAlleleCount;
    double locusAvgSubstitutionEffect;

    std::vector<size_t> locusGenotypeSizes;
    std::vector<double> locusGenotypeBreedingValues;
    std::vector<double> locusGenotypeAdditiveExpectations;
    std::vector<double> locusGenotypeMeanGeneticValues;
    std::vector<double> locusGenotypeDominanceDeviations;

    std::vector<std::vector<size_t> > locusGenotypeEcotypeSizes;

    double locusGeneticVariance;
    double locusPhenotypicVariance;
    double locusAdditiveVariance;
    double locusDominanceVariance;
    double locusInteractionVariance;
    double locusNonAdditiveVariance;
    double locusEnvirVariance;

    std::vector<double> locusEcotypeMeanGeneticValues;
    std::vector<double> locusEcotypeMeanBreedingValues;
    std::vector<double> locusEcotypeMeanNonAdditiveDeviations;
    std::vector<double> locusEcotypePhenotypicVariances;
    std::vector<double> locusEcotypeGeneticVariances;
    std::vector<double> locusEcotypeAdditiveVariances;
    std::vector<double> locusEcotypeNonAdditiveVariances;
    std::vector<double> locusEcotypeAlleleFrequencies;

    double locusExpectedHeterozygosity;
    double locusObservedHeterozygosity;

    double locusFst;
    double locusPst;
    double locusGst;
    double locusQst;
    double locusCst;

};

bool compareAlongTradeOff(const std::pair<double, double>&, const std::pair<double, double>&);

#endif //EXPLICITGENOMESPECIATION_POPULATION_H
