#ifndef EXPLICITGENOMESPECIATION_POPULATION_H
#define EXPLICITGENOMESPECIATION_POPULATION_H

#include "ParameterSet.h"
#include "GeneticArchitecture.h"
#include "Random.h"
// #include <list>

class Individual;
class LocusVariables;
typedef Individual const * PInd;

class Population {

public:

    Population(const size_t&);
    ~Population() {};

    // Getters
    size_t getPopSize() const { return popSize; };

    void survive() {};

    // Population(const ParameterSet&, const GeneticArchitecture&);
    // Population(const std::vector<PInd>&, const ParameterSet&);

    // void massExtinction();

    // Getters
    // size_t getNResources() const;
    // std::vector<size_t> getHabitatVector() const;

    // High-level setters
    // void dispersal(const ParameterSet&);
    // void sortByHabitat();
    // void resourceDynamics(const size_t&, const double&);
    // void reproduction(const size_t&, const ParameterSet&,
    //  const GeneticArchitecture&);
    // void survival(const double&);
    // void setBurnin();
    // void endBurnin();
    // void assignEcotypes();
    // void getLocalAttackRates(std::list<std::pair<double, double> >&,
    //  const size_t&);
    // void setEcotypeBoundary(const std::list<std::pair<double, double> >&,
    //  const size_t&);
    // void findEcotypeBoundary(const size_t&);

    // Reset the habitat of an individual
    // void resetIndividualHabitat();

    // Variance decomposition
    // void decomposeVariance(const double&);
    // void initializeVarianceComponents();
    // void accumulateMoments(const size_t&);
    // void completeMoments(const size_t&, const double&);
    // void accumulateSingleLocusContributions();
    // void calcEcotypeDifferentations(const size_t&, const double&);

    // Genome-wide variance decomposition
    // void decomposeVarianceAlongGenome(const double&);


private:

    // Makers
    std::vector<PInd> populate(const size_t&);

    // The population
    std::vector<PInd> individuals;

    // Census
    size_t popSize;

    // std::vector<PInd> females;
    // std::vector<PInd> males;
    // std::vector<PInd> offspring;
    // std::vector<std::pair<size_t, size_t> > genderCounts;
    // std::vector<std::pair<std::list<PInd>::iterator,
    //  std::list<PInd>::iterator> > idHabitatBoundaries;
    // std::vector<size_t> ecotypeSizes;
    // size_t popSize;

    // Mating features
    // std::vector<double> maleSuccesses;
    // std::discrete_distribution<size_t> maleMarket;

    // Ecological state
    // size_t nAccessibleResources;

    // void initializeSizePopEcologicalMetrics(const size_t&);

    // std::vector<std::pair<double, double> > resourceCapacities;
    // std::vector<std::pair<double, double> > replenishRates;
    // std::vector<std::pair<double, double> > resourceConsumption;
    // std::vector<std::pair<double, double> > resourceEql;
    // std::vector<std::pair<double, double> > ecotypeBoundaries;

    // Low-level setters
    // void setResourceCapacities(const double&, const double&);
    // void setReplenishRates(const double&);
    // void setResourceConsumption(const size_t&);
    // void setResourceEquilibrium(const size_t&);
    // void assignFitnesses(const size_t&, const double&);
    // void classifyGenders(const bool&);
    // void setMaleFitnesses(const size_t&, const double&);
    // void birth(const PInd&, const ParameterSet&, const GeneticArchitecture&);
    // void emptyPopulation();

    // Variance components
    // std::vector<double> meanPhenotypes;
    // std::vector<double> meanGeneticValues;
    // std::vector<double> phenotypicVariances;
    // std::vector<double> geneticVariances;
    // std::vector<double> additiveVariances;
    // std::vector<double> dominanceVariances;
    // std::vector<double> interactionVariances;
    // std::vector<double> nonAdditiveVariances;

    // std::vector<double> varianceAlleleFrequencies;
    // std::vector<double> populationExpectedHeterozygosity;

    // std::vector<std::vector<double> > ecotypeMeanGeneticValues;
    // std::vector<std::vector<double> > ecotypePhenotypicVariances;
    // std::vector<std::vector<double> > ecotypeGeneticVariances;
    // std::vector<std::vector<double> > ecotypeAdditiveVariances;
    // std::vector<std::vector<double> > ecotypeNonAdditiveVariances;

    // std::vector<double> Fst;
    // std::vector<double> Pst;
    // std::vector<double> Gst;
    // std::vector<double> Qst;
    // std::vector<double> Cst;

    // std::vector<LocusVariables> locusVariables;

};

/*

// Locus-specific genetic variables
class LocusVariables: public Population {

public:

    // Single-locus variance decomposition
    void decomposeLocusVariance(const double&);

    void initializeLocusVariables();
    void accumulateLocusGeneticMoments();
    void completeLocusGeneticMoments(const double&);
    void accumulateLocusCensus(const size_t&, const size_t&);
    void accumulateLocusAlleleCounts(const size_t&, const size_t&);
    void accumulateLocusGeneticValues(const size_t&, const size_t&,
     const double&);
    void accumulateLocusGeneticValuesByAlleleCounts(const size_t&,
     const double&);
    void calcLocusPhenotypicVariances();
    void regressLocusPhenotypeAgainstGenotype();
    void calcLocusAdditiveVariance();
    void calcLocusDominanceVariance();
    void calcLocusEcotypeAdditiveVariances(const double&);
    void accumulateLocusIndividualResiduals();
    void completeLocusInteractionVariance(const double&);
    void completeLocusNonAdditiveVariances(const double&);
    void calcLocusHeterozygosities(const double&);
    void calcLocusEcotypeDifferentiations(const double&);

    size_t locus;
    size_t trait;

    double locusMeanGeneticValue;
    double locusMeanAlleleCount;
    double locusAlleleFrequency;
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

    double locusPopulationExpectedHeterozygosity;
    double locusEcotypeExpectedHeterozygosity;
    double locusVarianceAlleleFrequencies;

    double locusFst;
    double locusPst;
    double locusGst;
    double locusQst;
    double locusCst;

};

*/

// Accessory functions
// bool compareAlongTradeOff(const std::pair<double, double>&,
//  const std::pair<double, double>&);
// double calcLogisticResourceEq(const double&, const double&, const double&);
// double Xst(const double&, const std::vector<double>&, const size_t&, const std::vector<size_t>&, const double&);
// void sum2mean(double&, const size_t&);
// void sumsq2var(double&, const size_t&, const double&, const double&);
// void sumprod2cov(double&, const size_t&, const double&, const double&, const double&);
// void clipDown(double&, const double&, const double& = 0.0);

#endif //EXPLICITGENOMESPECIATION_POPULATION_H
