//
// Created by p278834 on 7-5-2019.
//

#include "Population.h"
#include "random.h"
#include "square.h"
#include <cassert>
#include <vector>

double Xst(const double &popVar, const std::vector<double> &groupVars, const size_t &popSize, const std::vector<size_t> &groupSizes, const double& tiny)
{
    if (popVar < tiny) {
        return 0.0;
    }
    else {

        // Compare observed and expected variance between groups
        double groupVar = (groupSizes[0u] * groupVars[0u] + groupSizes[1u] * groupVars[1u]) / popSize;
        double Xst = 1.0 - groupVar / popVar;

        // Population structure statistic is bounded between zero and one
        if (Xst < tiny) {
            Xst = 0.0;
        }
        if (Xst > 1.0 - tiny) {
            Xst = 1.0;
        }
        return Xst;
    }
}

// Constructor

Population::Population(const ParameterSet &parameters, const GeneticArchitecture &geneticArchitecture)
{

    // With genetic sequence provided
    if (parameters.sequence.size() == parameters.nBits) {
        for (size_t i = 0u; i < parameters.nIndividualInit; ++i) {
            individuals.push_back(new Individual(parameters.sequence, parameters, geneticArchitecture));
        }
    }

    // Or not
    else {
        for (size_t i = 0u; i < parameters.nIndividualInit; ++i) {
            individuals.push_back(new Individual(parameters, geneticArchitecture));
        }
    }

    // Initialize resource growth parameters
    setResourceCapacities(parameters.maxResourceCapacity, parameters.habitatAsymmetry);
    setReplenishRates(parameters.maxResourceGrowth);
}


// Getters

size_t Population::getPopSize() const
{
    size_t popSize = individuals.size();
    return popSize;
}

size_t Population::getNResources() const
{
    return nAccessibleResources;
}

size_t Population::getEcotypeSize(const size_t &ecotype) const
{
    return ecotypeSizes[ecotype];
}

// High-level functions

void Population::dispersal(const ParameterSet& parameters)
{

    // Loop through the whole population if dispersal is high, all individuals have high chances to migrate
    if (parameters.dispersalRate > 0.5) {
        for (PInd pInd : individuals) {
            if (rnd::bernoulli(parameters.dispersalRate)) {
                pInd->disperse(parameters.nHabitat);
            }
        }
    }

    // Otherwise sample first the number of migrants
    else {
        const size_t populationSize = individuals.size();
        size_t nMigrants = rnd::binomial(populationSize, parameters.dispersalRate);
        if (nMigrants == 0u) {
            return;
        }
        std::set<size_t> migrants;
        while (migrants.size() < nMigrants) {
            migrants.insert(rnd::random_int(populationSize));
        }
        auto itInd = individuals.cbegin();  // An iterator pointing to an individual pointer (this is a pointerCeption)
        size_t j = 0u;
        for (size_t i : migrants) {
            std::advance(itInd, i - j);
            (*itInd)->disperse(parameters.nHabitat);
            j = i;
        }
    }
}

void Population::sortByHabitat()
{
    const size_t habitat = 0u;
    auto iti = individuals.begin();

    for (auto itj = individuals.end(); iti != itj;) {

        // If the individual is from the wrong habitat
        if ((*iti)->getHabitat() != habitat) {

            // Move individuals from the wrong habitat to the end
            --itj;
            std::swap(*iti, *itj);
        }
        else {
            ++iti;
        }
    }

    // Get the index of the first individual from the second habitat
    const size_t lastSortedHabitat = (*iti)->getHabitat();
    if (lastSortedHabitat == 0u) {
        ++iti;
    }

    // Get the limit inhabitants of each habitat
    idHabitatBoundaries[0u].first = individuals.begin();
    idHabitatBoundaries[0u].second = iti.operator--();
    idHabitatBoundaries[1u].first = iti;
    idHabitatBoundaries[1u].second = individuals.end();
}

void Population::resourceDynamics(const size_t &habitat, const double &ecoSelCoeff)
{

    // Calculate total resource consumption in the current habitat
    setResourceConsumption(habitat);

    // Calculate equilibrium resource concentrations
    setResourceEquilibrium(habitat);

    // Assign individual fitnesses
    assignFitnesses(habitat, ecoSelCoeff);

}

void Population::reproduction(const size_t &habitat, const ParameterSet &parameters, const GeneticArchitecture &geneticArchitecture)
{

    // Gender classification
    classifyGenders(parameters.isFemaleHeterogamy);

    // Gender counts
    const size_t nFemales = genderCounts[habitat].first = females.size();
    const size_t nMales = genderCounts[habitat].second = males.size();

    // Terminate if there are no males or no females
    if (nFemales == 0u || nMales == 0u) {
        return;
    }

    // Male fitnesses
    setMaleFitnesses(nMales, parameters.tiny);

    // Male market
    std::discrete_distribution<size_t> maleMarket(maleSuccesses.begin(), maleSuccesses.end());

    // Set the length of the mating season
    const size_t matingSeasonEnd = rnd::geometric(parameters.mateEvaluationCost);

    // Mate choice and reproduction
    for (PInd female : females) {

        // Sample family size for female
        female->sampleClutchSize(parameters.birthRate);

        // Find fathers for her babies
        female->chooseMates(matingSeasonEnd, maleMarket, males);

        // Babies are born
        birth(female, parameters, geneticArchitecture);
    }
}

void Population::survival(const double &survivalProb)
{
    emptyPopulation();
    for (PInd pInd : males) {
        if (pInd->survive(survivalProb)) {
            individuals.push_back(pInd);
        }
    }
    for (PInd pInd : females) {
        if (pInd->survive(survivalProb)) {
            individuals.push_back(pInd);
        }
    }
    for (PInd pInd : offspring) {
        individuals.push_back(pInd);
    }
}

void Population::setBurnin()
{
    nAccessibleResources = 1u;
}

void Population::endBurnin()
{
    nAccessibleResources = 2u;
}

bool compareAlongTradeOff (const std::pair<double, double> &x, const std::pair<double, double> &y)
{
    bool yOnLeft = y.first < y.second;
    if(x.first < x.second) {
        if(yOnLeft) return (x.first < y.first);
        else return true;
    }
    else {
        if(yOnLeft) return false;
        else return (x.second < y.second);
    }
}

void Population::getLocalAttackRates(std::list<std::pair<double, double> > &listAttackRates, const size_t &habitat)
{
    for (auto itInd = idHabitatBoundaries[habitat].first;; ++itInd) {
        listAttackRates.push_back((*itInd)->getAttackRates());
        if (itInd == idHabitatBoundaries[habitat].second)
        {
            break;
        }
    }
}

void Population::setEcotypeBoundary(const std::list<std::pair<double, double> > &listAttackRates,
        const size_t &habitat)
{
    ecotypeBoundaries[habitat] = listAttackRates.back();
    std::pair<double, double> payoff;
    for (auto attackRates : listAttackRates) {

        // Record payoffs from feeding on the two resources
        payoff.first = attackRates.first * resourceEql[habitat].first;
        payoff.second = attackRates.second * resourceEql[habitat].second;

        // When the first resource becomes more advantageous, the boundary has been reached
        if (payoff.second < payoff.first) {
            ecotypeBoundaries[habitat] = attackRates;
            break;
        }
    }
}

void Population::findEcotypeBoundary(const size_t &habitat)
{

    // Get attack rates
    std::list<std::pair<double, double> > listAttackRates;
    getLocalAttackRates(listAttackRates, habitat);

    // Sort attack rates
    listAttackRates.sort(compareAlongTradeOff);

    // Loop through individuals to find the ecotype boundary
    setEcotypeBoundary(listAttackRates, habitat);

}

void Population::assignEcotypes()
{

    // Find ecotype boundaries within each habitat
    findEcotypeBoundary(0u);
    findEcotypeBoundary(1u);

    // Set ecotypes
    for (PInd pInd : individuals) {
        pInd->setEcotype(ecotypeBoundaries[pInd->getHabitat()]);

        // Update census within each ecotype
        ++ecotypeSizes[pInd->getEcotype()];
    }
}

void sum2mean(double &mean, const size_t &nobs)
{
    mean /= nobs;
}

void sumsq2var(double &variance, const size_t &nobs, const double &mean)
{
    variance -= nobs * sqr(mean);
    variance /= nobs;
}

void sumprod2cov(double &covariance, const size_t &nobs, const double &firstmean, const double &secondmean)
{
    covariance -= nobs * firstmean * secondmean;
    covariance /= nobs;
}

void clipDown(double &value, const double &tiny, const double &lowerbound = 0.0)
{
    value = value > tiny ? value : lowerbound;
}

void Population::decomposeVarianceAlongGenome()
{
    for (LocusVariables * locus : locusVariables) {
        locus->decomposeLocusVariance();
    }
}

void Population::initializeVarianceComponents()
{
    additiveVariances = {0.0, 0.0, 0.0};
    dominanceVariances = {0.0, 0.0, 0.0};
    interactionVariances = {0.0, 0.0, 0.0};
    nonAdditiveVariances = {0.0, 0.0, 0.0};
    ecotypeAdditiveVariances = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    ecotypeNonAdditiveVariances = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};

    // meanPhenotypes = {0.0, 0.0, 0.0};
    meanGeneticValues = {0.0, 0.0, 0.0};
    phenotypicVariances = {0.0, 0.0, 0.0};
    geneticVariances = {0.0, 0.0, 0.0};
    ecotypeMeanGeneticValues = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    ecotypePhenotypicVariances = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
    ecotypeGeneticVariances = {{0.0, 0.0}, {0.0, 0.0}, {0.0, 0.0}};
}

void Population::accumulateSingleLocusContributions()
{
    // Accumulate other variance components from single locus contributions
    size_t trait;
    for (LocusVariables * locus :  locusVariables) {
        trait = locus->trait;
        additiveVariances[trait] += locus->locusAdditiveVariance;
        dominanceVariances[trait] += locus->locusDominanceVariance;
        interactionVariances[trait] += locus->locusInteractionVariance;
        nonAdditiveVariances[trait] += locus->locusNonAdditiveVariance;

        for (size_t ecotype = 0u; ecotype < 2u; ++ecotype) {
            ecotypeAdditiveVariances[trait][ecotype] += locus->locusEcotypeAdditiveVariances[ecotype];
            ecotypeNonAdditiveVariances[trait][ecotype] += locus->locusEcotypeNonAdditiveVariances[ecotype];
        }
    }
}

void Population::accumulateMoments(const size_t &trait)
{
    double phenotype;
    double geneticValue;
    size_t ecotype;

    // Compute genetic and phenotypic variance
    for (PInd pInd : individuals) {
        ecotype = pInd->getEcotype();
        phenotype = pInd->getPhenotypes()[trait];
        geneticValue = pInd->getGeneticValues()[trait];

        // meanPhenotypes[trait] += phenotype;  // Approximately equal to meanGeneticValues because E(environmental effects) = 0
        meanGeneticValues[trait] += geneticValue;
        phenotypicVariances[trait] += sqr(phenotype);
        geneticVariances[trait] += sqr(geneticValue);

        ecotypeMeanGeneticValues[trait][ecotype] += geneticValue;
        ecotypePhenotypicVariances[trait][ecotype] += sqr(phenotype);
        ecotypeGeneticVariances[trait][ecotype] += sqr(geneticValue);
    }
}

void Population::completeMoments(const size_t &trait)
{
    // Complete moments
    // sum2mean(meanPhenotypes[trait], popSize);
    sum2mean(meanGeneticValues[trait], popSize);
    sumsq2var(phenotypicVariances[trait], popSize, meanGeneticValues[trait]);
    sumsq2var(geneticVariances[trait], popSize, meanGeneticValues[trait]);
    for (size_t ecotype = 0u; ecotype < 2u; ++ecotype) {
        sum2mean(ecotypeMeanGeneticValues[trait][ecotype], ecotypeSizes[ecotype]);
        sumsq2var(ecotypePhenotypicVariances[trait][ecotype], ecotypeSizes[ecotype], ecotypeMeanGeneticValues[trait][ecotype]);
        sumsq2var(ecotypeGeneticVariances[trait][ecotype], ecotypeSizes[ecotype], ecotypeMeanGeneticValues[trait][ecotype]);
    }
}



void Population::calcEcotypeDifferentations(const size_t &trait, const double &tiny)
{
    Pst[trait] = Xst(phenotypicVariances, ecotypePhenotypicVariances, popSize, ecotypeSizes, tiny);
    Gst[trait] = Xst(geneticVariances, ecotypeGeneticVariances, popSize, ecotypeSizes, tiny);
    Qst[trait] = Xst(additiveVariances, ecotypeAdditiveVariances, popSize, ecotypeSizes, tiny);
    Cst[trait] = Xst(nonAdditiveVariances, ecotypeNonAdditiveVariances, popSize, ecotypeSizes, tiny);
}

void Population::decomposeVariance(const double &tiny)
{

    initializeVarianceComponents();

    for (size_t trait = 0u; trait < 3u; ++trait) {
        accumulateMoments(trait);
        completeMoments(trait);
    }

    accumulateSingleLocusContributions();

    // Calculate differentiation
    for (size_t trait = 0u; trait < 3u; ++trait) {
        calcEcotypeDifferentations(trait, tiny);
    }

}




// Accessory functions

double calcLogisticResourceEq(const double &resourceCapacity, const double &replenishRate, const double &consumption)
{
    double resource = resourceCapacity * (1.0 - consumption / replenishRate);
    return resource;
}


// Single-locus variance decomposition

void LocusVariables::decomposeLocusVariance()
{

    initializeLocusVariables();
    accumulateLocusGeneticMoments();
    completeLocusGeneticMoments();
    calcLocusPhenotypicVariances();
    regressLocusPhenotypeAgainstGenotype();
    calcLocusAdditiveVariance();
    calcLocusDominanceVariance();
    calcLocusEcotypeAdditiveVariances();
    completeLocusInteractionVariance();
    completeLocusNonAdditiveVariances();
    calcLocusHeterozygosities();
    calcLocusEcotypeDifferentiations();

}

void LocusVariables::initializeLocusVariables()
{
    // Initialize components
    locusGenotypeSizes = {0u, 0u, 0u};
    locusGenotypeEcotypeSizes = {{0u, 0u}, {0u, 0u}, {0u, 0u}};
    locusMeanGeneticValue = 0.0;
    locusGeneticVariance = 0.0;
    locusMeanAlleleCount = 0.0;
    locusVarAlleleCount = 0.0;
    locusCovGeneticValueAlleleCount = 0.0;
    locusGenotypeMeanGeneticValues = {0.0, 0.0, 0.0};
    locusEcotypeMeanGeneticValues = {0.0, 0.0};
    locusEcotypeGeneticVariances = {0.0, 0.0};
    locusEcotypeAlleleFrequencies = {0.0, 0.0};
    locusDominanceVariance = 0.0;
    locusEcotypeAdditiveVariances = {0.0, 0.0};
    locusEcotypeMeanBreedingValues = {0.0, 0.0};
    locusInteractionVariance = 0.0;
    locusNonAdditiveVariance = 0.0;
    locusEcotypeNonAdditiveVariances = {0.0, 0.0};
    locusEcotypeMeanNonAdditiveDeviations = {0.0, 0.0};
    locusObservedHeterozygosity = 0.0;
}

void LocusVariables::accumulateLocusCensus(const size_t &genotype, const size_t &ecotype)
{
    // Accumulate genotype group sizes
    ++locusGenotypeSizes[genotype];
    ++locusGenotypeEcotypeSizes[genotype][ecotype];
}

void LocusVariables::accumulateLocusAlleleCounts(const size_t &genotype, const size_t &ecotype)
{
    // Accumulate allele frequencies
    locusMeanAlleleCount += genotype;
    locusEcotypeAlleleFrequencies[ecotype] += genotype;
    locusVarAlleleCount += sqr(genotype);
}

void LocusVariables::accumulateLocusGeneticValues(const size_t &genotype, const size_t &ecotype, const double &geneticValue)
{
    // Accumulate genetic values and their squares for mean and variances
    locusMeanGeneticValue += geneticValue;
    locusEcotypeMeanGeneticValues[ecotype] += geneticValue;
    locusGenotypeMeanGeneticValues[genotype] += geneticValue;
    locusGeneticVariance += sqr(geneticValue);
    locusEcotypeGeneticVariances[ecotype] += sqr(geneticValue);
}

void LocusVariables::accumulateLocusGeneticValuesByAlleleCounts(const size_t &genotype, const double &geneticValue)
{
    // Accumulate product for covariance
    locusCovGeneticValueAlleleCount += genotype * geneticValue;
}

void LocusVariables::calcLocusHeterozygosities()
{
    for(size_t ecotype = 0u; ecotype < 2u; ++ecotype) {
        locusObservedHeterozygosity += ecotypeSizes[ecotype] * 2.0 * locusEcotypeAlleleFrequencies[ecotype] * (1.0 - locusEcotypeAlleleFrequencies[ecotype]);
    }
    locusObservedHeterozygosity /= popSize;
    locusExpectedHeterozygosity = locusMeanAlleleCount * (1.0 - 0.5 * locusMeanAlleleCount);
}

void LocusVariables::completeLocusGeneticMoments()
{
    // Complete population moments
    sum2mean(locusMeanGeneticValue, popSize);
    sum2mean(locusMeanAlleleCount, popSize);
    sumsq2var(locusGeneticVariance, popSize, locusMeanGeneticValue);
    sumsq2var(locusVarAlleleCount, popSize, locusMeanAlleleCount);
    sumprod2cov(locusCovGeneticValueAlleleCount, popSize, locusMeanGeneticValue, locusMeanAlleleCount);

    // Complete ecotype moments
    for (size_t ecotype = 0u; ecotype < 2u; ++ecotype) {
        sum2mean(locusEcotypeMeanGeneticValues[ecotype], ecotypeSizes[ecotype]);
        sum2mean(locusEcotypeAlleleFrequencies[ecotype], static_cast<size_t>(2.0 * ecotypeSizes[ecotype]));
        sumsq2var(locusEcotypeGeneticVariances[ecotype], ecotypeSizes[ecotype], locusEcotypeMeanGeneticValues[ecotype]);
    }

    // Complete genotype moments
    for (size_t genotype = 0u; genotype < 3u; ++genotype) {
        sum2mean(locusGenotypeMeanGeneticValues[genotype], locusGenotypeSizes[genotype]);
    }
}

void LocusVariables::calcLocusPhenotypicVariances()
{
    locusPhenotypicVariance = locusGeneticVariance + locusEnvirVariance;
    for (size_t ecotype = 0u; ecotype < 2u; ++ecotype) {
        locusEcotypePhenotypicVariances[ecotype] = locusEcotypeGeneticVariances[ecotype] + locusEnvirVariance;
    }
}

void LocusVariables::regressLocusPhenotypeAgainstGenotype()
{
    // Regression analysis at the whole population level
    locusAvgSubstitutionEffect = locusCovGeneticValueAlleleCount / locusVarAlleleCount;
    for (size_t genotype = 0u; genotype < 3u; ++genotype) {
        locusGenotypeBreedingValues[genotype] = locusAvgSubstitutionEffect * (genotype - locusMeanAlleleCount);
        locusGenotypeAdditiveExpectations[genotype] = locusGenotypeMeanGeneticValues[genotype] + locusGenotypeBreedingValues[genotype];
        locusGenotypeDominanceDeviations[genotype] = locusGenotypeMeanGeneticValues[genotype] - locusGenotypeAdditiveExpectations[genotype];
    }
}

void LocusVariables::calcLocusAdditiveVariance()
{
    locusAdditiveVariance = locusAvgSubstitutionEffect * locusVarAlleleCount;
}

void LocusVariables::calcLocusDominanceVariance()
{
    for (size_t genotype = 0u; genotype < 3u; ++genotype) {
        locusDominanceVariance += locusGenotypeSizes[genotype] * sqr(locusGenotypeDominanceDeviations[genotype]);
    }
    locusDominanceVariance /= popSize;
}

void LocusVariables::calcLocusEcotypeAdditiveVariances()
{
    // Calculate within-ecotype additive variance
    for (size_t ecotype = 0u; ecotype < 2u; ++ecotype) {
        for (size_t genotype = 0u; genotype < 3u; ++genotype) {
            locusEcotypeMeanBreedingValues[ecotype] += locusGenotypeEcotypeSizes[genotype][ecotype] * locusGenotypeBreedingValues[genotype];
            locusEcotypeAdditiveVariances[ecotype] += locusGenotypeEcotypeSizes[genotype][ecotype] * sqr(locusGenotypeBreedingValues[genotype]);
        }
        locusEcotypeMeanBreedingValues[ecotype] /= ecotypeSizes[ecotype];
        sumsq2var(locusEcotypeAdditiveVariances[ecotype], ecotypeSizes[ecotype], locusEcotypeMeanBreedingValues[ecotype]);
    }
}

void LocusVariables::completeLocusInteractionVariance()
{
    sumsq2var(locusInteractionVariance, popSize, 0.0);
}

void LocusVariables::completeLocusNonAdditiveVariances()
{
    sumsq2var(locusNonAdditiveVariance, popSize, 0.0);
    for (size_t ecotype = 0u; ecotype < 2u; ++ecotype) {
        locusEcotypeMeanNonAdditiveDeviations[ecotype] /= ecotypeSizes[ecotype];
        sumsq2var(locusEcotypeNonAdditiveVariances[ecotype], ecotypeSizes[ecotype], locusEcotypeMeanNonAdditiveDeviations[ecotype]);
    }
}

void LocusVariables::calcLocusEcotypeDifferentiations()
{
    // Calculate heterogeneity statistics
    locusFst = 1.0 - locusObservedHeterozygosity / locusExpectedHeterozygosity;
    locusPst = Xst(locusPhenotypicVariance, locusEcotypePhenotypicVariances, popSize, ecotypeSizes);
    locusGst = Xst(locusGeneticVariance, locusEcotypeGeneticVariances, popSize, ecotypeSizes);
    locusQst = Xst(locusAdditiveVariance, locusEcotypeAdditiveVariances, popSize, ecotypeSizes);
    locusCst = Xst(locusNonAdditiveVariance, locusEcotypeNonAdditiveVariances, popSize, ecotypeSizes);
}

void LocusVariables::accumulateLocusGeneticMoments()
{
    size_t genotype;
    size_t ecotype;
    double geneticValue;

    // Loop through individuals
    for (PInd pInd : individuals) {

        // Get individual information
        ecotype = pInd->getEcotype();
        genotype = pInd->getLocus(locus).alleleCount;
        geneticValue = pInd->getLocus(locus).locusGeneticValue;

        accumulateCensus(genotype, ecotype);
        accumulateAlleleCounts(genotype, ecotype);
        accumulateGeneticValues(genotype, ecotype, geneticValue);
        accumulateGeneticValuesByAlleleCounts(genotype, geneticValue);

    }
}

void LocusVariables::accumulateLocusIndividualResiduals()
{
    size_t ecotype;
    size_t genotype;
    double geneticValue;

    // Deviation due to epistasis, we need to loop again now that the regression is done
    for (PInd pInd : individuals) {

        // Get individual information
        ecotype = pInd->getEcotype();
        genotype = pInd->getLocus(locus).alleleCount;
        geneticValue = pInd->getLocus(locus).locusGeneticValue;

        // Measure squared deviations
        locusInteractionVariance += sqr(geneticValue - locusGenotypeMeanGeneticValues[genotype]);
        locusNonAdditiveVariance += sqr(geneticValue - locusGenotypeAdditiveExpectations[genotype]);
        locusEcotypeNonAdditiveVariances[ecotype] += sqr(geneticValue - locusGenotypeAdditiveExpectations[genotype]);
        locusEcotypeMeanNonAdditiveDeviations[ecotype] += geneticValue - locusGenotypeAdditiveExpectations[genotype];

    }
}











// Low-level functions

void Population::setResourceCapacities(const double &maxResourceCapacity, const double &habitatAsymmetry)
{
    resourceCapacities[0u].first = maxResourceCapacity;
    resourceCapacities[0u].second = (1.0 - habitatAsymmetry) * maxResourceCapacity;
    resourceCapacities[1u].first = (1.0 - habitatAsymmetry) * maxResourceCapacity;
    resourceCapacities[1u].first = maxResourceCapacity;
}

void Population::setReplenishRates(const double &maxResourceGrowth)
{
    replenishRates[0u].first = maxResourceGrowth;
    replenishRates[0u].second = maxResourceGrowth;
    replenishRates[1u].first = maxResourceGrowth;
    replenishRates[1u].second = maxResourceGrowth;
}

void Population::setResourceConsumption(const size_t &habitat)
{
    // Initialize consumption
    resourceConsumption[habitat].first = resourceConsumption[habitat].second = 0.0;

    // Loop through individuals
    for (auto itInd = idHabitatBoundaries[habitat].first;; ++itInd) {

        // The individual should be of the right habitat
        assert((*itInd)->getHabitat() == habitat);

        // Record attack rates
        std::pair<double, double> attackRates = (*itInd)->getAttackRates();

        // Are we in the burn-in period?
        if (nAccessibleResources < 2u) {
            attackRates.second = 0.0;
        }

        // Accumulate consumption
        resourceConsumption[habitat].first += attackRates.first;
        resourceConsumption[habitat].second += attackRates.second;

        if (itInd == idHabitatBoundaries[habitat].second)
        {
            break;
        }
    }
}

void Population::setResourceEquilibrium(const size_t &habitat)
{
    resourceEql[habitat].first = calcLogisticResourceEq(resourceCapacities[habitat].first,
                                                    replenishRates[habitat].first,
                                                    resourceConsumption[habitat].first);
    resourceEql[habitat].second = calcLogisticResourceEq(resourceCapacities[habitat].second,
                                                     replenishRates[habitat].second,
                                                     resourceConsumption[habitat].second);
}

void Population::assignFitnesses(const size_t &habitat, const double &ecoSelCoeff)
{
    for (PInd pInd : individuals) {

        // Reinforce stabilizing selection during burn-in period
        if (nAccessibleResources < 2u) {
            pInd->setBurninFitness(resourceEql[habitat], ecoSelCoeff);
        }

        // Or apply regular fitness function
        else {
            pInd->setFitness(resourceEql[habitat]);
        }
    }
}

void Population::classifyGenders(const bool &isFemaleHeterogamy)
{
    for (PInd pInd : individuals)
    {
        if (pInd->isFemale(isFemaleHeterogamy))
        {
            females.push_back(pInd);
        }
        else {
            males.push_back(pInd);
        }
    }
}

void Population::setMaleFitnesses(const size_t &nMales, const double &tiny)
{
    double sumMaleSuccesses = 0.0;
    for (size_t i = 0u; i < nMales; ++i) {

        maleSuccesses[i] = males[i]->getFitness();

        // Accumulate mating successes
        sumMaleSuccesses += maleSuccesses[i];

    }

    // All males have equal success if successes are too small
    if (sumMaleSuccesses < tiny) {
        maleSuccesses = std::vector<double>(nMales, 1.0);
    }
}

void Population::birth(const PInd &female, const ParameterSet &parameters, const GeneticArchitecture &geneticArchitecture)
{
    for (size_t idFather : female->getMates()) {
        offspring.push_back(new Individual(female, males[idFather], parameters, geneticArchitecture));
    }
}

void Population::emptyPopulation()
{
    individuals.erase(individuals.begin(), individuals.end());
}
