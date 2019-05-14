//
// Created by p278834 on 7-5-2019.
//

#include "Population.h"
#include "random.h"
#include <cassert>

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
    classifyGenders(parameters.isFemaleHeteroGamety);

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


// Accessory functions

double calcLogisticResourceEq(const double &resourceCapacity, const double &replenishRate, const double &consumption)
{
    double resource = resourceCapacity * (1.0 - consumption / replenishRate);
    return resource;
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
        assert((*ind)->getHabitat() == habitat);

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

void Population::classifyGenders(const bool &isFemaleHeteroGamety)
{
    for (PInd pInd : individuals)
    {
        if (pInd->isFemale(isFemaleHeteroGamety))
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



// Function to classify ecotypes
// Loop through the population and record attack rates
// Within local habitat
// Sort the attack rates along trade-off line with the function trade-off compare
// Loop through attack rates until the alternative resource becomes more advantageous
// That is the break-even point that delimits ecotypes