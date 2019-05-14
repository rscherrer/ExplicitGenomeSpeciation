//
// Created by p278834 on 7-5-2019.
//

#include "Population.h"
#include "random.h"
#include "queue"
#include <cassert>

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

void Population::dispersal(const ParameterSet& parameters)
{

    // Loop through the whole population if dispersal is high, all individuals have high chances to migrate
    if (parameters.dispersalRate > 0.5) {
        for (PInd pInd : individuals) {
            if (rnd::bernoulli(parameters.dispersalRate)) {
                pInd -> disperse(parameters.nHabitat);
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
        auto it = individuals.cbegin();
        size_t j = 0u;
        for (size_t i : migrants) {
            std::advance(it, i - j);
            (*it) -> disperse(parameters.nHabitat);
            j = i;
        }
    }
}

// This function needs massive refactoring
// In the new version of the model we use type II resource utilization and type II mate choice
// Which are much simpler than type I
// Of course we keep the scripts for implementing the type I in another branch of this project
// So they are not lost!!
// One difference though:
// If we remove T2 RU then we remove ecotypes, but ecotypes are needed for calculating ecological differentiation stats
// So we have to keep them
// So we have to sort individuals along the trade off line

// Maybe start by removing all instances of type I RU


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
    idFirstAndLast[0u].first = individuals.begin();
    idFirstAndLast[0u].second = iti.operator--();
    idFirstAndLast[1u].first = iti;
    idFirstAndLast[1u].second = individuals.end();
}

void Population::setResourceConsumption(const size_t &habitat)
{
    // Initialize consumption
    resourceConsumption[habitat].first = resourceConsumption[habitat].second = 0.0;

    // Loop through individuals
    for (auto ind = idFirstAndLast[habitat].first;; ++ind) {

        // The individual should be of the right habitat
        assert((*ind)->getHabitat() == habitat);

        // Record attack rates
        TradeOffPt attackRates = (*ind)->getAttackRates();

        // Are we in the burn-in period?
        if (nAccessibleResource < 2u) {
            attackRates.second = 0.0;
        }

        // Accumulate consumption
        resourceConsumption[habitat].first += attackRates.first;
        resourceConsumption[habitat].second += attackRates.second;

        if (ind == idFirstAndLast[habitat].second)
        {
            break;
        }
    }
}

double logisticResourceEq(const double &resourceCapacity, const double &replenishRate, const double &consumption)
{
    double resource = resourceCapacity * (1.0 - consumption / replenishRate);
    return resource;
}

void Population::setResourceEquilibrium(const size_t &habitat)
{
    resourceEql[habitat].first = logisticResourceEq(resourceCapacities[habitat].first,
            replenishRates[habitat].first, resourceConsumption[habitat].first);
    resourceEql[habitat].second = logisticResourceEq(resourceCapacities[habitat].second,
            replenishRates[habitat].second, resourceConsumption[habitat].second);
}

void Population::assignFitnesses(const size_t &habitat, const double &ecoSelCoeff)
{
    for (auto iti : individuals) {

        // Reinforce stabilizing selection during burn-in period
        if (nAccessibleResource < 2u) {
            iti->setBurninFitness(resourceEql[habitat], ecoSelCoeff);
        }

        // Or apply regular fitness function
        else {
            iti->setFitness(resourceEql[habitat]);
        }
    }
}

// Function to classify ecotypes
// Loop through the population and record attack rates
// Within local habitat
// Sort the attack rates along trade-off line with the function trade-off compare
// Loop through attack rates until the alternative resource becomes more advantageous
// That is the break-even point that delimits ecotypes

void Population::resourceDynamics(const size_t &habitat, const double &ecoSelCoeff)
{

    // Calculate total resource consumption in the current habitat
    setResourceConsumption(habitat);

    // Calculate equilibrium resource concentrations
    setResourceEquilibrium(habitat);

    // Assign individual fitnesses
    assignFitnesses(habitat, ecoSelCoeff);

}

void Population::classifyGenders(const bool &isFemaleHeteroGamety)
{
    for (iti : individuals)
    {
        if ((*iti)->isFemale(isFemaleHeteroGamety))
        {
            females.push(*iti);
        }
        else {
            males.push_back(*iti);
        }
    }
}

void Population::emptyPopulation()
{
    individuals.erase(individuals.begin(), individuals.end());
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

void Population::birth(const PInd &female, const ParameterSet &parameters, const GeneticArchitecture &geneticArchitecture)
{
    for (size_t idFather : female->getMates()) {
        offspring.push_back(new Individual(female, males[idFather], parameters, geneticArchitecture));
    }
}


void Population::survival(const double &survivalProb)
{
    emptyPopulation();
    for (PInd ind : males) {
        if (ind->survive(survivalProb)) {
            individuals.push_back(ind);
        }
    }
    for (PInd ind : females) {
        if (ind->survive(survivalProb)) {
            individuals.push_back(ind);
        }
    }
    for (PInd ind : offspring) {
        individuals.push_back(ind);
    }
}