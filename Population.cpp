//
// Created by p278834 on 7-5-2019.
//

#include "Population.h"
#include "random.h"
#include "queue"
#include <cassert>

Population::Population(const ParameterSet &parameters, const GeneticArchitecture &geneticArchitecture)
{
    if (parameters.sequence.size() == parameters.nBits) {
        for (size_t i = 0u; i < parameters.nIndividualInit; ++i) {
            individuals.push_back(new Individual(parameters.sequence, parameters, geneticArchitecture));
        }
    }
    else {
        for (size_t i = 0u; i < parameters.nIndividualInit; ++i) {
            individuals.push_back(new Individual(parameters, geneticArchitecture));
        }
    }
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


void Population::resourceDynamics(const size_t habitat, const ParameterSet &parameters)
{

    // Initialize consumption
    resourceConsumption[habitat].first = resourceConsumption[habitat].second = 0.0;

    auto iti = individuals.begin();
    std::list<TradeOffPt> vecAttackRates;

    // Loop through individuals from the end
    for (auto itj = individuals.end(); iti != itj;) {

        // If the individual lives in the current habitat
        if ((*iti)->getHabitat() == habitat) {

            // Record attack rates
            TradeOffPt attackRates = (*iti)->getAttackRates();

            // Are we in the burnin period?
            if (nAccessibleResource < 2u) {
                attackRates.second = 0.0;
            }

            // Accumulate
            vecAttackRates.push_back(attackRates);

            // For the moment, assume that the individual utilises the first resource
            resourceConsumption[habitat].first += attackRates.first;

            // But sum the attack rates on the second resource anyway if type II resource utilisation
            resourceConsumption[habitat].second += attackRates.second;


        }

        // Move individuals from the wrong habitat to the end
        else {
            --itj;
            std::swap(*iti, *itj);
        }
    }
    individuals.erase(individuals.begin(), iti);

    // Determine equilibrium scaled resource densities and assign final ecotype to the individual

    // Sort individuals along the trade-off line
    pts.sort(tradeOffCompare);
    breakEvenPoint = pts.back();
    breakEvenPoint.first *= 0.5;
    breakEvenPoint.second = 0.0;

    // Initialize resource dynamics
    resourceEql[habitat].first = (habitat == 0u ? 1.0 : 1.0 - parameters.habitatAsymmetry) / (1.0 + parameters.alpha * resourceConsumption[hab].first);
    resourceEql[habitat].second = (habitat == 1u ? 1.0 : 1.0 - parameters.habitatAsymmetry) / (1.0 + parameters.alpha * resourceConsumption[hab].second);

    // Find resource equilibrium and break-even point by looping along the trade-off line
    // (used only for ecotype classification in type II resource utilisation)
    for (const TradeOffPt &pt : pts) {
        const bool isResource1MoreAdvantageous = pt.first * resourceEql[habitat].first > pt.second * resourceEql[habitat].second;
        if (isResource1MoreAdvantageous) {
            // Set break-even point to be used in later ecotype classification
            breakEvenPoint = pt;
            break;
        }
    }


}

void Population::reproduction(const size_t hab, const ParameterSet &parameters)
{
    // Initialize the state of affairs of the population
    std::queue<PInd> females;
    std::vector<PInd> males;

    /*
     * // Classify gender (this could move to another part)
            if ((*iti)->isFemale(parameters.isFemaleHeteroGamety)) {
                females.push(*iti);
            }
            else {
                males.push_back(*iti);
            }
            ++iti;
     * */

    // Mate choice and offspring production
    const size_t nFemales = genderCounts[hab].first = females.size();
    const size_t nMales = genderCounts[hab].second = males.size();

    // Terminate if there are no males or no females
    if (nFemales == 0u || nMales == 0u) {
        return;
    }

    // Compute reproductive success for males
    std::vector<double> maleSuccess(nMales);
    double sumMaleSuccess = 0.0;

    // Loop through males
    for(size_t i = 0u; i < nMales; ++i) {

        // Pick the resource that yields the highest payoff (not if type II resource utilisation)
        TradeOffPt pt = males[i]->getAttackRate();
        if (nAccessibleResource < 2u) pt.second = 0.0;  // If burnin period
        maleSuccess[i] = parameters.isTypeIIResourceUtilisation ? pt.first * resourceEql[hab].first + pt.second * resourceEql[hab].second : std::max(pt.first * resourceEql[hab].first, pt.second * resourceEql[hab].second);

        // Add stabilising selection on mating trait during burn-in period
        if (nAccessibleResource < 2u) {
            maleSuccess[i] *= males[i]->getBurnInRpSc(parameters.ecoSelCoeff);
        }

        // Accumulate mating successes
        sumMaleSuccess += maleSuccess[i];

    }

    // In case of equal mating successes
    if (sumMaleSuccess < parameters.tiny) {
        maleSuccess = std::vector<double>(nMales, 1.0);
    }

    // Male market
    std::discrete_distribution<size_t> maleMarket(maleSuccess.begin(), maleSuccess.end());

    // Set the length of the mating season
    const size_t seasonEnd = rnd::geometric(parameters.mateEvaluationCost);

    // Mate choice and reproduction
    while (!females.empty()) {

        PInd fem = females.front();
        females.pop();

        // Compute female reproductive success
        TradeOffPt pt = fem->getAttackRate();
        if (nAccessibleResource < 2u) {
            pt.second = 0.0;
        }
        double femaleSuccess = parameters.isTypeIIResourceUtilisation ? pt.first * resourceEql[hab].first + pt.second * resourceEql[hab].second : std::max(pt.first * resourceEql[hab].first, pt.second * resourceEql[hab].second);
        if (nAccessibleResource < 2u) {
            femaleSuccess *= fem->getBurnInRpSc(parameters.ecoSelCoeff);
        }

        // Sample family size for female
        size_t nOffspring = rnd::poisson(parameters.beta * femaleSuccess);

        // Mate choice for each clutch
        fem->prepareChoice();
        for (size_t t = 0u; nOffspring && t < seasonEnd; ++t) {

            // Sample a male
            const size_t j = maleMarket(rnd::rng);

            // If mating is successful
            if (fem->acceptMate(males[j], parameters)) {

                // Add offspring to the population
                individuals.push_back(new Individual(fem, males[j], parameters, genome));

                // Check if the offspring survives development
                if (parameters.costIncompat > 0.0) {
                    if (rnd::bernoulli(individuals.back()->getViability())) {
                        individuals.pop_back();
                    }
                }

                --nOffspring;
            }
        }

        // Female survival
        if (rnd::bernoulli(parameters.survivalProb)) {
            individuals.push_back(fem);
        }
        else {
            delete fem;
        }
    }

    // Male survival
    for (size_t i = 0u; i < nMales; ++i) {
        if (rnd::bernoulli(parameters.survivalProb)) {
            individuals.push_back(males[i]);
        }
        else {
            delete males[i];
        }
        males[i] = nullptr;
    }
}