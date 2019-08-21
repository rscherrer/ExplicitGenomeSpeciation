#include "Population.h"
#include "Individual.h"
#include "utils.h"
#include "Random.h"
#include <cassert>
#include <iostream>


typedef std::discrete_distribution<size_t> Discrete;


/// Constructor
Population::Population(const size_t &popsize,
 const Genome &genome, const MultiNet &networks) :
    individuals(populate(popsize, genome, networks)),
    females({ }),
    males({ }),
    offspring({ }),
    survivors({ }),
    capacity({100.0, 100.0}),
    replenish({1.0, 1.0}),
    resources(capacity)
{
}


/// Function to initialize a population of individuals
Crowd Population::populate(const size_t &popsize,
 const Genome &genome, const MultiNet &networks)
{

    Crowd indivs;

    for (size_t ind = 0u; ind < popsize; ++ind)
        indivs.push_back(new Individual(genome, networks));

    return indivs;

}


Crowd Population::emigrate(const double &rate)
{
    // Loop through individuals in the pop
    // Everyone has a change to migrate

    Crowd migrants;
    Crowd residents;

    for (auto ind : individuals)
        if (rnd::bernoulli(rate))
            migrants.push_back(ind);
        else
            residents.push_back(ind);

    assert(residents.size() + migrants.size() == individuals.size());
    assert(migrants.size() <= individuals.size());
    assert(residents.size() <= individuals.size());

    individuals.clear();
    individuals = residents;

    return migrants;

}


/// Resource consumption
void Population::consume()
{

    // cap is the carrying capacity (K)
    // repl is the rate of replenishment of the resource (r)
    // The resource dynamics are as follows:
    // dR/dt = r (1 - R / K) - C R, where C is the consumption

    // Calculate the total amount of food consumed
    dVector consumed{0.0, 0.0};
    for (auto ind : individuals) {
        dVector rates = ind->getFeedingRates();
        for (size_t res = 0u; res < 2u; ++res)
            consumed[res] += rates[res];
    }

    assert(consumed[0u] >= 0.0);
    assert(consumed[1u] >= 0.0);

    // Update the resource equilibrium state
    for (size_t res = 0u; res < 2u; ++res)
        resources[res] = capacity[res] * (1.0 - consumed[res] / replenish[res]);

    assert(resources[0u] >= 0.0);
    assert(resources[1u] >= 0.0);

    // Split the resource among the individuals
    for (auto ind : individuals)
        ind->feed(resources);

}


/// Asexual reproduction function
void Population::reproduceAsexual(const double &birth,
 const Genome &genome, const MultiNet &networks)
{
    // Everybody gets a chance to produce babies
    size_t nAdults = individuals.size();
    while (nAdults) {
        size_t nOffspring = rnd::poisson(birth);
        while (nOffspring) {
            offspring.push_back(new Individual(genome, networks));
            --nOffspring;
        }
        --nAdults;
    }
}


/// Sexual reproduction function
void Population::reproduce(const double &birth, const double &strength,
 const Genome &genome, const MultiNet &networks, const double &cost)
{

    // Sort out moms and dads
    for (auto ind : individuals)
        if (ind->getGender())
            females.push_back(ind);
        else
            males.push_back(ind);

    // Prepare a weighted lottery based on male mating successes
    dVector successes;
    for (auto male : males)
        successes.push_back(male->getFitness());

    assert(successes.size() == males.size());

    Discrete maleMarket(successes.begin(), successes.end());

    // Sample the duration of the mating season this year
    const size_t seasonEnd = rnd::geometric(cost);

    // Every mom gets a chance to produce babies
    for (auto mom : females) {

        size_t nOffspring = rnd::poisson(birth * mom->getFitness());

        Haplotype egg = mom->recombine(genome.locations, genome.chromosomes);
        mom->mutate(egg);

        size_t time = 0u;

        // The mating season begins...
        while (nOffspring && time < seasonEnd) {

            // Sample a male
            const size_t encounter = maleMarket(rnd::rng);
            assert(encounter < males.size());
            auto dad = males[encounter];

            Haplotype sperm = dad->recombine(genome.locations,
             genome.chromosomes);
            dad->mutate(sperm);

            if (mom->acceptMate(dad->getEcoTrait(), strength)) {
                offspring.push_back(new Individual(genome, networks, egg,
                 sperm));
                --nOffspring;
            }

            ++time;
        }

    }

    females.clear();
    males.clear();
}


/// Function to make it to the next generation
bool Population::survive(const double &survival)
{

    // Sample life or death for every adult
    for (auto ind : individuals)
        if (rnd::bernoulli(survival))
            survivors.push_back(ind);

    individuals.clear();
    assert(individuals.size() == 0u);

    const size_t nSurvivors = survivors.size();
    const bool isAlive = nSurvivors != 0u;

    // Survivors make it to the next generation
    if (isAlive) {
        for (auto ind : survivors)
            individuals.push_back(ind);
        survivors.clear();
    }

    // Offspring make it to the next generation
    const size_t nOffspring = offspring.size();
    if (nOffspring != 0u) {
        for (auto ind : offspring)
            individuals.push_back(ind);
        offspring.clear();
    }

    assert(survivors.size() == 0u);
    assert(offspring.size() == 0u);
    assert(individuals.size() == nSurvivors + nOffspring);

    return isAlive;
}

