#include "Population.h"
#include "Individual.h"
#include "utils.h"
#include "Random.h"
#include <cassert>
#include <iostream>


typedef std::discrete_distribution<size_t> Discrete;

/// Function to initialize a population of individuals
Crowd Population::populate(const size_t &n, const double &ecosel,
 const double &maxfeeding, const GeneticArchitecture &arch)
{

    Crowd indivs;

    for (size_t ind = 0u; ind < n; ++ind) {
        auto indiv = new Individual(arch, ecosel, maxfeeding, arch.snpFreq);
        indivs.push_back(indiv);
    }

    return indivs;

}


Crowd Population::emigrate(const double &rate)
{
    // Loop through individuals in the pop
    // Everyone has a change to migrate

    Crowd migrants;
    Crowd residents;

    size_t nmigrants = rnd::binomial(individuals.size(), rate);
    vecDbl probs = ones(nmigrants);
    Haplotype whom = falses(individuals.size());

    while (nmigrants) {
        const size_t mig = rnd::sample(probs);
        whom[mig] = true;
        probs[mig] = 0.0;
        --nmigrants;
    }

    for (size_t ind = 0u; ind < individuals.size(); ++ind) {
        if(whom[ind])
            migrants.push_back(individuals[ind]);
        else
            residents.push_back(individuals[ind]);
    }

    assert(residents.size() + migrants.size() == individuals.size());
    assert(migrants.size() <= individuals.size());
    assert(residents.size() <= individuals.size());

    individuals.clear();
    individuals = residents;

    return migrants;

}


void Population::immigrate(const Crowd &newcomers)
{
    for (auto ind : newcomers)
        individuals.push_back(ind);
}


/// Resource consumption
void Population::consume()
{

    // cap is the carrying capacity (K)
    // repl is the rate of replenishment of the resource (r)
    // The resource dynamics are as follows:
    // dR/dt = r (1 - R / K) - C R, where C is the consumption

    // Calculate the total amount of food consumed
    vecDbl consumed{0.0, 0.0};
    for (auto ind : individuals) {
        vecDbl rates = ind->getFeedingRates();
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


void Population::burninConsume()
{
    // Calculate the total amount of food consumed
    double consumed = 0.0;
    for (auto ind : individuals) {
        double rate = ind->getFeedingRates()[0u];
        consumed += rate;
    }

    assert(consumed >= 0.0);

    // Update the resource equilibrium state
    resources[0u] = capacity[0u] * (1.0 - consumed / replenish[0u]);
    resources[1u] = 0.0;

    assert(resources[0u] >= 0.0);

    // Split the resource among the individuals
    for (auto ind : individuals)
        ind->feed(resources);
}

void Population::sortSexes()
{

    females.clear();
    males.clear();

    // Sort out moms and dads
    for (auto ind : individuals)
        if (ind->getGender())
            females.push_back(ind);
        else
            males.push_back(ind);
}

/// Sexual reproduction function
void Population::reproduce(const double &birth, const double &sexsel,
 const double &cost, const double &ecosel, const double &maxfeeding,
  const GeneticArchitecture &arch)
{

    if (!(females.size() > 0u) || !(males.size() > 0u)) return;

    // Prepare a weighted lottery based on male mating successes
    vecDbl successes;
    for (auto male : males)
        successes.push_back(male->getFitness());

    assert(successes.size() == males.size());

    Discrete maleMarket(successes.begin(), successes.end());

    // Sample the duration of the mating season this year
    const size_t seasonEnd = rnd::geometric(cost);

    // Every mom gets a chance to produce babies
    for (auto mom : females) {

        size_t nOffspring = rnd::poisson(birth * mom->getFitness());

        Haplotype egg = mom->recombine(arch);
        mom->mutate(egg);

        size_t time = 0u;

        // The mating season begins...
        while (nOffspring && time < seasonEnd) {

            // Sample a male
            const size_t encounter = maleMarket(rnd::rng);
            assert(encounter < males.size());
            auto dad = males[encounter];

            Haplotype sperm = dad->recombine(arch);
            dad->mutate(sperm);

            if (mom->acceptMate(dad->getEcoTrait(), sexsel)) {
                auto off = new Individual(arch, egg, sperm, ecosel, maxfeeding);
                offspring.push_back(off);
                --nOffspring;
            }

            ++time;
        }
    }
}


void Population::burninReproduce(const double &birth, const double &sexsel,
 const double &cost, const double &ecosel, const double &maxfeeding,
  const GeneticArchitecture &arch)
{
    if (!(females.size() > 0u) || !(males.size() > 0u)) return;

    // Prepare a weighted lottery based on male mating successes
    vecDbl successes;
    for (auto male : males) {
        double burninFactor = exp(- ecosel * sqr(male->getMatePref()));
        double success = male->getFitness() * burninFactor;
        successes.push_back(success);
    }

    assert(successes.size() == males.size());

    Discrete maleMarket(successes.begin(), successes.end());

    // Sample the duration of the mating season this year
    const size_t seasonEnd = rnd::geometric(cost);

    // Every mom gets a chance to produce babies
    for (auto mom : females) {

        double burninFactor = - ecosel * sqr(mom->getMatePref());
        double success = mom->getFitness() * burninFactor;
        size_t nOffspring = rnd::poisson(birth * success);

        Haplotype egg = mom->recombine(arch);
        mom->mutate(egg);

        size_t time = 0u;

        // The mating season begins...
        while (nOffspring && time < seasonEnd) {

            // Sample a male
            const size_t encounter = maleMarket(rnd::rng);
            assert(encounter < males.size());
            auto dad = males[encounter];

            Haplotype sperm = dad->recombine(arch);
            dad->mutate(sperm);

            if (mom->acceptMate(dad->getEcoTrait(), sexsel)) {
                auto off = new Individual(arch, egg, sperm, ecosel, maxfeeding);
                offspring.push_back(off);
                --nOffspring;
            }

            ++time;
        }
    }
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

void Population::resetEcoTraits(const double &value, const double &sel,
 const double &max)
{
    for (auto ind : individuals)
        ind->setEcoTrait(value, sel, max);
}

void Population::resetMatePrefs(const double &value)
{
    for (auto ind : individuals)
        ind->setMatePref(value);
}

void Population::resetGenders(const bool &sex)
{
    for (auto ind : individuals)
        ind->setGender(sex);
}

void Population::resetEcotypes(const size_t &ecotype)
{
    for (auto ind : individuals)
        ind->setEcotype(ecotype);
}
