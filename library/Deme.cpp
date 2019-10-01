#include "Deme.h"

typedef std::discrete_distribution<size_t> Discrete;


// Initialization

Crowd Deme::populate(const size_t &n, const double &ecosel,
 const double &maxfeeding, const GenArch &arch)
{

    // Crowd indivs;

    Crowd indivs;

    for (size_t ind = 0u; ind < n; ++ind) {
        auto indiv = new Individual(arch, ecosel, maxfeeding, arch.snpFreq);
        indivs.push_back(indiv);
    }

    return indivs;

}


// Life cycle

void Deme::exitBurnIn()
{
    burnin = false;
}

void Deme::sortSexes()
{

    females.clear();
    males.clear();
    females.shrink_to_fit();
    males.shrink_to_fit();

    // Sort out moms and dads
    for (auto ind : individuals)
        if (ind->getGender())
            females.push_back(ind);
        else
            males.push_back(ind);
}

Crowd Deme::emigrate(const double &rate)
{
    // Loop through individuals in the pop
    // Everyone has a change to migrate

    Crowd migrants;
    Crowd residents;

    size_t nmigrants = rnd::binomial(individuals.size(), rate);
    vecDbl probs = utl::ones(individuals.size());

    // All individuals that moved
    while (nmigrants) {
        const size_t mig = rnd::sample(probs);
        migrants.push_back(individuals[mig]);
        probs[mig] = 0.0;
        --nmigrants;
    }

    // All individuals that stayed
    for (size_t ind = 0u; ind < individuals.size(); ++ind)
        if (probs[ind] == 1.0)
            residents.push_back(individuals[ind]);

    assert(residents.size() + migrants.size() == individuals.size());
    assert(migrants.size() <= individuals.size());
    assert(residents.size() <= individuals.size());

    individuals.clear();
    individuals.shrink_to_fit();
    individuals = residents;

    return migrants;

}

void Deme::immigrate(const Crowd &newcomers)
{
    for (auto ind : newcomers)
        individuals.push_back(ind);
}

void Deme::consume()
{

    // cap is the carrying capacity (K)
    // repl is the rate of replenishment of the resource (r)
    // The resource dynamics are as follows:
    // dR/dt = r (1 - R / K) - C R, where C is the consumption

    // Calculate the total amount of food consumed
    vecDbl consumed{0.0, 0.0};
    for (auto ind : individuals) {
        vecDbl rates = ind->getFeedingRates();
        if (burnin) rates[1u] = 0.0;
        for (size_t res = 0u; res < 2u; ++res)
            consumed[res] += rates[res];
    }

    assert(consumed[0u] >= 0.0);
    assert(consumed[1u] >= 0.0);

    // Update the resource equilibrium state
    for (size_t res = 0u; res < 2u; ++res)
        resources[res] = capacity[res] * (1.0 - consumed[res] / replenish[res]);

    if (burnin) resources[1u] = 0.0;

    assert(resources[0u] >= 0.0);
    assert(resources[1u] >= 0.0);

    // Split the resource among the individuals
    for (auto ind : individuals) {
        ind->feed(resources);
    }

}

void Deme::reproduce(const double &birth, const double &sexsel,
 const double &cost, const double &ecosel, const double &maxfeeding,
  const GenArch &arch)
{

    if (!(females.size() > 0u) || !(males.size() > 0u)) return;

    // Prepare a weighted lottery based on male mating successes
    vecDbl successes;
    for (auto male : males) {
        double fit = male->getFitness();
        if (burnin) fit *= exp(- ecosel * utl::sqr(male->getMatePref()));
        successes.push_back(fit);
    }    

    assert(successes.size() == males.size());

    Discrete maleMarket(successes.begin(), successes.end());

    // Sample the duration of the mating season this year
    const size_t seasonEnd = rnd::geometric(cost);

    // Every mom gets a chance to produce babies
    for (auto mom : females) {

        // Produce a certain number of offspring
        double fecundity = birth * mom->getFitness();
        if (burnin) fecundity *= exp(- ecosel * utl::sqr(mom->getMatePref()));
        size_t nOffspring = rnd::poisson(fecundity);

        size_t time = 0u;

        // Should a female have offspring from different fathers?
        // Or from the same one?
        // Option 2 would increase the intensity of sexual selection
        // Option 1 would be closer to species doing extra pair mating?
        // Can be discussed during lab meeting

        // The mating season begins...
        while (time < seasonEnd) {

            // Sample a male
            const size_t encounter = maleMarket(rnd::rng);
            assert(encounter < males.size());
            auto dad = males[encounter];

            // Upon mating...
            if (mom->acceptMate(dad->getEcoTrait(), sexsel)) {

                while (nOffspring) {

                    // Produce gametes
                    Haplotype egg = mom->recombine(arch);
                    mom->mutate(egg);

                    Haplotype sperm = dad->recombine(arch);
                    dad->mutate(sperm);

                    // And a baby
                    auto off = new Individual(arch, egg, sperm, ecosel,
                     maxfeeding);
                    offspring.push_back(off);
                    --nOffspring;

                }

                break; // end the mating season
            }

            ++time;
        }
    }
}

bool Deme::survive(const double &survival)
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


// Getters

size_t Deme::getEcotype(const size_t &i) const
{
    return individuals[i]->getEcotype();
}

double Deme::getGenValue(const size_t &i, const size_t &t) const
{
    return individuals[i]->getGenValue(t);
}

double Deme::getTraitValue(const size_t &i, const size_t &t) const
{
    return individuals[i]->getTraitValue(t);
}

double Deme::getLocusValue(const size_t &i, const size_t &l) const
{
    return individuals[i]->getLocusValue(l);
}

size_t Deme::getZygosity(const size_t &i, const size_t &l) const
{
    return individuals[i]->getZygosity(l);
}


// Resetters used in testing

void Deme::resetEcoTraits(const double &value, const double &sel,
 const double &max)
{
    for (auto ind : individuals)
        ind->setEcoTrait(value, sel, max);
}

void Deme::resetMatePrefs(const double &value)
{
    for (auto ind : individuals)
        ind->setMatePref(value);
}

void Deme::resetGenders(const bool &sex)
{
    for (size_t ind = 0u; ind < individuals.size(); ++ind) {
        individuals[ind]->setGender(sex);
    }
}

void Deme::resetEcotypes(const size_t &ecotype)
{
    for (size_t ind = 0u; ind < individuals.size(); ++ind)
        individuals[ind]->resetEcotype(ecotype);
}

