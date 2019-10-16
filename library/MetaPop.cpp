#include "MetaPop.h"


// Initialization

Crowd MetaPop::populate(const Param &p, const GenArch &arch)
{

    // Generate a pool of individuals
    const size_t n = utl::sumu(p.demesizes);
    size_t n0 = p.demesizes[0u];
    assert(n0 <= n);
    Crowd indivs;
    indivs.reserve(n); // seems to be causing trouble

    for (size_t ind = 0u; ind < n; ++ind) {
        indivs.push_back(Individual(p, arch));
        if (ind >= n0) indivs.back().disperse();
    }

    assert(indivs.size() == n);

    return indivs;

}


// Life cycle of the population

void MetaPop::cycle(const Param &p, const GenArch &arch)
{

    // Dispersal
    if (!isburnin) disperse(p);

    // Consumption
    consume(p);

    // Reproduction
    reproduce(p, arch);

    // Survival
    survive(p);

}

void MetaPop::disperse(const Param &p)
{
    // Sample migrants across the population
    // Change the habitat attribute of these migrants

    vecDbl probs = utl::ones(population.size());
    size_t nmigrants = rnd::binomial(population.size(), p.dispersal);
    while (nmigrants) {
        const size_t mig = rnd::sample(probs);
        probs[mig] = 0.0;
        population[mig].disperse();
        --nmigrants;
    }
}

void MetaPop::consume(const Param &p)
{

    // Calculate the sum of feeding efficiencies in each habitat
    Matrix sumfeed = utl::zeros(2u, 2u);
    for (size_t i = 0u; i < population.size(); ++i) {
        const size_t hab = population[i].getHabitat();
        sumfeed[hab][0u] += population[i].getFeeding(0u);

        // Feed only on resource 0 during burnin
        if (!isburnin) sumfeed[hab][1u] += population[i].getFeeding(1u);
    }

    // Convert sums of feeding efficiencies into relative consumed food (C)
    Matrix consumed = utl::zeros(2u, 2u);
    for (size_t hab = 0u; hab < 2u; ++hab) {
        for (size_t res = 0u; res < 2u; ++res) {
            consumed[hab][res] = 1.0;
            consumed[hab][res] -= exp(-p.maxfeed * sumfeed[hab][res]);
            assert(consumed[hab][res] >= 0.0);
            assert(consumed[hab][res] <= 1.0);
        }
    }

    // Resource capacity without consumption (K)
    resources = utl::zeros(2u, 2u);
    resources[0u][0u] = p.capacity;
    resources[0u][1u] = isburnin ? 0.0 : p.capacity * p.hsymmetry;
    resources[1u][0u] = p.capacity * p.hsymmetry;
    resources[1u][1u] = isburnin ? 0.0 : p.capacity;

    assert(resources[0u][0u] >= 0.0);
    assert(resources[0u][1u] >= 0.0);
    assert(resources[1u][0u] >= 0.0);
    assert(resources[1u][1u] >= 0.0);

    // Absolute amount of food consumed (C K / r)
    for (size_t hab = 0u; hab < 2u; ++hab) {
        for (size_t res = 0u; res < 2u; ++res) {
            consumed[hab][res] *= resources[hab][res] / p.replenish;
            assert(consumed[hab][res] >= 0.0);
            assert(consumed[hab][res] <= resources[hab][res]);
        }
    }

    // Assign individual fitness and ecotypes
    for (size_t i = 0u; i < population.size(); ++i) {
        const size_t hab = population[i].getHabitat();
        double fitness = 0.0;
        vecDbl food = utl::zeros(2u);
        for (size_t res = 0u; res < 2u; ++res) {
            const double feed = population[i].getFeeding(res);
            if (sumfeed[hab][res])
                food[res] = consumed[hab][res] * feed / sumfeed[hab][res];
            fitness += food[res];
        }
        size_t ecotype = food[1u] > food[0u];
        assert(fitness >= 0.0);
        assert(ecotype == 0u || ecotype == 1u);
        population[i].feed(fitness, ecotype);
    }

    // Update the resource equilibrium (K - C K / r)
    for (size_t hab = 0u; hab < 2u; ++hab) {
        for (size_t res = 0u; res < 2u; ++res) {
            resources[hab][res] -= consumed[hab][res];
        }
    }
}

void MetaPop::reproduce(const Param &p, const GenArch &arch)
{

    // Table to count the sexes in both habitats
    sexcounts = utl::uzeros(2u, 2u);

    // Males' fitness determine their encounter probabilities
    Matrix probs = utl::zeros(2u, population.size());
    for (size_t i = 0u; i < population.size(); ++i) {

        const size_t sex = population[i].getGender();
        const size_t hab = population[i].getHabitat();
        ++sexcounts[hab][sex];

        if (!sex) {

            probs[hab][i] = population[i].getFitness();

            // Modify mating success during burnin
            if (isburnin) {
                const double y = population[i].getMatePref();
                probs[hab][i] *= exp(-p.ecosel * utl::sqr(y));
            }
        }
    }

    if (!(sexcounts[0u][0u] + sexcounts[1u][0u])) return; // exit if no males
    if (!(sexcounts[0u][1u] + sexcounts[1u][1u])) return; // exit if no females

    // Determine the length of the mating season
    const size_t seasonend = rnd::geometric(p.matingcost);

    if (!seasonend) return;

    const size_t nparents = population.size();

    // For each parent...
    for (size_t mom = 0u; mom < nparents; ++mom) {

        const size_t sex = population[mom].getGender();
        const size_t hab = population[mom].getHabitat();

        // Is it a female and are there males around?
        if (sex && sexcounts[hab][0u]) {

            // Progress throughout the mating season
            size_t timeleft = seasonend;
            while (timeleft) {

                --timeleft;

                // And encounters males one at a time with replacement
                const size_t dad = rnd::sample(probs[hab]);
                const double maletrait = population[dad].getEcoTrait();

                // If the female accepts to mate                
                if (rnd::bernoulli(population[mom].mate(maletrait, p))) {

                    // Determine fecundity
                    double fecundity = p.birth * population[mom].getFitness();

                    // Modify fecundity during burnin
                    if (isburnin) {
                        const size_t y = population[mom].getMatePref();
                        fecundity *= exp(-p.ecosel * utl::sqr(y));
                    }

                    // Sample clutch size
                    size_t noffspring = rnd::poisson(fecundity);
                    while (noffspring) {

                        // Give birth
                        population.push_back(Individual(p, arch,
                         population[mom], population[dad]));
                        --noffspring;
                    }

                    // End the mating season if female has mated
                    break;
                }
            }
        }
    }
}

// Lambda for removing dead individuals
auto burry = [](Individual ind) -> bool
{
    return !ind.isalive();
};

void MetaPop::survive(const Param &p)
{
    // Sample survival for each individual
    size_t nsurv = 0u;
    for (size_t i = 0u; i < population.size(); ++i) {
        population[i].survive(rnd::bernoulli(p.survival));
        if (population[i].isalive()) ++nsurv;
    }

    // Remove dead individuals
    auto it = std::remove_if(population.begin(), population.end(), burry);
    population.erase(it, population.end());
    population.shrink_to_fit();
    assert(population.size() == nsurv);

}


// Others

bool MetaPop::isextinct() const
{
    return population.size() == 0u;
}

void MetaPop::exitburnin()
{
    isburnin = false;
}
