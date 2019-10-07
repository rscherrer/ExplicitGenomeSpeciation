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
    // Calculate the total amount of food consumed in each habitat
    Matrix consumed = utl::zeros(2u, 2u);
    for (size_t i = 0u; i < population.size(); ++i) {
        const size_t hab = population[i].getHabitat();
        consumed[hab][0u] += population[i].getFeeding(0u);
        if (!isburnin) consumed[hab][1u] += population[i].getFeeding(1u);
    }

    assert(consumed[0u][0u] >= 0.0);
    assert(consumed[0u][1u] >= 0.0);
    assert(consumed[1u][0u] >= 0.0);
    assert(consumed[1u][1u] >= 0.0);

    // Update the resource equilibria
    for (size_t hab = 0u; hab < 2u; ++hab) {
        for (size_t res = 0u; res < 2u; ++res) {
            if (res == 1u && isburnin) {
                resources[hab][res] = 0.0;
            }
            else {
                resources[hab][res] = p.capacity;
                if (hab == res) resources[hab][res] *= p.hsymmetry;
                resources[hab][res] *= (1.0 - consumed[hab][res] / p.replenish);
                if (resources[hab][res] < 0.0) resources[hab][res] = 0.0;
            }
        }
    }

    // Could the resources be negative? If too much consumption?
    assert(resources[0u][0u] >= 0.0) ;
    assert(resources[0u][1u] >= 0.0);
    assert(resources[1u][0u] >= 0.0);
    assert(resources[1u][1u] >= 0.0);

    // Assign individual fitness and ecotypes
    for (size_t i = 0u; i < population.size(); ++i)
        population[i].feed(resources[population[i].getHabitat()]);
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

    // For each individual...
    for (size_t mom = 0u; mom < population.size(); ++mom) {

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
                if (population[mom].accept(maletrait, p)) {

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
auto burry = [&](Individual ind) -> bool
{
    return !ind.isalive();
};

void MetaPop::survive(const Param &p)
{
    // Sample survival for each individual
    for (size_t i = 0u; i < population.size(); ++i)
        population[i].survive(rnd::bernoulli(p.survival));

    // Remove dead individuals
    auto it = std::remove_if(population.begin(), population.end(), burry);
    population.erase(it, population.end());
    population.shrink_to_fit();
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


// For analysis

MatUns MetaPop::matingtrials(const Param &p) const
{
    // Count homogamic and heterogamic crossings

    // Table of crossings
    MatUns m = utl::uzeros(2u, 2u);

    const size_t nmal = sexcounts[0u][0u] + sexcounts[1u][0u];
    const size_t nfem = sexcounts[0u][1u] + sexcounts[1u][1u];

    if (nfem == 0u || nmal == 0u) return m;

    vecUns males;
    males.reserve(nmal);
    for (size_t i = 0u; i < population.size(); ++i)
        if (!population[i].getGender()) males.push_back(i);

    assert(males.size() == nmal);

    // For each female in the metapopulation...
    for (size_t fem = 0u; fem < population.size(); ++fem) {
        if (population[fem].getGender()) {

            const size_t ecofem = population[fem].getEcotype();

            // Determine the number of males encountered
            size_t nencounters = rnd::poisson(1.0 / p.matingcost);

            // For each male encountered...
            while (nencounters) {

                // Sample the male with replacement
                const size_t candidate = males[rnd::random(nmal)];
                const size_t maletrait = population[candidate].getEcoTrait();
                const size_t ecomal = population[candidate].getEcotype();

                // Evaluate the male and update the crossings accordingly
                if (population[fem].accept(maletrait, p)) {
                    ++m[ecofem][ecomal];
                }

                --nencounters;
            }
        }
    }

    return m;
}
