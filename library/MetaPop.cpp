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
        if (ind >= n0) indivs.back().disperse(); // is this not making a copy? test it
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
    if (p.dispersal < 0.5) {
        auto hasmigrated = boost::dynamic_bitset<>(population.size());
        size_t nmigrants = rnd::binomial(population.size(), p.dispersal);
        size_t t = 0u;
        while (nmigrants) {
            const size_t mig = rnd::random(population.size());
            if (!hasmigrated.test(mig)) {
                hasmigrated.set(mig);
                population[mig].disperse();
                --nmigrants;
                t = 0u;
            }

            // If we run too long without finding an individual to migrate
            // Then probably everyone has migrated
            ++t;
            if (t > 1000u) break;
        }
    }
    else {
        for (size_t i = 0u; i < population.size(); ++i)
            if (rnd::bernoulli(p.dispersal)) population[i].disperse();
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

    resources = utl::zeros(2u, 2u);

    if (!p.rdynamics) {

        // LOGISTIC

        // Calculate the resource equilibrium = K (1 - C / r)
        resources[0u][0u] = p.capacity;
        resources[0u][1u] = isburnin ? 0.0 : p.capacity * p.hsymmetry;
        resources[1u][0u] = p.capacity * p.hsymmetry;
        resources[1u][1u] = isburnin ? 0.0 : p.capacity;
        for (size_t hab = 0u; hab < 2u; ++hab) {
            for (size_t res = 0u; res < 2u; ++res) {
                resources[hab][res] *= 1.0 - sumfeed[hab][res] / p.replenish;
                if (resources[hab][res] < 0.0) resources[hab][res] = 0.0;
                assert(resources[hab][res] >= 0.0);
            }
        }
    }
    else {

        // CHEMOSTAT

        // Calculate the resource equilibrium = 1 / (1 + tau C)
        resources[0u][0u] = 1.0;
        resources[0u][1u] = isburnin ? 0.0 : p.hsymmetry;
        resources[1u][0u] = p.hsymmetry;
        resources[1u][1u] = isburnin ? 0.0 : 1.0;
        for (size_t hab = 0u; hab < 2u; ++hab) {
            for (size_t res = 0u; res < 2u; ++res) {
                resources[hab][res] /= 1.0 + p.trenewal * sumfeed[hab][res];
                assert(resources[hab][res] >= 0.0);
            }
        }
    }

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

    const size_t nparents = population.size();

    // Discrete distributions for each habitat
    auto market0 = std::discrete_distribution<size_t>(probs[0u].cbegin(), probs[0u].cend());
    auto market1 = std::discrete_distribution<size_t>(probs[1u].cbegin(), probs[1u].cend());

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
                const size_t dad = hab ? market1(rnd::rng) : market0(rnd::rng);
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
