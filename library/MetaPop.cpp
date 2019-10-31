#include "MetaPop.h"


// Initialization

Crowd MetaPop::populate(const Param &p, const GenArch &arch)
{

    // Generate a pool of individuals
    const size_t n = utl::sumu(p.demesizes);
    size_t n0 = p.demesizes[0u];
    assert(n0 <= n);
    Crowd indivs;
    indivs.reserve(n);

    for (size_t ind = 0u; ind < n; ++ind) {
        indivs.push_back(Individual(p, arch));
        if (ind >= n0) indivs[ind].disperse();
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

    /*
    if (p.dispersal < 0.1) {

        // Use geometric if rare
        auto getnextmigrant = rnd::iotagap(p.dispersal);
        getnextmigrant.reset(0u);
        size_t mig = 0u;
        for (;;) {
            mig = getnextmigrant(rnd::rng);
            if (mig > getSize()) break;
            population[mig].disperse();
        }
    }
    else if (p.dispersal < 0.5) {

        // Use binomial if uncommon
        auto getnmigrants = rnd::binomial(p.dispersal);
        vecUns inds(getSize());
        std::iota(inds.begin(), inds.end(), 0);
        auto getmigrant = rnd::samplenr(inds.begin(), inds.end());
        size_t nmig = getnmigrants(rnd::rng);
        while (nmig) {
            const size_t migrant = getmigrant(rnd::rng);
            population[migrant].disperse();
            --nmig;
        }

    }
    else {

        // Use bernoulli if common
        auto ismigrant = rnd::bernoulli(p.dispersal);
        for (size_t i = 0u; i < getSize(); ++i)
            if (ismigrant(rnd::rng)) population[i].disperse();

    }
    */

    // Use bernoulli if common
    auto ismigrant = rnd::bernoulli(p.dispersal);
    for (size_t i = 0u; i < getSize(); ++i)
        if (ismigrant(rnd::rng)) population[i].disperse();

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
    sexcounts = utl::uzeros(2u, 2u); // per habitat per sex

    // Determine the length of the season and exit if zero
    // Loop through the population and extract habitat and gender
    // Record IDs in each sex and each habitat
    // Keep track of census in each sex and each habitat
    // For each habitat...
    // If no males or no females, no reproduction and move on to next habitat
    // Make a vector of probabilities for males (including burnin condition)
    // Set all probabilities to zero if all fitnesses are the same
    // Loop through females
    // Initialize a mutable discrete distribution with uniform zero-policy
    // Determine fecundity (including burnin condition)
    // Sample the number of offspring from a Poisson
    // While the season is not over and mating hasn't occured
    // Sample a male and evaluate
    // Produce offspring and add them to the population if accepted
    // Then move on to the next female
    // Or progress through mating season and mutate the distribution

    auto getseasonend = rnd::geometric(p.matingcost);
    const size_t seasonend = getseasonend(rnd::rng);

    if (!seasonend) return;

    std::vector<vecUns> females;
    std::vector<vecUns> males;

    for (size_t hab = 0u; hab < 2u; ++hab) {
        vecUns fem;
        vecUns mal;
        fem.reserve(getSize());
        mal.reserve(getSize());
        females.push_back(fem);
        males.push_back(mal);
    }

    for (size_t i = 0u; i < getSize(); ++i) {

        const size_t hab = population[i].getHabitat();
        const size_t sex = population[i].getGender();

        ++sexcounts[hab][sex];

        if (sex)
            females[hab].push_back(i);
        else
            males[hab].push_back(i);

    }

    for (size_t hab = 0u; hab < 2u; ++hab) {

        females[hab].shrink_to_fit();
        males[hab].shrink_to_fit();

        const size_t nm = sexcounts[hab][0u];
        const size_t nf = sexcounts[hab][1u];

        assert(nm == males[hab].size());
        assert(nf == females[hab].size());

        if (!nm) continue;
        if (!nf) continue;

        vecDbl fit(nm);

        size_t i = 0u;
        double sum = 0.0;
        double ssq = 0.0;

        for (size_t m : males[hab]) {

            fit[i] = population[m].getFitness();

            if (isburnin) {
                const double y = population[m].getMatePref();
                fit[i] *= exp(-p.ecosel * utl::sqr(y));
            }

            sum += fit[i];
            ssq += utl::sqr(fit[i]);

            ++i;

        }

        const double var = ssq / nm - utl::sqr(sum / nm);

        if (var < 1.0E-6) fit = utl::ones(nm);

        auto getmale = rnd::mdiscrete();

        for (size_t f : females[hab]) {

            double fecundity = p.birth * population[f].getFitness();

            if (isburnin) {
                const double y = population[f].getMatePref();
                fecundity *= exp(-p.ecosel * utl::sqr(y));
            }

            if (fecundity == 0.0) continue;

            auto getclutchsize = rnd::poisson(fecundity);
            size_t noffspring = getclutchsize(rnd::rng);

            if (!noffspring) continue;

            size_t t = seasonend;

            vecDbl probs = fit;

            while (t) {

                getmale.mutate(probs.cbegin(), probs.cend());

                const size_t idm = getmale(rnd::rng);
                const size_t m = males[hab][idm];
                const double xm = population[m].getEcoTrait();
                auto ismating = rnd::bernoulli(population[f].mate(xm, p));

                if (ismating(rnd::rng)) {

                    while (noffspring) {

                        population.push_back(Individual(p, arch,
                         population[f], population[m]));

                        --noffspring;

                    }

                    break;

                }

                probs[idm] = 0.0;

                --t;
            }
        }
    }

    ///////////////

    /*

    // Males' fitness determine their encounter probabilities
    Matrix probs = utl::zeros();
    vecDbl var = utl::zeros(2u);
    vecDbl avg = utl::zeros(2u);

    for (size_t i = 0u; i < population.size(); ++i) {

        const size_t sex = population[i].getGender();
        const size_t hab = population[i].getHabitat();
        ++sexcounts[hab][sex];

        // If it's a male...
        if (!sex) {

            // Add its fitness to the vector of probabilities
            probs[hab][i] = population[i].getFitness();

            // Modify mating success during burnin
            if (isburnin) {
                const double y = population[i].getMatePref();
                probs[hab][i] *= exp(-p.ecosel * utl::sqr(y));
            }

            // Keep track of the variance
            var[hab] += utl::sqr(probs[hab][i]);
            avg[hab] += probs[hab][i];
        }
    }

    if (!(sexcounts[0u][0u] + sexcounts[1u][0u])) return; // exit if no males
    if (!(sexcounts[0u][1u] + sexcounts[1u][1u])) return; // exit if no females

    // Determine the length of the mating season
    auto getseasonend = rnd::geometric(p.matingcost);
    const size_t seasonend = getseasonend(rnd::rng);

    if (!seasonend) return;

    const size_t nparents = getSize();

    // One discrete distribution per habitat
    std::vector<rnd::discrete> markets;
    std::vector<bool> isequal = { false, false };

    for (size_t hab = 0u; hab < 2u; ++hab) {

        // Probabilities are all zero if fitnesses are too similar
        var[hab] /= nparents;
        avg[hab] /= nparents;
        var[hab] -= utl::sqr(avg[hab]);
        isequal[hab] = var[hab] <= 1.0E-6;
        if (var[hab] <= 1.0E-6) probs[hab] = utl::zeros(nparents);

        // Make distribution
        markets.push_back(rnd::discrete(probs[hab].cbegin(), probs[hab].cend()));
    }

    // One random uniform distribution in case no need for discrete
    auto equalmales = rnd::random(0, nparents - 1u);

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
                size_t dad;
                if (isequal[hab])
                    dad = equalmales(rnd::rng);
                else
                    dad = markets[hab](rnd::rng);

                const double maletrait = population[dad].getEcoTrait();
                const double prob = population[mom].mate(maletrait, p);
                auto ismating = rnd::bernoulli(prob);

                // If the female accepts to mate
                if (ismating(rnd::rng)) {

                    // Determine fecundity
                    double fecundity = p.birth * population[mom].getFitness();

                    // Modify fecundity during burnin
                    if (isburnin) {
                        const size_t y = population[mom].getMatePref();
                        fecundity *= exp(-p.ecosel * utl::sqr(y));
                    }

                    // Sample clutch size
                    size_t noffspring = 0u;
                    if (fecundity > 0.0) {
                        auto getclutchsize = rnd::poisson(fecundity);
                        noffspring = getclutchsize(rnd::rng);
                    }
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

    */
}

// Lambda for removing dead individuals
auto burry = [&](Individual ind) -> bool
{
    return !ind.isalive();
};

void MetaPop::survive(const Param &p)
{
    // Sample survival for each individual
    size_t nsurv = 0u;
    auto issurvivor = rnd::bernoulli(p.survival);
    for (size_t i = 0u; i < population.size(); ++i) {
        population[i].survive(issurvivor(rnd::rng));
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
