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

    if (p.dispersal == 0.0) return;

    if (p.dispersal < 0.1) {

        // Use geometric if rare
        auto getnextmigrant = rnd::iotagap(p.dispersal);
        getnextmigrant.reset(0u);
        size_t mig = 0u;
        for (;;) {
            mig = getnextmigrant(rnd::rng);
            if (mig >= getSize()) break;
            population[mig].disperse();
        }
    }
    /*
    else if (p.dispersal < 0.5) {

        // Use binomial if uncommon
        auto getnmigrants = rnd::binomial(p.dispersal);
        vecUns inds(getSize());
        std::iota(inds.begin(), inds.end(), 0);
        auto getmigrant = rndutils::make_shuffle_sampler(inds.cbegin(), inds.cend());
        size_t nmig = getnmigrants(rnd::rng);
        while (nmig) {
            auto migrant = getmigrant(rnd::rng).second;
            std::clog << migrant << '\n';
            // population[migrant].disperse();
            --nmig;
        }

    }
    */
    else {

        // Use bernoulli if common
        auto ismigrant = rnd::bernoulli(p.dispersal);
        for (size_t i = 0u; i < getSize(); ++i)
            if (ismigrant(rnd::rng)) population[i].disperse();

    }
}

void MetaPop::consume(const Param &p)
{

    // Calculate the sum of feeding efficiencies in each habitat
    Matrix sumfeed = utl::zeros(2u, 2u);
    double sumx = 0.0;
    for (size_t i = 0u; i < population.size(); ++i) {

        sumx += population[i].getEcoTrait();

        const size_t hab = population[i].getHabitat();
        sumfeed[hab][0u] += population[i].getFeeding(0u);

        // Feed only on resource 0 during burnin
        if (!isburnin) sumfeed[hab][1u] += population[i].getFeeding(1u);
    }

    const double meanx = sumx / population.size();

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
    for (size_t i = 0u; i < population.size(); ++i) {
        population[i].feed(resources[population[i].getHabitat()]);
        population[i].classify(meanx);
    }


}

void MetaPop::reproduce(const Param &p, const GenArch &arch)
{

    // Table to count the sexes in both habitats
    sexcounts = utl::uzeros(2u, 2u); // per habitat per sex

    // Determine the length of the season and exit if zero
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

    // Loop through the population and extract habitat and gender
    for (size_t i = 0u; i < getSize(); ++i) {

        const size_t hab = population[i].getHabitat();
        const size_t sex = population[i].getGender();

        // Keep track of census in each sex and each habitat
        ++sexcounts[hab][sex];

        // Record IDs in each sex and each habitat
        if (sex)
            females[hab].push_back(i);
        else
            males[hab].push_back(i);

    }

    // For each habitat...
    for (size_t hab = 0u; hab < 2u; ++hab) {

        females[hab].shrink_to_fit();
        males[hab].shrink_to_fit();

        const size_t nm = sexcounts[hab][0u];
        const size_t nf = sexcounts[hab][1u];

        assert(nm == males[hab].size());
        assert(nf == females[hab].size());

        // If no males or no females, no reproduction and move on
        if (!nm) continue;
        if (!nf) continue;

        vecDbl fit(nm);

        size_t i = 0u;
        double sum = 0.0;
        double ssq = 0.0;

        // Make a vector of probabilities for males
        for (size_t m : males[hab]) {

            fit[i] = population[m].getFitness();

            // Modified in burn-in
            if (isburnin) {
                const double y = population[m].getMatePref();
                fit[i] *= exp(-p.ecosel * utl::sqr(y));
            }

            sum += fit[i];
            ssq += utl::sqr(fit[i]);

            ++i;

        }

        // Set all probabilities to one if all fitnesses are the same
        const double var = ssq / nm - utl::sqr(sum / nm);
        if (var < 1.0E-6) fit = utl::ones(nm);

        // Initialize a mutable discrete distribution with uniform zero-policy
        auto getmale = rnd::mdiscrete();

        // Loop through females
        for (size_t f : females[hab]) {

            // Determine fecundity
            double fecundity = p.birth * population[f].getFitness();

            // Modified during burn-in
            if (isburnin) {
                const double y = population[f].getMatePref();
                fecundity *= exp(-p.ecosel * utl::sqr(y));
            }

            if (fecundity == 0.0) continue;

            // Sample the number of offspring from a Poisson
            auto getclutchsize = rnd::poisson(fecundity);
            size_t noffspring = getclutchsize(rnd::rng);

            if (!noffspring) continue;

            size_t t = seasonend;

            vecDbl probs = fit;

            // While the season is not over and mating hasn't occured
            while (t) {

                getmale.mutate(probs.cbegin(), probs.cend());

                // Sample a male and evaluate
                const size_t idm = getmale(rnd::rng);
                const size_t m = males[hab][idm];
                const double xm = population[m].getEcoTrait();
                auto ismating = rnd::bernoulli(population[f].mate(xm, p));

                // Produce offspring and add them to the population if accepted
                if (ismating(rnd::rng)) {

                    while (noffspring) {

                        population.push_back(Individual(p, arch,
                         population[f], population[m]));

                        --noffspring;

                    }

                    break;

                }

                // Update the distribution to avoid resampling
                probs[idm] = 0.0;
                --t;
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
