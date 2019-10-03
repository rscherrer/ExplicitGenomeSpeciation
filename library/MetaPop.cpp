#include "MetaPop.h"


// Initialization

Crowd MetaPop::populate(const Param &p, const GenArch &arch)
{

    // Generate a pool of individuals
    const size_t n = utl::sumu(p.demesizes);
    Crowd indivs;
    indivs.reserve(n);

    for (size_t ind = 0u; ind < n; ++ind)
        indivs.push_back(Individual(p, arch));

    return indivs;

}

// Main function

void MetaPop::disperse(const Param &p) // only if not burnin
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
    Matrix consumed = utl::matzeros(2u, 2u);
    for (auto &ind : population) {
        const size_t hab = ind.getHabitat();
        consumed[hab][0u] += ind.getFeeding(0u);
        if (!isburnin) consumed[hab][1u] += ind.getFeeding(1u);
    }

    assert(consumed[0u][0u] >= 0.0);
    assert(consumed[0u][1u] >= 0.0);
    assert(consumed[1u][0u] >= 0.0);
    assert(consumed[1u][1u] >= 0.0);

    // Update the resource equilibria
    Matrix resources = utl::matzeros(2u, 2u);
    for (size_t hab = 0u; hab < 2u; ++hab) {
        for (size_t res = 0u; res < 2u; ++res) {
            if (res == 1u && isburnin) {
                resources[hab][res] == 0.0;
            }
            else {
                resources[hab][res] = p.capacity;
                if (hab == res) resources[hab][res] *= p.hsymmetry;
                resources[hab][res] *= (1.0 - consumed[hab][res] / p.replenish);
            }
        }
    }

    // Could the resources be negative? If too much consumption?
    assert(resources[0u][0u] >= 0.0);
    assert(resources[0u][1u] >= 0.0);
    assert(resources[1u][0u] >= 0.0);
    assert(resources[1u][1u] >= 0.0);

    // Assign individual fitness and ecotypes
    for (auto &ind : population) ind.feed(resources);
}

void MetaPop::reproduce(const Param &p)
{
    // Table to count the sexes in both habitats
    MatUns sexcounts = utl::matuzeros(2u, 2u);

    // Males' fitness determine their encounter probabilities
    Matrix probs = utl::matzeros(2u, population.size());
    for (size_t i = 0u; i < population.size(); ++i) {

        const size_t sex = population[i].getGender();
        const size_t hab = population[i].getHabitat();
        ++sexcounts[sex][hab];

        if (!sex) {

            probs[hab][i] = population[i].getFitness();

            // Modify mating success during burnin
            if (isburnin) {
                const double y = population[i].getMatePref();
                probs[hab][i] *= exp(-ecosel * utl::sqr(y));
            }
        }
    }

    if (!utl::sumu(sexcounts[0u])) return; // exit if no males
    if (!utl::sumu(sexcounts[1u])) return; // exit if no females

    // Determine the length of the mating season
    const size_t seasonend = rnd::geometric(p.matingcost);

    if (!seasonend) return;

    // For each individuals...
    for (auto &fem : population) {

        // Is it a female?
        if (!fem.getGender()) continue;

        // Are there males around?
        if (!sexcounts[0u][fem.getHabitat()]) continue;

        // Progress throughout the mating season
        size_t timeleft = seasonend;
        while (timeleft) {

            --timeleft;

            // And encounters males one at a time with replacement
            const size_t pick = rnd::sample(probs[fem.getHabitat()]);
            auto &candidate = population[pick];

            // If the female accepts to mate
            if (fem.accept(candidate)) {

                // Determine fecundity
                double fecundity = p.birth * fem.getFitness();

                // Modify fecundity during burnin
                if (isburnin) {
                    const size_t y = fem.getMatePref();
                    fecundity *= exp(-ecosel * utl::sqr(y));
                }

                // Sample clutch size
                size_t noffspring = rnd::poisson(fecundity);
                while (noffspring) {

                    // Give birth
                    population.push_back(Individual(fem, candidate));
                    --noffspring;
                }

                // End the mating season if female has mated
                break;
            }
        }
    }
}

void MetaPop::survive()
{
    // Sample survival for each individual
    size_t nsurvivors = 0u;
    for (auto &ind : population) {
        const bool isalive = rnd::bernoulli(survival);
        ind.survives(isalive);
        if (isalive) ++nsurvivors;
    } // offspring always

    // Sort the population to put dead individuals at the end
    std::sort (population.begin(), population.end(), burry);

    // Resize the population to get rid of dead ones
    population.resize(nsurvivors);
}

void MetaPop::cycle() // don't forget the burnin
{
    // Dispersal
    if (!isburnin) disperse();

    // Consumption
    consume();

    // Reproduction
    reproduce();

    // Survival
    survive();

}

void MetaPop::analyze(const GenArch& arch)
{
    stats.reset(t, arch);
    stats.analyze(pops, arch);
    stats.setEcoIsolation();
    stats.setSpatialIsolation(pops);
    stats.setMatingIsolation(pops, matingcost, sexsel);
}


// Getters

size_t MetaPop::getDemeSize(const size_t &d) const
{
    size_t n = 0u;
    for (auto &ind : population)
        if (ind.getHabitat() == d) ++n;
    return n;
}

size_t MetaPop::getNOffspring(const size_t &p) const
{
    return pops[p].getNOffspring();
}

size_t MetaPop::getSumEcotypes(const size_t &p) const
{
    return pops[p].getSumEcotypes();
}

size_t MetaPop::getSumFemEcotypes() const
{
    double sum = 0.0;
    for (auto &pop : pops) sum += pop.getSumFemEcotypes();
    return sum;
}

double MetaPop::getResource(const size_t &p, const size_t &r) const
{
    return pops[p].getResource(r);
}

double MetaPop::getVarP(const size_t &trait, const size_t &eco) const
{
    return stats.getVarP(trait, eco);
}

double MetaPop::getSsqPhe(const size_t &trait, const size_t &eco) const
{
    return stats.getSsqPhe(trait, eco);
}

double MetaPop::getSumPhe(const size_t &trait, const size_t &eco) const
{
    return stats.getSumPhe(trait, eco);
}

double MetaPop::getSumTrait(const size_t &trait, const size_t &pop) const
{
    return pops[pop].getSumTrait(trait);
}


// Resetters used in tests

void MetaPop::consume() // metapopulation-level fitness/ecotype resetter
{
    for (auto &pop : pops)
        pop.consume();
}

void MetaPop::sortSexes() // metapopulation level sex-vectors resetter
{
    for (auto &pop : pops)
        pop.sortSexes();
}

void MetaPop::resetEcoTraits(const size_t &p, const double &x)
{
    pops[p].resetEcoTraits(x, ecosel, maxfeed);
}

void MetaPop::resetMatePrefs(const size_t &p, const double &y)
{
    pops[p].resetMatePrefs(y);
}

void MetaPop::resetEcotypes(const size_t &p, const size_t &e)
{
    pops[p].resetEcotypes(e);
}

void MetaPop::resetGenders(const size_t &p, const bool &sex)
{
    pops[p].resetGenders(sex);
    pops[p].sortSexes();
}
