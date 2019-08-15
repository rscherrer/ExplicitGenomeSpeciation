#include "GeneticArchitecture.h"
#include "utils.h"
#include <cassert>
#include <algorithm>

/// Genome constructor
Genome::Genome(const size_t &nTraits, const std::vector<size_t> &nLociPerTrait,
 const size_t &nLoci, const double &shape, const double &scale) :
    traits(makeEncodedTraits(nTraits, nLociPerTrait)),
    locations(std::vector<double> { 0.0 }),
    effects(std::vector<double> { 0.0 }),
    dominances(std::vector<double> { 0.0 })
{

    locations.pop_back();
    effects.pop_back();
    dominances.pop_back();

    // Sample locations, effect sizes and dominance coefficients across the
    // genome
    setLocationsEffectSizesAndDominance(nTraits, nLoci, shape, scale);

    assert(traits.size() == nLoci);
    assert(effects.size() == nLoci);
    assert(dominances.size() == nLoci);
    assert(locations.size() == nLoci);
}

/// Function to randomly assign loci to their encoded traits
std::vector<size_t> Genome::makeEncodedTraits(const size_t &nTraits,
 const std::vector<size_t> &nLociPerTrait) const noexcept
{

    std::vector<size_t> encoded;

    // Make an ordered vector of trait indices
    for (size_t trait = 0u; trait < nTraits; ++trait) {
        for (size_t locus = 0u; locus < nLociPerTrait[trait]; ++locus) {
            encoded.push_back(trait);
            assert(encoded.back() <= 2u);
        }
    }

    // Shuffle encoded traits randomly
    std::shuffle(encoded.begin(), encoded.end(), rnd::rng);

    return encoded;

}
/// Function to sample locations, effect sizes and dominance across the genome
void Genome::setLocationsEffectSizesAndDominance(const size_t &nTraits,
 const size_t &nLoci, const double &shape, const double &scale)
{

    // Prepare squared roots of sums of squared effect sizes and dominance
    // coefficients
    std::vector<double> sqrtsumsqEffectSizes {0.0, 0.0, 0.0};
    std::vector<double> sqrtsumsqDominanceCoeffs {0.0, 0.0, 0.0};

    // Sample locations, effect sizes and dominance coefficients
    for (size_t locus = 0u; locus < nLoci; ++locus)
    {

        // Locations are sampled uniformly
        locations.push_back(rnd::uniform(1.0));

        assert(locations.back() > 0.0);
        assert(locations.back() < 1.0);

        // Effect sizes are sampled from a two-sided Gamma distribution
        double effectsize = std::gamma_distribution<double>(shape, scale)(
         rnd::rng);
        effectsize = rnd::bernoulli(0.5) ? effectsize * -1.0 : effectsize;
        effects.push_back(effectsize);

        // Squared effect sizes are accumulated for normalizing
        sqrtsumsqEffectSizes[traits[locus]] += sqr(effectsize);

        // Dominance coefficients are sampled from a one-sided normal
        // distribution
        const double dominance = fabs(rnd::normal(0.0, 1.0));
        assert(dominance > 0.0);
        dominances.push_back(dominance);

        // Squared dominance coefficients are accumulated for normalizing
        sqrtsumsqDominanceCoeffs[traits[locus]] += sqr(dominance);

    }

    // Sort gene locations by increasing order
    std::sort(locations.begin(), locations.end());

    for (size_t locus = 1u; locus < nLoci; ++locus)
        assert(locations[locus] > locations[locus - 1u]);

    // Take the square roots of the sums of squares for normalization
    for (size_t trait = 0u; trait < nTraits; ++trait)
    {
        sqrtsumsqEffectSizes[trait] = sqrt(sqrtsumsqEffectSizes[trait]);
        sqrtsumsqDominanceCoeffs[trait] = sqrt(sqrtsumsqDominanceCoeffs[trait]);

        assert(sqrtsumsqEffectSizes[trait] > 0.0);
        assert(sqrtsumsqDominanceCoeffs[trait] > 0.0);
    }

    // Normalize effect sizes and dominance coefficients across the genome
    for (size_t locus = 0u; locus < nLoci; ++locus)
    {
        effects[locus] /= sqrtsumsqEffectSizes[traits[locus]];
        dominances[locus] /= sqrtsumsqDominanceCoeffs[traits[locus]];
    }
}
