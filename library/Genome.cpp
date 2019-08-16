#include "GeneticArchitecture.h"
#include "utils.h"
#include <cassert>
#include <algorithm>

/// Genome constructor
Genome::Genome(const std::vector<size_t> &nLociPerTrait,
 const size_t &nLoci, const size_t &nchrom, const double &shape,
  const double &scale) :
    nloci(nLoci),
    chromosomes(makeChromosomes(nchrom)),
    traits(makeEncodedTraits(nLociPerTrait)),
    locations(std::vector<double> { }),
    effects(std::vector<double> { }),
    dominances(std::vector<double> { })
{

    setLocationsEffectSizesAndDominance(nLoci, shape, scale);
    assert(chromosomes.size() == nchrom);
    assert(traits.size() == nloci);
    assert(effects.size() == nloci);
    assert(dominances.size() == nloci);
    assert(locations.size() == nloci);
}


/// Function to make a vector of chromosome sizes
std::vector<double> Genome::makeChromosomes(const size_t &nchrom)
{

    std::vector<double> chromends;

    // Chromosomes all have the same size
    for (size_t chrom = 0u; chrom < nchrom; ++chrom)
        chromends.push_back((chrom + 1.0) / nchrom);

    return chromends;

}


/// Function to randomly assign loci to their encoded traits
std::vector<size_t> Genome::makeEncodedTraits(const std::vector<size_t>
 &nLociPerTrait)
{

    std::vector<size_t> encoded;

    // Make an ordered vector of trait indices
    for (size_t trait = 0u; trait < 3u; ++trait)
        for (size_t locus = 0u; locus < nLociPerTrait[trait]; ++locus)
            encoded.push_back(trait);

    assert(encoded.size() == nloci);

    // Shuffle encoded traits randomly
    std::shuffle(encoded.begin(), encoded.end(), rnd::rng);

    assert(encoded.size() == nloci);

    std::vector<size_t> nvertices {0u, 0u, 0u};
    for (size_t locus = 0u; locus < nloci; ++locus)
        ++nvertices[encoded[locus]];

    for (size_t trait = 0u; trait < 3u; ++trait)
        assert(nvertices[trait] == nLociPerTrait[trait]);

    return encoded;

}

/// Function to sample locations, effect sizes and dominance across the genome
void Genome::setLocationsEffectSizesAndDominance(const size_t &nLoci,
 const double &shape, const double &scale)
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
    for (size_t trait = 0u; trait < 3u; ++trait)
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
