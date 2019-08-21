#include "GeneticArchitecture.h"
#include "utils.h"
#include <cassert>
#include <algorithm>
#include <iostream>

/// Genome constructor
Genome::Genome(const std::vector<size_t> &nLociPerTrait,
 const size_t &nLoci, const size_t &nchrom, const double &shape,
  const double &scale, const double &domvar) :
    nloci(nLoci),
    chromosomes(makeChromosomes(nchrom)),
    traits(makeEncodedTraits(nLociPerTrait)),
    locations(makeLocations()),
    effects(makeEffects(shape, scale)),
    dominances(makeDominances(domvar))
{
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


dVector Genome::makeLocations()
{
    dVector positions;

    for (size_t locus = 0u; locus < nloci; ++locus) {
        const double pos = rnd::uniform(1.0);
        positions.push_back(pos);
        assert(pos >= 0.0);
        assert(pos <= 1.0);
    }

    std::sort(positions.begin(), positions.end());

    for (size_t locus = 1u; locus < nloci; ++locus)
        assert(positions[locus] > positions[locus - 1u]);

    return positions;
}


dVector Genome::makeEffects(const double &shape, const double &scale)
{

    if (shape == 0.0 || scale == 0.0) return zeros(nloci);

    dVector effectsizes;
    dVector sss = {0.0, 0.0, 0.0}; // square rooted sum of squares

    for (size_t locus = 0u; locus < nloci; ++locus) {

        const double effect = rnd::bigamma(shape, scale);
        effectsizes.push_back(effect);
        sss[traits[locus]] += sqr(effect);
    }

    for (size_t trait = 0u; trait < 3u; ++trait)
        sss[trait] = sqrt(sss[trait]);

    for (size_t locus = 0u; locus < nloci; ++locus)
        effectsizes[locus] /= sss[traits[locus]];

    return effectsizes;
}


dVector Genome::makeDominances(const double &var)
{

    if (var == 0.0) return zeros(nloci);

    dVector coefficients;
    dVector sss = {0.0, 0.0, 0.0}; // square rooted sum of squares

    for (size_t locus = 0u; locus < nloci; ++locus) {
        const double dom = rnd::hnormal(var);
        coefficients.push_back(dom);
        sss[traits[locus]] += sqr(dom);
    }

    for (size_t trait = 0u; trait < 3u; ++trait)
        sss[trait] = sqrt(sss[trait]);

    for (size_t locus = 0u; locus < nloci; ++locus)
        coefficients[locus] /= sss[traits[locus]];

    return coefficients;
}
