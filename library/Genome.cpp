#include "GeneticArchitecture.h"
#include "Genome.h"
#include "utils.h"
#include <cassert>
#include <algorithm>
#include <iostream>

vecDbl Genome::makeChromosomes(const size_t &nchrom)
{

    vecDbl chromends;

    // Chromosomes all have the same size
    for (size_t chrom = 0u; chrom < nchrom; ++chrom)
        chromends.push_back((chrom + 1.0) / nchrom);

    return chromends;

}

vecUns Genome::makeEncodedTraits(const std::vector<size_t>
 &nLociPerTrait)
{

    vecUns encoded;

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


vecDbl Genome::makeLocations()
{
    vecDbl positions;

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


vecDbl Genome::makeEffects(const double &shape, const double &scale)
{

    if (shape == 0.0 || scale == 0.0) return zeros(nloci);

    vecDbl effectsizes;
    vecDbl sss = {0.0, 0.0, 0.0}; // square rooted sum of squares

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


vecDbl Genome::makeDominances(const double &var)
{

    if (var == 0.0) return zeros(nloci);

    vecDbl coefficients;
    vecDbl sss = {0.0, 0.0, 0.0}; // square rooted sum of squares

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
