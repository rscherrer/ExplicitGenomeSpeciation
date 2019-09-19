#include "GeneticArchitecture.h"
#include "ParameterSet.h"
#include "Random.h"
#include "Network.h"
#include "Genome.h"
#include "utils.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>

class Network;
class Genome;

/// Function to make a vector of chromosome sizes
vecDbl GeneticArchitecture::makeChromosomes()
{

    vecDbl chromsizes;

    // Chromosomes all have the same size
    for (size_t i = 0u; i < nChromosomes; ++i)
        chromsizes.push_back((i + 1.0) / nChromosomes);

    return chromsizes;

}


/// Function from architecture to call the Genome constructor
Genome GeneticArchitecture::makeGenome()
{
    const Genome gen = Genome(nLociPerTrait, nLoci, nChromosomes,
     effectSizeShape, effectSizeScale, dominanceVariance, femHeterogamy,
      recombinationRate);
    return gen;
}


/// Function to make a vector of interacting partner loci for each trait
MultiNet GeneticArchitecture::makeNetworks()
{
    MultiNet multinet;

    // For each trait
    for (size_t trait = 0u; trait < 3u; ++trait) {

        const size_t nloci = nLociPerTrait[trait];
        const size_t nedges = nEdgesPerTrait[trait];
        const double skew = skewnesses[trait];
        const double intshape = interactionWeightShape;
        const double intscale = interactionWeightScale;

        Network network = Network(trait, nloci, nedges, skew, intshape,
         intscale, genome);

        multinet.push_back(network);
    }

    // The indices in these network maps are indices among the loci underlying
    // a given trait,
    // not absolute loci indices across the genome

    assert(multinet.size() == 3u);

    for (size_t trait = 0u; trait < 3u; ++trait) {
        assert(multinet[trait].map.size() == nEdgesPerTrait[trait]);
    }

    return multinet;
}

vecUns GeneticArchitecture::makeEncodedTraits()
{

    vecUns encoded;

    // Make an ordered vector of trait indices
    for (size_t trait = 0u; trait < 3u; ++trait)
        for (size_t locus = 0u; locus < nLociPerTrait[trait]; ++locus)
            encoded.push_back(trait);

    assert(encoded.size() == nLoci);

    // Shuffle encoded traits randomly
    std::shuffle(encoded.begin(), encoded.end(), rnd::rng);

    assert(encoded.size() == nLoci);

    std::vector<size_t> nvertices {0u, 0u, 0u};
    for (size_t locus = 0u; locus < nLoci; ++locus)
        ++nvertices[encoded[locus]];

    for (size_t trait = 0u; trait < 3u; ++trait)
        assert(nvertices[trait] == nLociPerTrait[trait]);

    return encoded;

}



vecDbl GeneticArchitecture::makeLocations()
{
    vecDbl positions;

    for (size_t locus = 0u; locus < nLoci; ++locus) {
        const double pos = rnd::uniform(1.0);
        positions.push_back(pos);
        assert(pos >= 0.0);
        assert(pos <= 1.0);
    }

    std::sort(positions.begin(), positions.end());

    for (size_t locus = 1u; locus < nLoci; ++locus)
        assert(positions[locus] > positions[locus - 1u]);

    return positions;
}



vecDbl GeneticArchitecture::makeEffects()
{

    const double shape = effectSizeShape;
    const double scale = effectSizeScale;

    if (shape == 0.0 || scale == 0.0) return zeros(nLoci);

    vecDbl effectsizes;
    vecDbl sss = zeros(3u); // square rooted sum of squares

    for (size_t locus = 0u; locus < nLoci; ++locus) {

        const double effect = rnd::bigamma(shape, scale);
        effectsizes.push_back(effect);
        sss[traits[locus]] += sqr(effect);
    }

    for (size_t trait = 0u; trait < 3u; ++trait)
        sss[trait] = sqrt(sss[trait]);

    for (size_t locus = 0u; locus < nLoci; ++locus)
        effectsizes[locus] /= sss[traits[locus]];

    return effectsizes;
}



vecDbl GeneticArchitecture::makeDominances()
{

    if (dominanceVariance == 0.0) return zeros(nLoci);

    vecDbl coefficients;
    vecDbl sss = zeros(3u); // square rooted sum of squares

    for (size_t locus = 0u; locus < nLoci; ++locus) {
        const double dom = rnd::hnormal(dominanceVariance);
        coefficients.push_back(dom);
        sss[traits[locus]] += sqr(dom);
    }

    for (size_t trait = 0u; trait < 3u; ++trait)
        sss[trait] = sqrt(sss[trait]);

    for (size_t locus = 0u; locus < nLoci; ++locus)
        coefficients[locus] /= sss[traits[locus]];

    return coefficients;
}


