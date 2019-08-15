#include "GeneticArchitecture.h"
#include "ParameterSet.h"
#include "Random.h"
#include "utils.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>

typedef std::vector<Network> MultiNet;

/// Constructor of genetic architecture
GeneticArchitecture::GeneticArchitecture(const ParameterSet &pars) :
    nChromosomes(pars.getNChromosomes()),
    nLoci(pars.getNLoci()),
    nLociPerTrait(pars.getNLociPerTrait()),
    nEdgesPerTrait(pars.getNEdgesPerTrait()),
    skewnesses(pars.getSkewnesses()),
    effectSizeShape(pars.getEffectSizeShape()),
    effectSizeScale(pars.getEffectSizeScale()),
    interactionWeightShape(pars.getInteractionWeightShape()),
    interactionWeightScale(pars.getInteractionWeightScale()),
    chromosomeSizes(makeChromosomes()),
    genome(makeGenome()),
    networks(makeNetworks())
{
}


/// Function to make a vector of chromosome sizes
std::vector<double> GeneticArchitecture::makeChromosomes()
{

    std::vector<double> chromsizes;

    // Chromosomes all have the same size
    for (size_t i = 0u; i < nChromosomes; ++i)
        chromsizes.push_back((i + 1.0) / nChromosomes);

    return chromsizes;

}


/// Function to make a vector of interacting partner loci for each trait
MultiNet GeneticArchitecture::makeNetworks()
{
    MultiNet multinet;

    // For each trait
    for (size_t trait = 0u; trait < 3u; ++trait)
    {

        Network network = Network(trait, nLociPerTrait[trait],
         nEdgesPerTrait[trait], skewnesses[trait], interactionWeightShape,
          interactionWeightShape, genome);

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


/// Function from architecture to call the Genome constructor
Genome GeneticArchitecture::makeGenome()
{
    const Genome gen = Genome(nLociPerTrait, nLoci, effectSizeShape,
     effectSizeScale);
    return gen;
}
