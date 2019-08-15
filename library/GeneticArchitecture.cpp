#include "GeneticArchitecture.h"
#include "ParameterSet.h"
#include "Random.h"
#include "utils.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>


/// Constructor of genetic architecture
GeneticArchitecture::GeneticArchitecture(const ParameterSet &pars) :
    nTraits(pars.getNTraits()),
    nChromosomes(pars.getNChromosomes()),
    nLoci(pars.getNLoci()),
    nLociPerTrait(pars.getNLociPerTrait()),
    nEdgesPerTrait(pars.getNEdgesPerTrait()),
    skewnesses(pars.getSkewnesses()),
    effectSizeShape(pars.getEffectSizeShape()),
    effectSizeScale(pars.getEffectSizeScale()),
    interactionWeightShape(pars.getInteractionWeightShape()),
    interactionWeightScale(pars.getInteractionWeightScale()),
    chromosomeSizes(makeChromosomeSizes()),
    traitNetworks(makeTraitNetworks()),
    genome(makeGenome())
{}


/// Function to make a vector of chromosome sizes
std::vector<double> GeneticArchitecture::makeChromosomeSizes() const noexcept
{

    std::vector<double> chromsizes;

    // Chromosomes all have the same size
    for (size_t i = 0u; i < nChromosomes; ++i)
        chromsizes.push_back((i + 1.0) / nChromosomes);

    return chromsizes;

}


/// Function to make a vector of interacting partner loci for each trait
std::vector<Network> GeneticArchitecture::makeTraitNetworks() const
 noexcept
{
    std::vector<Network> networks;

    // For each trait
    for (size_t trait = 0u; trait < nTraits; ++trait)
    {

        // Make a network map (a vector of edges) for the current trait using
        // the preferential
        // attachment algorithm
        Network network = Network(nLociPerTrait[trait], nEdgesPerTrait[trait],
         skewnesses[trait], interactionWeightShape, interactionWeightShape);

        networks.push_back(network);
    }

    // The indices in these network maps are indices among the loci underlying
    // a given trait,
    // not absolute loci indices across the genome

    assert(networks.size() == nTraits);
    for (size_t trait = 0u; trait < nTraits; ++trait)
        assert(networks[trait].map.size() == nEdgesPerTrait[trait]);

    return networks;
}


/// Function from architecture to call the Genome constructor
Genome GeneticArchitecture::makeGenome() const noexcept
{
    const Genome gen = Genome(nTraits, nLociPerTrait, nLoci, effectSizeShape,
     effectSizeScale);
    return gen;
}
