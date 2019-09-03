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
    dominanceVariance(pars.getDominanceVariance()),
    genome(makeGenome()),
    networks(makeNetworks())
{
    assert(networks.size() == 3u);
}


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
     effectSizeShape, effectSizeScale, dominanceVariance);
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



