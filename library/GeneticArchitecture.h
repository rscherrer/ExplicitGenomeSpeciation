#ifndef EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H
#define EXPLICITGENOMESPECIATION_GENETICARCHITECTURE_H

#include "ParameterSet.h"
#include "Random.h"
#include "Genome.h"
#include "Network.h"
#include "types.h"
#include <vector>
#include <cstddef>
#include <cassert>
#include <stddef.h>


typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Network> MultiNet;


class GeneticArchitecture {

public:

    GeneticArchitecture(const ParameterSet&);

    Genome getGenome() const { return genome; }
    MultiNet getNetworks() const { return networks; }

private:

    size_t nChromosomes;
    size_t nLoci;
    vecUns nLociPerTrait;
    vecUns nEdgesPerTrait;
    vecDbl skewnesses;
    double effectSizeShape;
    double effectSizeScale;
    double interactionWeightShape;
    double interactionWeightScale;
    double dominanceVariance;
    bool femHeterogamy;

    Genome genome;
    MultiNet networks;

    /// Makers
    vecDbl makeChromosomes();
    Genome makeGenome();
    MultiNet makeNetworks();

};

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
    femHeterogamy(pars.getIsFemaleHeterogamy()),
    genome(makeGenome()),
    networks(makeNetworks())
{
    assert(networks.size() == 3u);
}



#endif
