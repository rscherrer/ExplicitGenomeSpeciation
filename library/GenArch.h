#ifndef EXPLICITGENOMESPECIATION_GENARCH_H
#define EXPLICITGENOMESPECIATION_GENARCH_H

#include "Param.h"
#include "Random.h"
#include "Network.h"
#include "Types.h"
#include "Utilities.h"
#include <cassert>
#include <stddef.h>
#include <cstddef>

typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Network> MultiNet;

// The genetic architecture contains locus-specific details about the
// genotype-phenotype map. It is created using the parameters, but contains
// large arrays of values across the whole genome, and is therefore
// larger than the Param class.

class GenArch {

    friend class Collector;

public:

    GenArch(const Param &pars) :
        chromosomes(makeChromosomes(pars)),
        traits(makeEncodedTraits(pars)),
        locations(makeLocations(pars)),
        effects(makeEffects(pars)),
        dominances(makeDominances(pars)),
        networks(makeNetworks(pars))
    {
        assert(chromosomes.size() == pars.nchrom);
        assert(traits.size() == pars.nloci);
        assert(effects.size() == pars.nloci);
        assert(dominances.size() == pars.nloci);
        assert(locations.size() == pars.nloci);
        assert(networks.size() == 3u);
    }

    vecDbl chromosomes; // per chromosome
    vecUns traits; // per locus
    vecDbl locations; // per locus
    vecDbl effects; // per locus
    vecDbl dominances; // per locus
    MultiNet networks; // per trait

    // Getters called from outside
    size_t getNetworkSize(const size_t &trait) const
    {
        return networks[trait].map.size();
    }
    Edge getEdge(const size_t &trait, const size_t &edge) const
    {
        return networks[trait].map[edge];
    }
    double getSumWeights(const size_t &trait) const
    {
        double sum = 0.0;
        for (size_t edge = 0u; edge < networks[trait].weights.size(); ++edge)
            sum += networks[trait].weights[edge];
        return sum;
    }

private:

    MultiNet makeNetworks(const Param&);
    vecDbl makeChromosomes(const Param&);
    vecUns makeEncodedTraits(const Param&);
    vecDbl makeLocations(const Param&);
    vecDbl makeEffects(const Param&);
    vecDbl makeDominances(const Param&);

};

#endif
