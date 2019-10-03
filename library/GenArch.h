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
//#include <algorithm>

typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Network> MultiNet;

// The genetic architecture contains locus-specific details about the
// genotype-phenotype map. It is created using the parameters, but contains
// large arrays of values across the whole genome, and is therefore
// larger than the Param class.

class GenArch {

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

    Network getNetwork(const size_t &i) const { return networks[i]; }

    vecDbl chromosomes;
    vecUns traits;
    vecDbl locations;
    vecDbl effects;
    vecDbl dominances;
    MultiNet networks;

private:

    MultiNet makeNetworks(const Param&);
    vecDbl makeChromosomes(const Param&);
    vecUns makeEncodedTraits(const Param&);
    vecDbl makeLocations(const Param&);
    vecDbl makeEffects(const Param&);
    vecDbl makeDominances(const Param&);

};




#endif
