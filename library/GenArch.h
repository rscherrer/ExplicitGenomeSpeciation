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
        isreset(resetseed(pars.archseed)),
        chromosomes(makeChromosomes(pars)),
        traits(makeEncodedTraits(pars)),
        locations(makeLocations(pars)),
        effects(makeEffects(pars)),
        dominances(makeDominances(pars)),
        networks(makeNetworks(pars))
    {

        // If we want to use a specific genetic architecture
        // We can supply a seed just for the architecture
        // At the end we need to reset the seed back to the original

        rnd::rng.seed(pars.archseed);

        assert(utl::sumu(pars.nvertices) == pars.nloci);
        assert(chromosomes.size() == pars.nchrom);
        assert(traits.size() == pars.nloci);
        assert(effects.size() == pars.nloci);
        assert(dominances.size() == pars.nloci);
        assert(locations.size() == pars.nloci);
        assert(networks.size() == 3u);
        assert(isreset);

        rnd::rng.seed(pars.seed);

    }

    bool isreset;

    vecDbl chromosomes;     // per chromosome
    vecUns traits;          // per locus
    vecDbl locations;       // per locus
    vecDbl effects;         // per locus
    vecDbl dominances;      // per locus
    MultiNet networks;      // per trait

    // Getters called from tests
    size_t getNetworkSize() const
    {
        size_t nedges = 0u;
        for (size_t trait = 0u; trait < 3u; ++trait) {
            nedges += getNetworkSize(trait);
        }
        return nedges;
    }
    size_t getNetworkSize(const size_t &trait) const
    {
        return networks[trait].map.size();
    }
    Edge getEdge(const size_t &trait, const size_t &edge) const
    {
        return networks[trait].map[edge];
    }
    size_t getSumTraits() const
    {
        size_t sum = 0u;
        for (auto x : traits) sum += x;
        return sum;
    }
    double getSumEffects() const
    {
        double sum = 0.0;
        for (auto x : effects) sum += x;
        return sum;
    }
    double getSumDominances() const
    {
        double sum = 0.0;
        for (auto x : dominances) sum += x;
        return sum;
    }
    double getSumWeights(const size_t &trait) const
    {
        double sum = 0.0;
        for (size_t edge = 0u; edge < networks[trait].weights.size(); ++edge)
            sum += networks[trait].weights[edge];
        return sum;
    }

private:

    bool resetseed(const size_t&);
    MultiNet makeNetworks(const Param&);
    vecDbl makeChromosomes(const Param&);
    vecUns makeEncodedTraits(const Param&);
    vecDbl makeLocations(const Param&);
    vecDbl makeEffects(const Param&);
    vecDbl makeDominances(const Param&);

};

#endif
