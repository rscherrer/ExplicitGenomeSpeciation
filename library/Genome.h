#ifndef EXPLICITGENOMESPECIATION_GENOME_H
#define EXPLICITGENOMESPECIATION_GENOME_H

#include "Types.h"
#include <vector>
#include <cassert>
#include <stddef.h>

// Remove this file

struct Genome
{

    Genome(const vecUns &nLociPerTrait, const size_t &nLoci,
     const size_t &nchrom, const double &shape, const double &scale,
      const double &domvar, const double &recombination,
       const bool &heterogamy) :
        nloci(nLoci),
        chromosomes(makeChromosomes(nchrom)),
        traits(makeEncodedTraits(nLociPerTrait)),
        locations(makeLocations()),
        effects(makeEffects(shape, scale)),
        dominances(makeDominances(domvar)),
        recombrate(recombination),
        femgamy(heterogamy)
    {
        assert(chromosomes.size() == nchrom);
        assert(traits.size() == nloci);
        assert(effects.size() == nloci);
        assert(dominances.size() == nloci);
        assert(locations.size() == nloci);
    }

    size_t nloci;
    vecDbl chromosomes;

    vecUns traits;
    vecDbl locations;
    vecDbl effects;
    vecDbl dominances;

    double recombrate;
    bool femgamy;

    vecDbl makeChromosomes(const size_t&);
    vecUns makeEncodedTraits(const vecUns&);
    vecDbl makeLocations();
    vecDbl makeEffects(const double&, const double&);
    vecDbl makeDominances(const double& = 1.0);

};

#endif
