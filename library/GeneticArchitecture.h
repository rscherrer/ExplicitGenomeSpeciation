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

    GeneticArchitecture(const ParameterSet &pars) :
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
        snpFreq(pars.getSNPFreq()),
        femHeterogamy(pars.getIsFemaleHeterogamy()),
        recombinationRate(pars.getRecombinationRate()),
        scaleA(pars.getScaleA()),
        scaleD(pars.getScaleD()),
        scaleI(pars.getScaleI()),
        scaleE(pars.getScaleE()),
        chromosomes(makeChromosomes()),
        traits(makeEncodedTraits()),
        locations(makeLocations()),
        effects(makeEffects()),
        dominances(makeDominances()),
        networks(makeNetworks())
    {
        assert(chromosomes.size() == nChromosomes);
        assert(traits.size() == nLoci);
        assert(effects.size() == nLoci);
        assert(dominances.size() == nLoci);
        assert(locations.size() == nLoci);
        assert(networks.size() == 3u);
    }

    MultiNet getNetworks() const { return networks; }

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
    double snpFreq;
    bool femHeterogamy;
    double recombinationRate;

    vecDbl scaleA;
    vecDbl scaleD;
    vecDbl scaleI;
    vecDbl scaleE;

    vecDbl chromosomes;
    vecUns traits;
    vecDbl locations;
    vecDbl effects;
    vecDbl dominances;

    MultiNet networks;

private:

    MultiNet makeNetworks();
    vecDbl makeChromosomes();
    vecUns makeEncodedTraits();
    vecDbl makeLocations();
    vecDbl makeEffects();
    vecDbl makeDominances();

};




#endif
