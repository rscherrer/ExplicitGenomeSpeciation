#ifndef EXPLICITGENOMESPECIATION_GENOME_H
#define EXPLICITGENOMESPECIATION_GENOME_H

#include "types.h"
#include <vector>
#include <stddef.h>

typedef std::vector<double> vecDbl;
typedef std::vector<size_t> vecUns;

/// A container of constant genetic features across the genome
struct Genome
{

    Genome(const vecUns&, const size_t&, const size_t&,
     const double&, const double&, const double&, const bool&);

    size_t nloci;
    vecDbl chromosomes;

    vecUns traits;
    vecDbl locations;
    vecDbl effects;
    vecDbl dominances;

    bool femgamy;

    // A vector of locus epistatic interaction partners and weights
    // std::vector<std::vector<std::pair<size_t, double> > > interactions;

    // Member functions
    vecDbl makeChromosomes(const size_t&);
    vecUns makeEncodedTraits(const vecUns&);
    vecDbl makeLocations();
    vecDbl makeEffects(const double&, const double&);
    vecDbl makeDominances(const double& = 1.0);

    void setLocationsEffectSizesAndDominance(const size_t&,
     const double&, const double&);


};

#endif
