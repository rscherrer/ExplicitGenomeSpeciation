#ifndef EXPLICITGENOMESPECIATION_NETWORK_H
#define EXPLICITGENOMESPECIATION_NETWORK_H

#include "Param.h"
#include "Types.h"
#include "Utilities.h"
#include "Random.h"
#include <cassert>
#include <stddef.h>

typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Edge> vecEdg;

// A class for a gene regulatory network. One network underlies one of the
// three traits of the simulation. Networks are generated using a
// preferential attachment algorithm, and are part of the bigger class GenArch.
// Could it be a nested class?

class Network
{

public:

    Network(const size_t &character, const Param &pars, const vecDbl &locations,
     const vecUns &traits) :
        trait(character),
        map(makeMap()),
        loci(makeUnderlyingLoci(locations, traits)),
        edges(makeEdges()),
        weights(makeWeights(pars.effectshape, pars.effectscale))
    {
        assert(map.size() == pars.nedges[character]);
        assert(loci.size() == pars.nvertices[character]);
        assert(edges.size() == pars.nedges[character]);
        assert(weights.size() == pars.nedges[character]);
    }

    size_t trait;

    vecEdg map;
    vecUns loci;
    vecEdg edges;
    vecDbl weights;

    // Makers
    vecEdg makeMap();
    vecUns makeUnderlyingLoci(const vecDbl&, const vecUns&);
    vecEdg makeEdges();
    vecDbl makeWeights(const double&, const double&);

};



#endif
