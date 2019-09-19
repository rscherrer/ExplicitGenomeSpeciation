#ifndef EXPLICITGENOMESPECIATION_NETWORK_H
#define EXPLICITGENOMESPECIATION_NETWORK_H

#include "types.h"
#include <vector>
#include <cassert>
#include <stddef.h>

typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Edge> vecEdg;

/// A container for a gene regulatory network
struct Network
{

    Network(const size_t &character, const size_t &nVertices,
     const size_t &nEdges, const double &skew, const double &shape,
      const double &scale, const vecDbl &locations, const vecUns &traits) :
        trait(character),
        nvertices(nVertices),
        nedges(nEdges),
        skewness(skew),
        map(makeMap()),
        loci(makeUnderlyingLoci(locations, traits)),
        edges(makeEdges()),
        weights(makeWeights(shape, scale))
    {
        assert(map.size() == nedges);
        assert(loci.size() == nvertices);
        assert(edges.size() == nedges);
        assert(weights.size() == nedges);
    }

    size_t trait;
    size_t nvertices;
    size_t nedges;
    double skewness;

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
