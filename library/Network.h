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

    Network(const size_t &character, const Param &pars, const vecUns &traits) :
        trait(character),
        map(makeMap(pars)),
        loci(makeUnderlyingLoci(pars, traits)),
        edges(makeEdges(pars)),
        weights(makeWeights(pars))
    {

        assert(loci.size() == pars.nvertices[character]);

        const size_t realnedges = map.size();

        assert(map.size() == realnedges);
        assert(edges.size() == realnedges);
        assert(weights.size() == realnedges);

        assert(realnedges <= pars.nedges[character]);

        if (realnedges < pars.nedges[character]) {
            throw std::runtime_error("The requested number of edges was not realized for trait " + character);
        }
    }

    size_t trait;

    vecEdg map;
    vecUns loci;
    vecEdg edges;
    vecDbl weights;

    // Makers
    vecEdg makeMap(const Param&) const;
    vecUns makeUnderlyingLoci(const Param&, const vecUns&) const;
    vecEdg makeEdges(const Param&) const;
    vecDbl makeWeights(const Param&) const;

};

#endif
