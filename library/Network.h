#ifndef EXPLICITGENOMESPECIATION_NETWORK_H
#define EXPLICITGENOMESPECIATION_NETWORK_H

#include "Param.h"

#include "Utilities.h"
#include "Random.h"
#include <cassert>
#include <stddef.h>

typedef std::pair<size_t, size_t> Edge;

// A class for a gene regulatory network. One network underlies one of the
// three traits of the simulation. Networks are generated using a
// preferential attachment algorithm, and are part of the bigger class GenArch.
// Could it be a nested class?

class Network
{

public:

    Network(const size_t &character, const Param &pars, const std::vector<size_t> &traits) :
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

    std::vector<Edge> map;
    std::vector<size_t> loci;
    std::vector<Edge> edges;
    std::vector<double> weights;

    // Makers
    std::vector<Edge> makeMap(const Param&) const;
    std::vector<size_t> makeUnderlyingLoci(const Param&, const std::vector<size_t>&) const;
    std::vector<Edge> makeEdges(const Param&) const;
    std::vector<double> makeWeights(const Param&) const;

};

#endif
