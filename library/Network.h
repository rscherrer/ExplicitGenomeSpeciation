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

class Network
{

public:

    Network(const size_t &character, const Param &pars,
            const std::vector<size_t> &traits) :
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
            std::string msg = "The requested number of edges was not realized for trait ";
            msg += static_cast<char>(character);
            throw std::runtime_error(msg);
        }
    }

    size_t trait;

    std::vector<Edge> map;
    std::vector<size_t> loci;
    std::vector<Edge> edges;
    std::vector<double> weights;

    // Getters
    bool isConnected() const
    {
        // Get the degree of each locus
        std::vector<size_t> degrees(loci.size(), 0u);

        for (size_t e = 0u; e < map.size(); ++e) {
            ++degrees[map[e].first];
            ++degrees[map[e].second];
        }

        // Verify that no locus has a degree of zero
        for (size_t l = 0u; l < degrees.size(); ++l)
            if (!degrees[l]) return false;

        return true;
    }

    // Makers
    std::vector<Edge> makeMap(const Param&) const;
    std::vector<size_t> makeUnderlyingLoci(
            const Param&,
            const std::vector<size_t>&
    ) const;
    std::vector<Edge> makeEdges(const Param&) const;
    std::vector<double> makeWeights(const Param&) const;

};

#endif
