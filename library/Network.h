#ifndef EXPLICITGENOMESPECIATION_NETWORK_H
#define EXPLICITGENOMESPECIATION_NETWORK_H

#include "types.h"
#include <vector>
#include <stddef.h>

typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Edge> vecEdg; // vector of pairs
typedef std::vector<size_t> vecUns;
typedef std::vector<double> vecDbl;
class Genome;

/// A container for a gene regulatory network
struct Network
{

    Network(const size_t&, const size_t&, const size_t&, const double&,
     const double&, const double&, const Genome&);

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
    vecUns makeLoci(const Genome&);
    vecEdg makeEdges();
    vecDbl makeWeights(const double&, const double&);

};


#endif
