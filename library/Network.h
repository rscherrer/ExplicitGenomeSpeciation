#ifndef EXPLICITGENOMESPECIATION_NETWORK_H
#define EXPLICITGENOMESPECIATION_NETWORK_H

#include <vector>

typedef std::pair<size_t, size_t> Edge;
typedef std::vector<Edge> pVector; // vector of pairs
typedef std::vector<size_t> uVector;
typedef std::vector<double> dVector;
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

    pVector map;
    uVector loci;
    pVector edges;
    dVector weights;

    // Makers
    pVector makeMap();
    uVector makeLoci(const Genome&);
    pVector makeEdges();
    dVector makeWeights(const double&, const double&);

};


#endif
