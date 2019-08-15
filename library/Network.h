#ifndef EXPLICITGENOMESPECIATION_NETWORK_H
#define EXPLICITGENOMESPECIATION_NETWORK_H

typedef std::pair<size_t, size_t> Edge;

/// A container for a gene regulatory network
struct Network
{

    Network(const size_t&, const size_t&, const double&, const double&,
     const double&);

    // Number of vertices, edges and skewness
    size_t nVertices;
    size_t nEdges;
    double skewness;

    // A map of interacting genes
    std::vector<Edge> map;

    // A vector of interaction weights
    std::vector<double> weights;

    // Functions to make the network
    std::vector<double> makeWeights(const double&, const double&) const
     noexcept;
    std::vector<Edge> makeNetwork(size_t) const noexcept;
    void initializeNetwork(std::vector<Edge>&, size_t&, std::vector<size_t>&)
     const noexcept;
    void growNetwork(std::vector<Edge>&, size_t&, std::vector<size_t>&)
     const noexcept;
    void sortNetwork(std::vector<Edge>&, const std::vector<size_t>&) const
     noexcept;

};


#endif
