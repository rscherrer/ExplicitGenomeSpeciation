#ifndef EXPLICITGENOMESPECIATION_NETWORK_H
#define EXPLICITGENOMESPECIATION_NETWORK_H

typedef std::pair<size_t, size_t> Edge;

/// A container for a gene regulatory network
struct Network
{

    Network(const size_t&, const size_t&, const size_t&, const double&,
     const double&, const double&, const Genome&);

    size_t trait;
    size_t nvertices;
    size_t nedges;
    double skewness;

    std::vector<Edge> map;
    std::vector<size_t> loci;
    std::vector<Edge> edges;
    std::vector<double> weights;

    // Makers
    std::vector<Edge> makeMap();
    std::vector<size_t> makeLoci(const Genome&);
    std::vector<Edge> makeEdges();
    std::vector<double> makeWeights(const double&, const double&);

    // Builders
    void initializeNetwork(std::vector<Edge>&, size_t&, std::vector<size_t>&)
     const noexcept;
    void growNetwork(std::vector<Edge>&, size_t&, std::vector<size_t>&)
     const noexcept;
    void sortNetwork(std::vector<Edge>&, const std::vector<size_t>&) const
     noexcept;

};


#endif
