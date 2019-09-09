    //assert(nvertices > 1u);

    //if (n_edges == 0u) return pairs;

    // Start the network
    //std::vector<size_t> degrees(nvertices, 0u);
    //initializeNetwork(network, n_edges, degrees);

    // Grow network by linking preferentially to well-connected nodes
    //growNetwork(network, n_edges, degrees);

    // Relabel node indices after sorting with respect to degree
    //sortNetwork(network, degrees);

//-------------

/// Function to compare edges in a network with respect to their degrees
bool edgeCompare(const Edge &x, const Edge &y) noexcept
{
    if (x.first == y.first) {
        return (x.second < y.second);
    }
    else {
        return (x.first < y.first);
    }
}

/// Function to sort the edges in the network by degree
void Network::sortNetwork(std::vector<Edge> &network,
 const std::vector<size_t> &degrees) const noexcept
{

    // Compute the ranks of all vertices with respect to their degrees
    std::vector<size_t> ranks(nvertices, 0u);
    for (size_t vertex = 0u; vertex < nvertices - 1u; ++vertex) {
        for (size_t othervertex = vertex + 1u; othervertex < nvertices;
             ++othervertex)
        {
            if (degrees[othervertex] > degrees[vertex])
            {
                ++ranks[vertex];
            }
            else
            {
                ++ranks[othervertex];
            }
        }
    }

    // For each edge, swap partners to make sure the first partner always have
    // the lower rank
    for (Edge &edge : network) {
        edge.first = ranks[edge.first];
        edge.second = ranks[edge.second];
        if (edge.first > edge.second)
        {
            std::swap(edge.first, edge.second);
        }
    }

    // Sort the edges
    std::sort(network.begin(), network.end(), edgeCompare);

}

