#include "GeneticArchitecture.h"
#include "utils.h"
#include <cassert>
#include <algorithm>
#include <iostream>

typedef std::discrete_distribution<size_t> Discrete;

/// Network constructor
Network::Network(const size_t &character, const size_t &nVertices,
 const size_t &nEdges, const double &skew, const double &shape,
  const double &scale, const Genome &genome) :
    trait(character),
    nvertices(nVertices),
    nedges(nEdges),
    skewness(skew),
    map(makeMap()),
    loci(makeLoci(genome)),
    edges(makeEdges()),
    weights(makeWeights(shape, scale))
{
    assert(map.size() == nedges);
    assert(loci.size() == nvertices);
    assert(edges.size() == nedges);
    assert(weights.size() == nedges);
}

/// // Make a map of pairwise connexions
std::vector<Edge> Network::makeMap()
{

    // There is a total number of vertices
    // There is a total number of edges that can be made
    // The network starts with one edge between 0 and 1
    // Then each new vertex comes in
    // The number of partners of that vertex is sampled
    // The partners of that vertex are sampled without replacement
    // The number of edges still to be made is updated

    assert(nvertices > 1u);
    std::vector<Edge> connexions;
    if (!nedges) return connexions;
    connexions.push_back(std::make_pair(0u, 1u));
    size_t nleft = nedges - 1;

    // For each vertex
    for (size_t vertex = 2u; nleft && vertex < nvertices; ++vertex) {

        // Number of partners
        size_t npartners = nleft;
        const double prob = 1.0 / (nvertices - vertex);
        assert(prob >= 0.0);
        assert(prob <= 1.0);
        if (vertex == nvertices - 1u)
            npartners = rnd::binomial(nleft, prob);

        // Attachment probabilities
        std::vector<double> probs(vertex);
        for (size_t node = 0u; node < vertex; ++node)
            probs[node] = 1.0;

        // For each edge of that vertex
        for (size_t edge = 0u; nleft && edge < npartners; ++edge) {

            // Choose partner without replacement
            Discrete attachment(probs.begin(), probs.end());
            const size_t partner = attachment(rnd::rng);
            assert(partner < vertex);
            connexions.push_back(std::make_pair(partner, vertex));
            probs[partner] = 0.0;
            --nleft;

        }
    }

    assert(nleft == 0u);
    assert(connexions.size() == nedges);

    return connexions;

    //assert(nvertices > 1u);

    //if (n_edges == 0u) return pairs;

    // Start the network
    //std::vector<size_t> degrees(nvertices, 0u);
    //initializeNetwork(network, n_edges, degrees);

    // Grow network by linking preferentially to well-connected nodes
    //growNetwork(network, n_edges, degrees);

    // Relabel node indices after sorting with respect to degree
    //sortNetwork(network, degrees);

}

/// Function to detect the loci underlying a trait
std::vector<size_t> Network::makeLoci(const Genome& genome)
{
    std::vector<size_t> underlying;

    // The current trait must be a field of the network
    // Loop throughout the genome's vector of encoded traits
    // Record all loci that encode the current trait

    for (size_t locus = 0u; locus < genome.nloci; ++locus)
        if (genome.traits[locus] == trait)
            underlying.push_back(locus);

    // Genome has nloci = zero so this function goes nuts
    // There is a problem in how the genome is initialized

    if (underlying.size() != nvertices)
        std::cout << "\nProblematic genome, nloci = " << genome.nloci <<
         ", length of effect sizes = " << genome.effects.size() << "\n\n";

    assert(underlying.size() == nvertices);

    return underlying;
}

/// Function to map the network to the genome
std::vector<Edge> Network::makeEdges()
{

    std::vector<Edge> mapped;

    // Loop through the pairs in the map
    // Use the map id to find the loci in the vector of underlying loci
    // Add theses id to the vector of mapped edges

    for (Edge partners : map) {
        const size_t locus1 = loci[partners.first];
        const size_t locus2 = loci[partners.second];
        mapped.push_back(std::make_pair(locus1, locus2));
    }

    return mapped;

}

/// Function to create a new network
void Network::initializeNetwork(std::vector<Edge> &network, size_t &n_edges,
                                std::vector<size_t> &degrees) const noexcept
{

    if (!network.empty()) network.clear();

    // Create initial network
    network.emplace_back(Edge {0u, 1u});
    degrees[0u] = degrees[1u] = 1u;
    --n_edges;

}


/// Function to iteratively grow a network based on the preferential attachment
/// algorithm
void Network::growNetwork(std::vector<Edge> &network, size_t &n_edges,
 std::vector<size_t> &degrees) const noexcept
{

    // For each vertex in the network...
    for (size_t vertex = 2u; vertex < nvertices && n_edges > 0u; ++vertex) {

        // Assign probability weights to all potential partners
        std::vector<double> probs(vertex);

        for (size_t partner = 0u; partner < vertex; ++partner) {
            probs[partner] = pow(degrees[partner], skewness);
        }

        // Sample the number of attachments of the current vertex
        size_t nAttachments = (vertex == nvertices - 1u ? nedges :
         rnd::binomial(nedges, 1.0 / (nvertices - vertex)));

        // For each new attachment...
        while (nAttachments) {

            // Compute the sum of weights, quit the loop if too small
            double sumProbs = 0.0;
            for (size_t partner = 0u; partner < vertex; ++partner)
                sumProbs += probs[partner];
            if (sumProbs < 1.0) break;

            // Sample the new partner without replacement
            std::discrete_distribution<size_t> attachmentProbs(probs.begin(),
             probs.end());
            size_t partner = attachmentProbs(rnd::rng);

            // Add the new partner to the network
            network.emplace_back(Edge {vertex, partner});
            probs[partner] = 0.0;
            ++degrees[partner];
            ++degrees[vertex];
            --nAttachments;
            --n_edges;
        }
    }
}


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


/// Function to sample interaction weights across edges of a gene regulatory
/// network
std::vector<double> Network::makeWeights(const double &shape,
 const double &scale)
{
    std::vector<double> interweights;
    double sqrtsumsqWeights = 0.0;

    // For each edge in the network...
    for (size_t edge = 0u; edge < nedges; ++edge) {

        // Sample the weight from a two-sided Gamma distribution
        double weight = std::gamma_distribution<double>(shape, scale)(rnd::rng);
        weight = rnd::bernoulli(0.5) ? weight * -1.0 : weight;
        interweights.push_back(weight);

        // Accumulate square rooted sum of squared interaction weights for
        // later normalization
        sqrtsumsqWeights += sqr(weight);
    }

    // Square root the normalizing factor
    sqrtsumsqWeights = sqrtsumsqWeights > 0.0 ? sqrt(sqrtsumsqWeights) : 1.0;

    // Normalize at the end
    for (size_t edge = 0u; edge < nedges; ++edge)
        interweights[edge] /= sqrtsumsqWeights;

    return interweights;
}
