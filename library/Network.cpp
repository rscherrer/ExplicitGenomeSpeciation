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


/// Make a map of pairwise connexions
pVector Network::makeMap()
{

    // There is a total number of vertices
    // There is a total number of edges that can be made
    // The network starts with one edge between 0 and 1
    // Then each new vertex comes in
    // The number of partners of that vertex is sampled
    // The partners of that vertex are sampled without replacement
    // The number of edges still to be made is updated

    assert(nvertices > 1u);

    pVector connexions;
    if (!nedges) return connexions;
    uVector degrees = uzeros(nvertices);

    // First connexion
    connexions.push_back(std::make_pair(0u, 1u));
    ++degrees[0u];
    ++degrees[1u];

    size_t nleft = nedges - 1; // number edges left to make

    // For each vertex
    for (size_t vertex = 2u; nleft && vertex < nvertices; ++vertex) {

        // Sample number of partners
        size_t npartners = nleft;
        const double prob = 1.0 / (nvertices - vertex);
        assert(prob >= 0.0);
        assert(prob <= 1.0);
        if (vertex == nvertices - 1u)
            npartners = rnd::binomial(nleft, prob);

        // Assign attachment probabilities
        dVector probs(vertex);
        for (size_t node = 0u; node < vertex; ++node)
            probs[node] = pow(degrees[node], skewness);

        // For each edge of that vertex
        for (size_t edge = 0u; nleft && edge < npartners; ++edge) {

            if (sum(probs) < 1.0) break;

            // Make a bond without replacement and update degree distribution
            const size_t partner = rnd::sample(probs);
            assert(partner < vertex);
            connexions.push_back(std::make_pair(partner, vertex));
            probs[partner] = 0.0;
            ++degrees[vertex];
            ++degrees[partner];
            assert(degrees[vertex] < nedges);
            assert(degrees[partner] < nedges);
            --nleft;
        }
    }

    assert(nleft == 0u);
    assert(connexions.size() == nedges);

    return connexions;

}


/// Function to detect the loci underlying a trait
uVector Network::makeLoci(const Genome& genome)
{
    uVector underlying;

    // The current trait must be a field of the network
    // Loop throughout the genome's vector of encoded traits
    // Record all loci that encode the current trait

    for (size_t locus = 0u; locus < genome.nloci; ++locus)
        if (genome.traits[locus] == trait)
            underlying.push_back(locus);

    assert(underlying.size() == nvertices);

    return underlying;
}


/// Function to map the network to the genome
pVector Network::makeEdges()
{

    pVector mapped;

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


/// Sample interaction weights
dVector Network::makeWeights(const double &shape,
 const double &scale)
{
    if (shape == 0.0 || scale == 0.0) return zeros(nedges);

    dVector intweights;
    double sss = 0.0; // square rooted sum of squares

    // For each edge in the network...
    for (size_t edge = 0u; edge < nedges; ++edge) {

        // Two-sided Gamma distribution
        const double weight = rnd::bigamma(shape, scale);
        intweights.push_back(weight);
        sss += sqr(weight);
    }

    // Normalize
    sss = sss > 0.0 ? sqrt(sss) : 1.0;
    for (size_t edge = 0u; edge < nedges; ++edge)
        intweights[edge] /= sss;

    return intweights;
}


