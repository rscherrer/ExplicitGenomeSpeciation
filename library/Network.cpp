#include "Network.h"

typedef std::discrete_distribution<size_t> Discrete;

vecEdg Network::makeMap(const Param& p) const
{

    // There is a total number of vertices
    // There is a total number of edges that can be made
    // The network starts with one edge between 0 and 1
    // Then each new vertex comes in
    // The number of partners of that vertex is sampled
    // The partners of that vertex are sampled without replacement
    // The number of edges still to be made is updated

    assert(p.nvertices[trait] > 1u);

    vecEdg connexions;
    if (!p.nedges[trait]) return connexions;
    connexions.reserve(p.nedges[trait]);
    vecUns degrees = utl::uzeros(p.nvertices[trait]);

    // First connexion
    connexions.push_back(std::make_pair(0u, 1u));
    ++degrees[0u];
    ++degrees[1u];

    size_t nleft = p.nedges[trait] - 1u; // number edges left to make

    // For each vertex...
    for (size_t vertex = 2u; nleft && vertex < p.nvertices[trait]; ++vertex) {

        // Sample number of partners
        size_t npartners = nleft;
        const double prob = 1.0 / (p.nvertices[trait] - vertex);
        assert(prob >= 0.0);
        assert(prob <= 1.0);
        if (vertex == p.nvertices[trait] - 1u)
            npartners = rnd::binomial(nleft, prob);

        // Assign attachment probabilities
        vecDbl probs(vertex);
        for (size_t node = 0u; node < vertex; ++node)
            probs[node] = pow(degrees[node], p.skews[trait]);

        // For each edge of that vertex
        for (size_t edge = 0u; nleft && edge < npartners; ++edge) {

            if (utl::sum(probs) < 1.0) break;

            // Make a bond without replacement and update degree distribution
            const size_t partner = rnd::sample(probs);
            assert(partner < vertex);
            connexions.push_back(std::make_pair(partner, vertex));
            probs[partner] = 0.0;
            ++degrees[vertex];
            ++degrees[partner];
            assert(degrees[vertex] < p.nedges[trait]);
            assert(degrees[partner] < p.nedges[trait]);
            --nleft;
        }
    }

    assert(nleft == 0u);
    assert(connexions.size() == p.nedges[trait]);

    return connexions;

}

vecUns Network::makeUnderlyingLoci(const Param &p, const vecUns &traits) const
{
    vecUns underlying;
    underlying.reserve(p.nvertices[trait]);

    // The current trait must be a field of the network
    // Loop throughout the genome's vector of encoded traits
    // Record all loci that encode the current trait

    for (size_t locus = 0u; locus < p.nloci; ++locus)
        if (traits[locus] == trait)
            underlying.push_back(locus);

    assert(underlying.size() == p.nvertices[trait]);

    return underlying;
}

vecEdg Network::makeEdges(const Param &p) const
{

    vecEdg mapped;
    mapped.reserve(p.nedges[trait]);

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

vecDbl Network::makeWeights(const Param &p) const
{
    if (p.interactionshape == 0.0 || p.interactionscale == 0.0)
        return utl::zeros(p.nedges[trait]);

    vecDbl intweights;
    intweights.reserve(p.nedges[trait]);
    double sss = 0.0; // square rooted sum of squares

    // For each edge in the network...
    for (size_t edge = 0u; edge < p.nedges[trait]; ++edge) {

        // Two-sided Gamma distribution
        const double w = rnd::bigamma(p.interactionshape, p.interactionscale);
        intweights.push_back(w);
        sss += utl::sqr(w);
    }

    // Normalize
    sss = sss > 0.0 ? sqrt(sss) : 1.0;
    for (size_t edge = 0u; edge < p.nedges[trait]; ++edge)
        intweights[edge] /= sss;

    return intweights;
}


