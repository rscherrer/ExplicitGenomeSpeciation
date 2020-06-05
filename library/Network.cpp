#include "Network.h"

bool checkedges(const size_t &e, const size_t &n)
{
    return e <= n * (n - 1u) / 2u; // max number of edges given number of nodes
}

std::vector<Edge> Network::makeMap(const Param& p) const
{

    // There is a total number of vertices
    // There is a total number of edges that can be made
    // The network starts with one edge between 0 and 1
    // Then each new vertex comes in
    // The number of partners of that vertex is sampled
    // The partners of that vertex are sampled without replacement
    // The number of edges still to be made is updated

    assert(p.nvertices[trait] > 1u);
    assert(checkedges(p.nedges[trait], p.nvertices[trait]));

    std::vector<Edge> connexions;
    connexions.reserve(p.nedges[trait]);
    std::vector<size_t> degrees = std::vector<size_t>(p.nvertices[trait], 0u);

    // First connexion
    connexions.push_back(std::make_pair(0u, 1u));
    ++degrees[0u];
    ++degrees[1u];

    size_t eleft = p.nedges[trait] - 1u; // number of edges left to make
    if (!p.nedges[trait]) return connexions;

    // For each vertex...
    for (size_t vertex = 2u; eleft && vertex < p.nvertices[trait]; ++vertex) {

        // Sample number of partners from a binomial
        // This procedure conditions on the average of the degree distribution
        // being nedges / nvertices

        const size_t nleft = p.nvertices[trait] - vertex;
        const double prob = 1.0 / nleft;
        assert(prob >= 0.0);
        assert(prob <= 1.0);

        size_t npartners;
        if (vertex == p.nvertices[trait] - 1u) {
            npartners = eleft;
        }
        else {
            auto seekpartners = rnd::binomial(eleft - nleft, prob);
            npartners = 1u + seekpartners(rnd::rng);
        }

        // Assign attachment probabilities
        std::vector<double> probs(vertex);
        for (size_t node = 0u; node < vertex; ++node)
            probs[node] = pow(degrees[node], p.skews[trait]);

        // For each edge of that vertex
        for (size_t edge = 0u; eleft && edge < npartners; ++edge) {

            if (utl::sum(probs) < 1.0) break;

            // Sampling multiple targets without replacement
            // Use Hanno's mutable discrete distribution

            auto getpartner = rnd::discrete(probs.cbegin(), probs.cend());
            const size_t partner = getpartner(rnd::rng);
            assert(partner < vertex);
            connexions.push_back(std::make_pair(partner, vertex));
            probs[partner] = 0.0;
            ++degrees[vertex];
            ++degrees[partner];
            assert(degrees[vertex] < p.nedges[trait]);
            assert(degrees[partner] < p.nedges[trait]);
            --eleft;
        }
    }

    // If the number of edges required is too large the algorithm may fail
    // to reach a big enough network -- break in this case

    if (connexions.size() != p.nedges[trait])
        throw std::runtime_error("Required sized network could not be made.");

    assert(eleft == 0u);

    return connexions;

}

std::vector<size_t> Network::makeUnderlyingLoci(const Param &p, const std::vector<size_t> &traits) const
{
    std::vector<size_t> underlying(p.nvertices[trait]);

    // The current trait must be a field of the network
    // Loop throughout the genome's vector of encoded traits
    // Record all loci that encode the current trait

    size_t i = 0u;
    for (size_t locus = 0u; locus < p.nloci; ++locus) {
        if (traits[locus] == trait) {
            underlying[i] = locus;
            ++i;
        }
    }

    std::shuffle(underlying.begin(), underlying.end(), rnd::rng);

    assert(underlying.size() == p.nvertices[trait]);

    return underlying;
}

std::vector<Edge> Network::makeEdges(const Param &p) const
{

    std::vector<Edge> mapped(p.nedges[trait]);

    // Loop through the pairs in the map
    // Use the map id to find the loci in the vector of underlying loci
    // Add theses id to the vector of mapped edges

    for (size_t e = 0u; e < map.size(); ++e) {
        const size_t locus1 = loci[map[e].first];
        const size_t locus2 = loci[map[e].second];
        mapped[e] = std::make_pair(locus1, locus2);
    }

    return mapped;

}

std::vector<double> Network::makeWeights(const Param &p) const
{
    if (p.interactionshape == 0.0)
        return std::vector<double>(p.nedges[trait], 0.0);
    if (p.interactionscale == 0.0)
        return std::vector<double>(p.nedges[trait], 1.0);

    std::vector<double> intweights(p.nedges[trait]);
    double sss = 0.0; // square rooted sum of squares

    // Interaction weights are sampled from a two-sided Gamma distribution
    auto getweight = rnd::gamma(p.interactionshape, p.interactionscale);
    auto isflipped = rnd::bernoulli(0.5);

    // For each edge in the network...
    for (size_t edge = 0u; edge < p.nedges[trait]; ++edge) {

        // Two-sided Gamma distribution
        double w = getweight(rnd::rng);
        if (isflipped(rnd::rng)) w *= -1.0;
        intweights[edge] = w;
        sss += utl::sqr(w);
    }

    // Normalize
    sss = sss > 0.0 ? sqrt(sss) : 1.0;
    assert(sss > 0.0);
    for (size_t edge = 0u; edge < p.nedges[trait]; ++edge)
        intweights[edge] /= sss;

    return intweights;
}


