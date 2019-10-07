#include "GenArch.h"

vecDbl GenArch::makeChromosomes(const Param &p)
{

    vecDbl chromends;
    chromends.reserve(p.nchrom);

    // Chromosomes all have the same size
    for (size_t chrom = 0u; chrom < p.nchrom; ++chrom)
        chromends.push_back((chrom + 1.0) / p.nchrom);

    return chromends;

}

vecUns GenArch::makeEncodedTraits(const Param &p)
{

    vecUns encoded;
    encoded.reserve(p.nloci);

    // Make an ordered vector of trait indices
    for (size_t trait = 0u; trait < 3u; ++trait)
        for (size_t locus = 0u; locus < p.nvertices[trait]; ++locus)
            encoded.push_back(trait);

    assert(encoded.size() == p.nloci);

    // Shuffle encoded traits randomly
    std::shuffle(encoded.begin(), encoded.end(), rnd::rng);

    assert(encoded.size() == p.nloci);

    vecUns nvertices {0u, 0u, 0u};
    for (size_t locus = 0u; locus < p.nloci; ++locus)
        ++nvertices[encoded[locus]];

    for (size_t trait = 0u; trait < 3u; ++trait)
        assert(nvertices[trait] == p.nvertices[trait]);

    return encoded;

}



vecDbl GenArch::makeLocations(const Param &p)
{
    vecDbl positions;
    positions.reserve(p.nloci);

    for (size_t locus = 0u; locus < p.nloci; ++locus)
        positions.push_back(rnd::uniform(1.0));

    std::sort(positions.begin(), positions.end());

    for (size_t locus = 1u; locus < p.nloci; ++locus)
        assert(positions[locus] > positions[locus - 1u]);

    return positions;
}


vecDbl GenArch::makeEffects(const Param &p)
{

    if (p.effectshape == 0.0 || p.effectscale == 0.0)
        return utl::zeros(p.nloci);

    vecDbl effectsizes;
    effectsizes.reserve(p.nloci);
    vecDbl sss = utl::zeros(3u); // square rooted sum of squares

    for (size_t locus = 0u; locus < p.nloci; ++locus) {

        const double effect = rnd::bigamma(p.effectshape, p.effectscale);
        effectsizes.push_back(effect);
        sss[traits[locus]] += utl::sqr(effect);
    }

    for (size_t trait = 0u; trait < 3u; ++trait)
        sss[trait] = sqrt(sss[trait]);

    for (size_t locus = 0u; locus < p.nloci; ++locus)
        effectsizes[locus] /= sss[traits[locus]];

    return effectsizes;
}


vecDbl GenArch::makeDominances(const Param &p)
{

    if (p.dominancevar == 0.0) return utl::zeros(p.nloci);

    vecDbl coefficients;
    coefficients.reserve(p.nloci);
    vecDbl sss = utl::zeros(3u); // square rooted sum of squares

    for (size_t locus = 0u; locus < p.nloci; ++locus) {
        const double dom = rnd::hnormal(p.dominancevar);
        coefficients.push_back(dom);
        sss[traits[locus]] += utl::sqr(dom);
    }

    for (size_t trait = 0u; trait < 3u; ++trait)
        sss[trait] = sqrt(sss[trait]);

    for (size_t locus = 0u; locus < p.nloci; ++locus)
        coefficients[locus] /= sss[traits[locus]];

    return coefficients;
}

MultiNet GenArch::makeNetworks(const Param &p)
{
    MultiNet multinet;
    multinet.reserve(3u);

    for (size_t trait = 0u; trait < 3u; ++trait)
        multinet.push_back(Network(trait, p, traits));

    // The indices in these network maps are indices among the loci underlying
    // a given trait,
    // not absolute loci indices across the genome

    assert(multinet.size() == 3u);

    return multinet;
}
