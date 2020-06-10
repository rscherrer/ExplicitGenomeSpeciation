#include "Collector.h"

std::vector<Locus> Collector::emptyloci(const GenArch &arch) const
{
    std::vector<Locus> loci;
    loci.reserve(arch.traits.size());
    for (size_t locus = 0u; locus < arch.traits.size(); ++locus) {
        loci.push_back(Locus(locus, arch.traits[locus]));
    }
    return loci;
}

std::vector<Connexion> Collector::emptyconnexions(const GenArch &arch) const
{
    std::vector<Connexion> connexions;
    const size_t nedges = arch.getNetworkSize();

    if (!nedges) return connexions;

    connexions.reserve(nedges);
    for (size_t trait = 0u; trait < 3u; ++trait) {
        for (size_t edge = 0u; edge < arch.getNetworkSize(trait); ++edge) {
            const size_t i = arch.networks[trait].edges[edge].first;
            const size_t j = arch.networks[trait].edges[edge].second;
            connexions.push_back(Connexion(edge, i, j, trait));
        }
    }

    assert(connexions.size() == nedges);
    return connexions;
}

double Xst(const std::vector<double> &v, const std::vector<size_t> &n)
{

    assert(v.size() == 3u);
    assert(n.size() == 3u);

    // If there is only one ecotype, there is no divergence
    if (!n[0u]) return 0.0;
    if (!n[1u]) return 0.0;

    // Check that variances are positive
    assert(v[0u] >= 0.0);
    assert(v[1u] >= 0.0);
    assert(v[2u] >= 0.0);

    // If there is no variance, there is no divergence
    if (v[2u] == 0.0) return 0.0;

    double xst = 1.0 - (n[0u] * v[0u] + n[1u] * v[1u]) / (n[2u] * v[2u]);

    utl::correct(xst, 0.0, 1.0E-5); // make sure we do not have negative values

    assert(xst >= 0.0);
    assert(xst <= 1.0);    

    return xst;
}

void Collector::analyze(const MetaPop &m, const Param &p, const GenArch &a)
{

    // Indices for readability
    const size_t tot = 2u; // across ecotypes
    const size_t all = 3u; // across genotypes
    const size_t aa = 0u;
    const size_t Aa = 1u;
    const size_t AA = 2u;
    const size_t bb = aa;
    const size_t Bb = Aa;
    const size_t BB = AA;

    // Reset
    counts = utl::uzeros(3u, 3u); // per habitat per ecotype
    means = utl::zeros(3u, 3u, 3u); // per trait per habitat per ecotype
    varG = utl::zeros(3u, 3u); // per trait per ecotype
    varP = utl::zeros(3u, 3u); // per trait per ecotype
    varA = utl::zeros(3u, 3u); // per trait per ecotype
    varN = utl::zeros(3u, 3u); // per trait per ecotype
    varD = std::vector<double>(3u, 0.0); // per trait
    varI = std::vector<double>(3u, 0.0); // per trait
    varT = std::vector<double>(3u, 0.0); // per trait
    Pst = std::vector<double>(3u, 0.0); // per trait
    Gst = std::vector<double>(3u, 0.0); // per trait
    Qst = std::vector<double>(3u, 0.0); // per trait
    Cst = std::vector<double>(3u, 0.0); // per trait
    Fst = std::vector<double>(3u, 0.0); // per trait
    EI = 0.0;
    SI = 0.0;
    RI = 0.0;

    // Create tables for counts
    // per trait per habitat per ecotype
    std::vector<std::vector<std::vector<double> > > sumgen = utl::zeros(3u, 3u, 3u);
    std::vector<std::vector<std::vector<double> > > sumphe = utl::zeros(3u, 3u, 3u);
    std::vector<std::vector<std::vector<double> > > ssqgen = utl::zeros(3u, 3u, 3u);
    std::vector<std::vector<std::vector<double> > > ssqphe = utl::zeros(3u, 3u, 3u);

    // Calculate sums within ecotypes and habitats
    for (size_t i = 0u; i < m.population.size(); ++i) {

        // Count densities
        const size_t eco = m.getEcotype(i);
        const size_t hab = m.getHabitat(i);
        ++counts[hab][eco];

        // Accumulate trait and genetic values
        for (size_t trait = 0u; trait < 3u; ++trait) {
            const double gen = m.population[i].getGenValue(trait);
            const double phe = m.population[i].getTraitValue(trait);
            sumgen[trait][hab][eco] += gen;
            sumphe[trait][hab][eco] += phe;
            ssqgen[trait][hab][eco] += utl::sqr(gen);
            ssqphe[trait][hab][eco] += utl::sqr(phe);
        }
    }

    // Counts across ecotypes and across habitats
    utl::marginalize(counts);

    assert(counts[0u][0u] + counts[0u][1u] + counts[1u][0u] + counts[1u][1u] == m.getSize());

    // For each trait...
    for (size_t trait = 0u; trait < 3u; ++trait) {

        // Calculate sums across ecotypes and across habitats
        utl::marginalize(sumgen[trait]);
        utl::marginalize(sumphe[trait]);
        utl::marginalize(ssqgen[trait]);
        utl::marginalize(ssqphe[trait]);

        // Check sums of squares
        for (size_t hab = 0u; hab < 3u; ++hab) {
            for (size_t eco = 0u; eco < 3u; ++eco) {
                assert(ssqgen[trait][hab][eco] >= 0.0);
                assert(ssqphe[trait][hab][eco] >= 0.0);
            }
        }

        // Calculate means within and across ecotypes and habitats
        means[trait] = utl::dividemat(sumgen[trait], counts);

    }

    // From now use only counts across habitats
    std::vector<size_t> ecounts = counts[2u]; // per ecotype
    assert(ecounts[tot] = m.population.size());
    assert(ecounts[tot] == ecounts[0u] + ecounts[1u]);
    // per trait per ecotype
    std::vector<std::vector<double> > esumgen = utl::zeros(3u, 3u);
    std::vector<std::vector<double> > esumphe = utl::zeros(3u, 3u);
    std::vector<std::vector<double> > essqgen = utl::zeros(3u, 3u);
    std::vector<std::vector<double> > essqphe = utl::zeros(3u, 3u);

    for (size_t trait = 0u; trait < 3u; ++trait) {
        esumgen[trait] = sumgen[trait][2u];
        esumphe[trait] = sumphe[trait][2u];
        essqgen[trait] = ssqgen[trait][2u];
        essqphe[trait] = ssqphe[trait][2u];
    }

    // Initialize genome-wide variances in heterozygosity
    std::vector<double> varS = std::vector<double>(3u, 0.0); // per trait

    // Scan the genome and compute locus-specific variances
    for (size_t l = 0u; l < genomescan.size(); ++l) {

        // Reset locus
        genomescan[l].varG = std::vector<double>(3u, 0.0); // per ecotype
        genomescan[l].varP = std::vector<double>(3u, 0.0); // per ecotype
        genomescan[l].varA = std::vector<double>(3u, 0.0); // per ecotype
        genomescan[l].varN = std::vector<double>(3u, 0.0); // per ecotype
        genomescan[l].varD = 0.0;
        genomescan[l].varI = 0.0;
        genomescan[l].varZ = 0.0;
        genomescan[l].varX = 0.0;
        genomescan[l].Pst = 0.0;
        genomescan[l].Gst = 0.0;
        genomescan[l].Qst = 0.0;
        genomescan[l].Cst = 0.0;
        genomescan[l].Fst = 0.0;
        genomescan[l].alpha = 0.0;
        genomescan[l].beta = std::vector<double>(3u, 0.0); // per genotype
        genomescan[l].meang = 0.0;
        genomescan[l].freq = 0.0;

        const double locusvarE = p.locusE[genomescan[l].trait];

        // Create count tables
        std::vector<std::vector<size_t> > gcounts = utl::uzeros(3u, 4u); // per ecotype per genotype
        std::vector<std::vector<double> > gsumgen = utl::zeros(3u, 4u); // per ecotype per genotype
        std::vector<std::vector<double> > gssqgen = utl::zeros(3u, 4u); // per ecotype per genotype

        // Calculate sums within ecotypes and genotypes
        for (size_t i = 0u; i < m.population.size(); ++i) {

            const size_t eco = m.population[i].getEcotype();
            const double gen = m.population[i].getLocusValue(genomescan[l].id);
            const size_t zyg = m.population[i].getZygosity(genomescan[l].id);

            ++gcounts[eco][zyg];
            gsumgen[eco][zyg] += gen;
            gssqgen[eco][zyg] += utl::sqr(gen);

        }

        // Calculate sums over ecotypes and genotypes
        utl::marginalize(gcounts);
        utl::marginalize(gsumgen);
        utl::marginalize(gssqgen);

        // Check sums of squares
        for (size_t eco = 0u; eco < 3u; ++eco) {
            for (size_t zyg = 0u; zyg < 4u; ++zyg) {
                assert(gssqgen[eco][zyg] >= 0.0);
            }
        }

        // Calculate allele frequencies within and across ecotypes
        std::vector<double> allfreq = std::vector<double>(3u, 0.0);
        for (size_t eco = 0u; eco < 3u; ++eco) {
            if (ecounts[eco]) {
                allfreq[eco] = gcounts[eco][AA] + 0.5 * gcounts[eco][Aa];
                allfreq[eco] /= ecounts[eco];
            }
            assert(allfreq[eco] >= 0.0);
            assert(allfreq[eco] <= 1.0);
        }

        genomescan[l].freq = allfreq[tot];

        // Regress genetic values on allele counts

        // meanq = mean allele count

        const double meanq = 2.0 * allfreq[tot];
        assert(meanq >= 0.0);
        assert(meanq <= 2.0);

        // meang = mean genetic value

        const double meang = gsumgen[tot][all] / ecounts[tot];
        genomescan[l].meang = meang;

        // varq = variance in allele count

        // varq = E(q2) - E(q)2
        // using ssq(q) = 4 nAA + nAa

        double varq;
        varq = (4.0 * gcounts[tot][AA] + gcounts[tot][Aa]) / ecounts[tot];
        varq -= utl::sqr(meanq);
        assert(varq >= 0.0);

        genomescan[l].varZ = varq;

        // covqg = covariance between allele count and genetic value

        // covqg = E(qg) - E(q)E(g)
        // using sum(qg) = sum_Aa(g) + 2 sum_AA(g)

        double covqg;
        covqg = (2.0 * gsumgen[tot][AA] + gsumgen[tot][Aa]) / ecounts[tot];
        covqg -= meanq * meang;

        // alpha = average effect of a mutation

        const double alpha = varq == 0.0 ? 0.0 : covqg / varq;
        genomescan[l].alpha = alpha;

        // Variance in gene expression

        const double scaleD = p.scaleD[genomescan[l].trait];

        double varx = gcounts[tot][Aa] * utl::sqr(scaleD * a.dominances[l]);
        varx += gcounts[tot][AA];
        varx += gcounts[tot][aa];
        varx /= ecounts[tot];
        double meanx = gcounts[tot][Aa] * scaleD * a.dominances[l];
        meanx += gcounts[tot][AA];
        meanx -= gcounts[tot][aa];
        meanx /= ecounts[tot];
        varx -= utl::sqr(meanx);
        assert(varx >= 0.0);
        genomescan[l].varX = varx;

        // Calculate genotype-specific statistics

        // delta = deviation in genetic value due to dominance
        // expec = expected genetic value under additivity
        // beta = breeding value
        // = deviation from population mean genetic value due to additivity

        // delta(AA) = meang(AA) - expec(AA)
        // expec(AA) = meang(all) + beta(AA)
        // beta(AA) = alpha (AA - meanq)

        std::vector<double> gbeta = std::vector<double>(3u, 0.0); // per genotype
        std::vector<double> gexpec = std::vector<double>(3u, 0.0); // per genotype
        std::vector<double> gmeans = std::vector<double>(3u, 0.0); // per genotype
        std::vector<double> gdelta = std::vector<double>(3u, 0.0) ; // per genotype
        for (size_t zyg : { aa, Aa, AA }) {
            if (gcounts[tot][zyg]) {
                gbeta[zyg] = alpha * (zyg - meanq);
                gexpec[zyg] = meang + gbeta[zyg];
                gmeans[zyg] = gsumgen[tot][zyg] / gcounts[tot][zyg];
                gdelta[zyg] = gmeans[zyg] - gexpec[zyg];
            }
        }

        genomescan[l].beta = gbeta;

        // Calculate variance components within and across ecotypes
        for (size_t eco = 0u; eco < 3u; ++eco) {

            if (!ecounts[eco]) continue;

            // Genetic variance = V(g)
            genomescan[l].varG[eco] = gssqgen[eco][all] / ecounts[eco];
            const double emean = gsumgen[eco][all] / ecounts[eco];
            genomescan[l].varG[eco] -= utl::sqr(emean);
            utl::correct(genomescan[l].varG[eco], 0.0, 1.0E-5);
            assert(genomescan[l].varG[eco] >= 0.0);

            // Phenotypic variance
            genomescan[l].varP[eco] = genomescan[l].varG[eco] + locusvarE;
            utl::correct(genomescan[l].varP[eco], 0.0, 1.0E-5);
            assert(genomescan[l].varP[eco] >= 0.0);

            // Additive variance = V(beta)
            // Variance in breeding values
            genomescan[l].varA[eco] = gcounts[eco][aa] * utl::sqr(gbeta[aa]);
            genomescan[l].varA[eco] += gcounts[eco][Aa] * utl::sqr(gbeta[Aa]);
            genomescan[l].varA[eco] += gcounts[eco][AA] * utl::sqr(gbeta[AA]);
            genomescan[l].varA[eco] /= ecounts[eco];
            double meanb = gcounts[eco][aa] * gbeta[aa];
            meanb += gcounts[eco][Aa] * gbeta[Aa];
            meanb += gcounts[eco][AA] * gbeta[AA];
            meanb /= ecounts[eco];
            genomescan[l].varA[eco] -= utl::sqr(meanb);
            utl::correct(genomescan[l].varA[eco], 0.0, 1.0E-5);
            assert(genomescan[l].varA[eco] >= 0.0);

            // Simpler equation for the whole pop
            // Stopped using it because makes bugs harder to detect
            // (asserts may fail just because of this assumption)
            // genomescan[l].varA[tot] = utl::sqr(alpha) * varq;

            // Contribute to genome-wide additive variance
            varA[genomescan[l].trait][eco] += genomescan[l].varA[eco];

            // Non-additive variance = V(gamma)
            genomescan[l].varN[eco] = gexpec[aa] * gsumgen[eco][aa];
            genomescan[l].varN[eco] += gexpec[Aa] * gsumgen[eco][Aa];
            genomescan[l].varN[eco] += gexpec[AA] * gsumgen[eco][AA];
            genomescan[l].varN[eco] *= -2.0;
            genomescan[l].varN[eco] += gcounts[eco][aa] * utl::sqr(gexpec[aa]);
            genomescan[l].varN[eco] += gcounts[eco][Aa] * utl::sqr(gexpec[Aa]);
            genomescan[l].varN[eco] += gcounts[eco][AA] * utl::sqr(gexpec[AA]);
            genomescan[l].varN[eco] += gssqgen[eco][all];
            genomescan[l].varN[eco] /= ecounts[eco];
            double meandev = gsumgen[eco][all];
            meandev -= gcounts[eco][aa] * gexpec[aa];
            meandev -= gcounts[eco][Aa] * gexpec[Aa];
            meandev -= gcounts[eco][AA] * gexpec[AA];
            meandev /= ecounts[eco];
            genomescan[l].varN[eco] -= utl::sqr(meandev);
            utl::correct(genomescan[l].varN[eco], 0.0, 1.0E-5);
            assert(genomescan[l].varN[eco] >= 0.0);

            // Contribute to genome-wide non-additive variance
            varN[genomescan[l].trait][eco] += genomescan[l].varN[eco];

        }

        // Dominance variance (across ecotypes only) = V(delta)
        genomescan[l].varD = gcounts[tot][aa] * utl::sqr(gdelta[aa]);
        genomescan[l].varD += gcounts[tot][Aa] * utl::sqr(gdelta[Aa]);
        genomescan[l].varD += gcounts[tot][AA] * utl::sqr(gdelta[AA]);
        genomescan[l].varD /= ecounts[tot];
        utl::correct(genomescan[l].varD, 0.0, 1.0E-5);
        assert(genomescan[l].varD >= 0.0);

        // Contribute to genome-wide dominance variance
        varD[genomescan[l].trait] += genomescan[l].varD;

        // Interaction variance (across ecotypes only) = V(epsilon)
        genomescan[l].varI = gmeans[aa] * gsumgen[tot][aa];
        genomescan[l].varI += gmeans[Aa] * gsumgen[tot][Aa];
        genomescan[l].varI += gmeans[AA] * gsumgen[tot][AA];
        genomescan[l].varI *= -2.0;
        genomescan[l].varI += gcounts[tot][aa] * utl::sqr(gmeans[aa]);
        genomescan[l].varI += gcounts[tot][Aa] * utl::sqr(gmeans[Aa]);
        genomescan[l].varI += gcounts[tot][AA] * utl::sqr(gmeans[AA]);
        genomescan[l].varI += gssqgen[tot][all];
        genomescan[l].varI /= ecounts[tot];
        utl::correct(genomescan[l].varI, 0.0, 1.0E-5);
        assert(genomescan[l].varI >= 0.0);

        // Contribute to genome-wide interaction variance
        varI[genomescan[l].trait] += genomescan[l].varI;

        // Variance in within-ecotype heterozygosity
        double lvarS = ecounts[0u] * utl::sqr(allfreq[0u]);
        lvarS += ecounts[1u] * utl::sqr(allfreq[1u]);
        lvarS /= ecounts[tot];
        lvarS -= utl::sqr(allfreq[tot]);

        varS[genomescan[l].trait] += lvarS;

        // Average heterozygosity within ecotypes
        double h = ecounts[0u] * allfreq[0u] * (1.0 - allfreq[0u]);
        h += ecounts[1u] * allfreq[1u] * (1.0 - allfreq[1u]);
        h /= ecounts[tot];
        assert(h >= 0.0);
        assert(h <= 1.0);

        // Heterozygosity across ecotypes
        double H = allfreq[tot] * (1.0 - allfreq[tot]);
        assert(H >= 0.0);
        assert(H <= 1.0);

        varT[genomescan[l].trait] += H;

        // Calculate divergence metrics
        genomescan[l].Pst = Xst(genomescan[l].varP, ecounts);
        genomescan[l].Gst = Xst(genomescan[l].varG, ecounts);
        genomescan[l].Qst = Xst(genomescan[l].varA, ecounts);
        genomescan[l].Cst = Xst(genomescan[l].varN, ecounts);
        genomescan[l].Fst = H > 0.0 ? 1.0 - h / H : 0.0;

        // Check divergence metrics
        assert(genomescan[l].Pst <= 1.0);
        assert(genomescan[l].Gst <= 1.0);
        assert(genomescan[l].Qst <= 1.0);
        assert(genomescan[l].Cst <= 1.0);
        assert(genomescan[l].Fst <= 1.0);
        assert(genomescan[l].Pst >= 0.0);
        assert(genomescan[l].Gst >= 0.0);
        assert(genomescan[l].Qst >= 0.0);
        assert(genomescan[l].Cst >= 0.0);

    }

    // Calculate genome-wide variances for each trait
    for (size_t trait = 0u; trait < 3u; ++trait) {

        // Within and across ecotypes
        for (size_t eco = 0u; eco < 3u; ++eco) {

            if (!ecounts[eco]) continue;

            // Genetic variance
            varG[trait][eco] = essqgen[trait][eco] / ecounts[eco];
            varG[trait][eco] -= utl::sqr(esumgen[trait][eco] / ecounts[eco]);
            utl::correct(varG[trait][eco], 0.0, 1.0E-5);
            assert(varG[trait][eco] >= 0.0);

            // Phenotypic variance
            varP[trait][eco] = essqphe[trait][eco] / ecounts[eco];
            varP[trait][eco] -= utl::sqr(esumphe[trait][eco] / ecounts[eco]);
            utl::correct(varP[trait][eco], 0.0, 1.0E-5);
            assert(varP[trait][eco] >= 0.0);

        }

        // Divergence metrics
        Pst[trait] = Xst(varP[trait], ecounts);
        Gst[trait] = Xst(varG[trait], ecounts);
        Qst[trait] = Xst(varA[trait], ecounts);
        Cst[trait] = Xst(varN[trait], ecounts);
        if (varT[trait] != 0.0) Fst[trait] = varS[trait] / varT[trait];

        // Check divergence metrics
        assert(Pst[trait] <= 1.0);
        assert(Gst[trait] <= 1.0);
        assert(Qst[trait] <= 1.0);
        assert(Cst[trait] <= 1.0);
        assert(Fst[trait] <= 1.0);
        assert(Pst[trait] >= 0.0);
        assert(Gst[trait] >= 0.0);
        assert(Qst[trait] >= 0.0);
        assert(Cst[trait] >= 0.0);
    }

    // Speciation metrics

    // Ecological isolation
    EI = Pst[0u];

    // Spatial isolation    
    double norm = counts[0u][0u] + counts[0u][1u];
    norm *= counts[1u][0u] + counts[1u][1u];
    norm *= counts[0u][0u] + counts[1u][0u];
    norm *= counts[0u][1u] + counts[1u][1u];
    assert(norm >= 0.0);
    if (norm != 0.0) {
        SI = counts[0u][0u] * counts[1u][1u];
        SI -= counts[0u][1u] * counts[1u][0u];
        SI /= sqrt(norm);
    }
    assert(SI >= -1.0);
    assert(SI <= 1.0);

    // Mating isolation

    // Sort males and females in the population
    std::vector<size_t> males;
    std::vector<size_t> females;
    males.reserve(m.population.size());
    females.reserve(m.population.size());
    for (size_t i =0u; i < m.population.size(); ++i) {
        const size_t sex = m.population[i].getGender();
        if (sex) females.push_back(i); else males.push_back(i);
    }
    males.shrink_to_fit();
    females.shrink_to_fit();

    if (females.size() && males.size()) {

        std::vector<std::vector<size_t> > crosses = utl::uzeros(2u, 2u);
        size_t ntrials = p.ntrials;

        // Sample from a distribution of males and a distribution of females
        auto femalepool = rnd::random(0u, females.size() - 1u);
        auto malepool = rnd::random(0u, males.size() - 1u);

        // Sample many pairs of males and females with replacement
        while (ntrials) {

            const size_t fem = females[femalepool(rnd::rng)];
            const size_t mal = males[malepool(rnd::rng)];

            // See if the female accepts the male or not
            const double maletrait = m.population[mal].getTraitValue(0u);
            const double prob = m.population[fem].mate(maletrait, p);
            auto ismating = rnd::bernoulli(prob);

            const size_t ecof = m.population[fem].getEcotype();
            const size_t ecom = m.population[mal].getEcotype();

            // Count homogamic and heterogamic crosses
            if (ismating(rnd::rng)) ++crosses[ecof][ecom];

            --ntrials;
        }

        // Compute mating isolation statistic
        norm = crosses[0u][0u] + crosses[0u][1u];
        norm *= crosses[1u][0u] + crosses[1u][1u];
        norm *= crosses[0u][0u] + crosses[1u][0u];
        norm *= crosses[0u][1u] + crosses[1u][1u];
        assert(norm >= 0.0);

        if (norm != 0.0) {
            RI = crosses[0u][0u] * crosses[1u][1u];
            RI -= crosses[0u][1u] * crosses[1u][0u];
            RI /= sqrt(norm);
        }
    }

    assert(RI >= -1.0);
    assert(RI <= 1.0);

    // Network scan
    for (size_t e = 0u; e < a.getNetworkSize(); ++e) {

        // Reset the edge information
        networkscan[e].corgen = 0.0;
        networkscan[e].corbreed = 0.0;
        networkscan[e].corfreq = 0.0;
        networkscan[e].avgi = 0.0;
        networkscan[e].avgj = 0.0;

        // Compute the statistics

        const size_t i = networkscan[e].loc1;
        const size_t j = networkscan[e].loc2;

        double covgen = 0.0;
        std::vector<std::vector<size_t> > ggcounts = utl::uzeros(3u, 3u);

        // Loop through individuals
        for (size_t k = 0u; k < m.getSize(); ++k) {

            const double geni = m.population[k].getLocusValue(i);
            const double genj = m.population[k].getLocusValue(j);
            const size_t zygi = m.population[k].getZygosity(i);
            const size_t zygj = m.population[k].getZygosity(j);

            covgen += geni * genj;
            ++ggcounts[zygi][zygj];

        }

        // Correlation in genetic value between partners

        covgen /= ecounts[tot];
        const double meangi = genomescan[i].meang;
        const double meangj = genomescan[j].meang;
        covgen -= meangi * meangj;

        const double vargeni = genomescan[i].varG[tot];
        const double vargenj = genomescan[j].varG[tot];

        norm = sqrt(vargeni * vargenj);
        const double corgen = norm == 0.0 ? 0.0 : covgen / norm;

        networkscan[e].corgen = corgen;

        // Correlation in breeding values between partners

        const double varbreedi = genomescan[i].varA[tot];
        const double varbreedj = genomescan[j].varA[tot];

        const std::vector<double> breedi = genomescan[i].beta;
        const std::vector<double> breedj = genomescan[j].beta;

        double covbreed = 0.0;
        for (size_t zygi : { aa, Aa, AA }) {
            for (size_t zygj : {bb, Bb, BB}) {
                covbreed += ggcounts[zygi][zygj] * breedi[zygi] * breedj[zygj];
            }
        }
        covbreed /= ecounts[tot];

        norm = sqrt(varbreedi * varbreedj);
        const double corbreed = norm == 0.0 ? 0.0 : covbreed / norm;

        networkscan[e].corbreed = corbreed;

        // Correlation in allele counts between partners

        const double varzygi = genomescan[i].varZ;
        const double varzygj = genomescan[j].varZ;

        double covzyg = 4.0 * ggcounts[AA][BB];
        covzyg += 2.0 * ggcounts[Aa][BB];
        covzyg += 2.0 * ggcounts[AA][Bb];
        covzyg += ggcounts[Aa][Bb]; // aa and bb do not count
        covzyg /= ecounts[tot];
        covzyg -= 4.0 * genomescan[i].freq * genomescan[j].freq;

        norm = sqrt(varzygi * varzygj);
        const double corfreq = norm == 0.0 ? 0.0 : covzyg / norm;

        networkscan[e].corfreq = corfreq;

        // Variance in average effect due to epistasis

        const double varexpi = genomescan[i].varX;
        const double varexpj = genomescan[j].varX;

        const size_t id = networkscan[e].id;
        const size_t trait = networkscan[e].trait;

        const double weight = a.networks[trait].weights[id];

        const double scaleA = p.scaleA[trait];
        const double scaleI = p.scaleI[trait];

        const double additi = scaleA * a.effects[i];
        const double additj = scaleA * a.effects[j];

        const double ratioi = additi == 0.0 ? 0.0 : scaleI * weight / additi;
        const double ratioj = additj == 0.0 ? 0.0 : scaleI * weight / additj;

        const double alphai = genomescan[i].alpha;
        const double alphaj = genomescan[j].alpha;

        const double avgi = utl::sqr(utl::sqr(alphaj * ratioj) * varexpi);
        const double avgj = utl::sqr(utl::sqr(alphai * ratioi) * varexpj);

        networkscan[e].avgi = avgi; // for partner i
        networkscan[e].avgj = avgj; // for partner j

    }
}

std::vector<double> Collector::get_Fst() const
{
    std::vector<double> output(genomescan.size());
    for(size_t i = 0; i < genomescan.size(); ++i) {
        output[i] = genomescan[i].Fst;
    }
    return output;
}

std::vector<double> Collector::get_Gst() const
{
    std::vector<double> output(genomescan.size());
    for(size_t i = 0; i < genomescan.size(); ++i) {
        output[i] = genomescan[i].Gst;
    }
    return output;
}

std::vector<double> Collector::get_eco_trait(const MetaPop &m) const
{
   std::vector<double> output(m.population.size());
   for(size_t i = 0; i < m.population.size(); ++i) {
       output[i] = m.population[i].getTraitValue(0u);
   }
   return output;
}

std::vector<double> Collector::get_eco_trait_deme(const MetaPop &m,
                                       size_t deme) const
{
   std::vector<double> output;
   for(size_t i = 0; i < m.population.size(); ++i) {

       if(m.population[i].getHabitat() == deme) {
           output.push_back(m.population[i].getTraitValue(0u));
       }
   }
   return output;
}

std::vector<double> Collector::get_sex_trait(const MetaPop &m) const
{
   std::vector<double> output(m.population.size());
   for(size_t i = 0; i < m.population.size(); ++i) {
       output[i] = m.population[i].getTraitValue(1u);
   }
   return output;
}

std::vector<double> Collector::get_sex_trait_deme(const MetaPop &m,
                                       size_t deme) const
{
   std::vector<double> output;
   for(size_t i = 0; i < m.population.size(); ++i) {

       if(m.population[i].getHabitat() == deme) {
           output.push_back(m.population[i].getTraitValue(1u));
       }
   }
   return output;
}



std::vector<double> Collector::get_neu_trait(const MetaPop &m) const
{
   std::vector<double> output(m.population.size());
   for(size_t i = 0; i < m.population.size(); ++i) {
       output[i] = m.population[i].getTraitValue(2u);
   }
   return output;
}

std::vector<double> Collector::get_neu_trait_deme(const MetaPop &m,
                                       size_t deme) const
{
   std::vector<double> output;
   for(size_t i = 0; i < m.population.size(); ++i) {

       if(m.population[i].getHabitat() == deme) {
           output.push_back(m.population[i].getTraitValue(2u));
       }
   }
   return output;
}
