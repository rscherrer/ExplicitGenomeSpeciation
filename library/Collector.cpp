#include "Collector.h"

// Analytical functions of the Collector module
//---------------------------------------------

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

// Reset the fields within the collector
void Collector::reset()
{
    counts = utl::uzeros(3u, 3u); // per habitat per ecotype
    means = utl::zeros(3u, 3u, 3u); // per trait per habitat per ecotype
    varG = utl::zeros(3u, 3u); // per trait per ecotype
    varP = utl::zeros(3u, 3u); // per trait per ecotype
    varA = utl::zeros(3u, 3u); // per trait per ecotype
    varN = utl::zeros(3u, 3u); // per trait per ecotype
    varD = std::vector<double>(3u, 0.0); // per trait
    varI = std::vector<double>(3u, 0.0); // per trait
    varT = std::vector<double>(3u, 0.0); // per trait
    varS = std::vector<double>(3u, 0.0); // per trait
    Pst = std::vector<double>(3u, 0.0); // per trait
    Gst = std::vector<double>(3u, 0.0); // per trait
    Qst = std::vector<double>(3u, 0.0); // per trait
    Cst = std::vector<double>(3u, 0.0); // per trait
    Fst = std::vector<double>(3u, 0.0); // per trait
    EI = 0.0;
    SI = 0.0;
    RI = 0.0;

    // per trait per habitat per ecotype
    sumgen = utl::zeros(3u, 3u, 3u);
    sumphe = utl::zeros(3u, 3u, 3u);
    ssqgen = utl::zeros(3u, 3u, 3u);
    ssqphe = utl::zeros(3u, 3u, 3u);

    // per trait per ecotype
    esumgen = utl::zeros(3u, 3u);
    esumphe = utl::zeros(3u, 3u);
    essqgen = utl::zeros(3u, 3u);
    essqphe = utl::zeros(3u, 3u);

}

void Collector::fillMatrices(const MetaPop &m) {

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

    // Check the sums add up well
    assert(counts[0u][0u] + counts[0u][1u] + counts[1u][0u] +
     counts[1u][1u] == m.getSize());

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
    ecounts = counts[2u]; // per ecotype
    assert(ecounts[tot] = m.population.size());
    assert(ecounts[tot] == ecounts[0u] + ecounts[1u]);

    for (size_t trait = 0u; trait < 3u; ++trait) {
        esumgen[trait] = sumgen[trait][2u];
        esumphe[trait] = sumphe[trait][2u];
        essqgen[trait] = ssqgen[trait][2u];
        essqphe[trait] = ssqphe[trait][2u];
    }

}

void Locus::reset() {

    varG = std::vector<double>(3u, 0.0);
    varP = std::vector<double>(3u, 0.0);
    varA = std::vector<double>(3u, 0.0);
    varN = std::vector<double>(3u, 0.0);
    freqs = std::vector<double>(3u, 0.0);
    varD = 0.0;
    varI = 0.0;
    varQ = 0.0;
    varX = 0.0;
    varS = 0.0;
    covQG = 0.0;
    Pst = 0.0;
    Gst = 0.0;
    Qst = 0.0;
    Cst = 0.0;
    Fst = 0.0;
    alpha = 0.0;
    meanQ = 0.0;
    meanG = 0.0;
    h = 0.0;
    H = 0.0;
    hobs = std::vector<double>(2u, 0.0);

    gcounts = utl::uzeros(3u, 4u);
    gsumgen = utl::zeros(3u, 4u);
    gssqgen = utl::zeros(3u, 4u);

    gbeta = std::vector<double>(3u, 0.0);
    gexpec = std::vector<double>(3u, 0.0);
    gmeans = std::vector<double>(3u, 0.0);
    gdelta = std::vector<double>(3u, 0.0);

}

void Collector::resetLocus(const size_t &l) {
    genomescan[l].reset();
}

void Locus::fillMatrices(const MetaPop &m, const size_t &nloci) {

    // Calculate sums within ecotypes and genotypes
    for (size_t i = 0u; i < m.population.size(); ++i) {

        const size_t eco = m.population[i].getEcotype();
        const double gen = m.population[i].getLocusValue(id);
        const size_t zyg = m.population[i].getZygosity(id, nloci);

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

}

void Locus::calcAlleleFreqs(const std::vector<size_t> &ecounts) {

    // Calculate allele frequencies within and across ecotypes
    for (size_t eco = 0u; eco < 3u; ++eco) {
        if (ecounts[eco]) {
            freqs[eco] = gcounts[eco][AA] + 0.5 * gcounts[eco][Aa];
            freqs[eco] /= ecounts[eco];
        }
        assert(freqs[eco] >= 0.0);
        assert(freqs[eco] <= 1.0);
    }

}

// Calculate the mean allele count
void Locus::calcMeanQ() {
    meanQ = 2.0 * freqs[tot];
    assert(meanQ >= 0.0);
    assert(meanQ <= 2.0);
}

// Calculate the mean genetic value
void Locus::calcMeanG(const size_t &n) {
    meanG = gsumgen[tot][all] / n;
}

// Calculate the variance in allele counts
void Locus::calcVarQ(const size_t &n) {

    // varq = E(q2) - E(q)2
    // using ssq(q) = 4 nAA + nAa

    varQ = (4.0 * gcounts[tot][AA] + gcounts[tot][Aa]) / n;
    varQ -= utl::sqr(meanQ);
    assert(varQ >= 0.0);

}

// Calculate the covariance between allele count and genetic value
void Locus::calcCovQG(const size_t &n) {

    // covqg = E(qg) - E(q)E(g)
    // using sum(qg) = sum_Aa(g) + 2 sum_AA(g)

    covQG = (2.0 * gsumgen[tot][AA] + gsumgen[tot][Aa]) / n;
    covQG -= meanQ * meanG;

}

// Calculate the average mutational effect
void Locus::calcAvgMutEffect() {
    alpha = varQ == 0.0 ? 0.0 : covQG / varQ;
}

// Calculate the variance in gene expression
void Locus::calcVarX(const double &scaleD, const double &dom, const size_t &n) {

    varX = gcounts[tot][Aa] * utl::sqr(scaleD * dom);
    varX += gcounts[tot][AA];
    varX += gcounts[tot][aa];
    varX /= n;
    double meanx = gcounts[tot][Aa] * scaleD * dom;
    meanx += gcounts[tot][AA];
    meanx -= gcounts[tot][aa];
    meanx /= n;
    varX -= utl::sqr(meanx);
    assert(varX >= 0.0);

}

void Locus::calcGenotypeStats() {

    // delta = deviation in genetic value due to dominance
    // expec = expected genetic value under additivity
    // beta = breeding value
    // = deviation from population mean genetic value due to additivity

    // delta(AA) = meang(AA) - expec(AA)
    // expec(AA) = meang(all) + beta(AA)
    // beta(AA) = alpha (AA - meanq)

    for (size_t zyg : { aa, Aa, AA }) {
        if (gcounts[tot][zyg]) {
            gbeta[zyg] = alpha * (zyg - meanQ);
            gexpec[zyg] = meanG + gbeta[zyg];
            gmeans[zyg] = gsumgen[tot][zyg] / gcounts[tot][zyg];
            gdelta[zyg] = gmeans[zyg] - gexpec[zyg];
        }
    }

}

// Calculate the genetic variance
void Locus::calcVarG(const size_t &eco, const size_t &n) {
    // Genetic variance = V(g)
    varG[eco] = gssqgen[eco][all] / n;
    const double emean = gsumgen[eco][all] / n;
    varG[eco] -= utl::sqr(emean);
    utl::correct(varG[eco], 0.0, 1.0E-5);
    assert(varG[eco] >= 0.0);
}

// Calculate phenotypic variance from genetic variance
void Locus::calcVarP(const size_t &eco, const double &locusvarE) {

    varP[eco] = varG[eco] + locusvarE;
    utl::correct(varP[eco], 0.0, 1.0E-5);
    assert(varP[eco] >= 0.0);

}

// Calculate the additive variance
void Locus::calcVarA(const size_t &eco, const size_t &n) {

    // Simpler equation for the whole pop
    // Stopped using it because makes bugs harder to detect
    // (asserts may fail just because of this assumption)
    // genomescan[l].varA[tot] = utl::sqr(alpha) * varq;

    // Additive variance = V(beta)
    // Variance in breeding values
    varA[eco] = gcounts[eco][aa] * utl::sqr(gbeta[aa]);
    varA[eco] += gcounts[eco][Aa] * utl::sqr(gbeta[Aa]);
    varA[eco] += gcounts[eco][AA] * utl::sqr(gbeta[AA]);
    varA[eco] /= n;
    double meanb = gcounts[eco][aa] * gbeta[aa];
    meanb += gcounts[eco][Aa] * gbeta[Aa];
    meanb += gcounts[eco][AA] * gbeta[AA];
    meanb /= n;
    varA[eco] -= utl::sqr(meanb);
    utl::correct(varA[eco], 0.0, 1.0E-5);
    assert(varA[eco] >= 0.0);

}

// Calculate the non-additive variance within an ecotype
void Locus::calcVarN(const size_t &eco, const size_t &n) {

    // Non-additive variance = V(gamma)
    varN[eco] = gexpec[aa] * gsumgen[eco][aa];
    varN[eco] += gexpec[Aa] * gsumgen[eco][Aa];
    varN[eco] += gexpec[AA] * gsumgen[eco][AA];
    varN[eco] *= -2.0;
    varN[eco] += gcounts[eco][aa] * utl::sqr(gexpec[aa]);
    varN[eco] += gcounts[eco][Aa] * utl::sqr(gexpec[Aa]);
    varN[eco] += gcounts[eco][AA] * utl::sqr(gexpec[AA]);
    varN[eco] += gssqgen[eco][all];
    varN[eco] /= n;
    double meandev = gsumgen[eco][all];
    meandev -= gcounts[eco][aa] * gexpec[aa];
    meandev -= gcounts[eco][Aa] * gexpec[Aa];
    meandev -= gcounts[eco][AA] * gexpec[AA];
    meandev /= n;
    varN[eco] -= utl::sqr(meandev);
    utl::correct(varN[eco], 0.0, 1.0E-5);
    assert(varN[eco] >= 0.0);

}

// Calculate the dominance variance in the whole population
void Locus::calcVarD(const size_t &n) {

    // Dominance variance (across ecotypes only) = V(delta)
    varD = gcounts[tot][aa] * utl::sqr(gdelta[aa]);
    varD += gcounts[tot][Aa] * utl::sqr(gdelta[Aa]);
    varD += gcounts[tot][AA] * utl::sqr(gdelta[AA]);
    varD /= n;
    utl::correct(varD, 0.0, 1.0E-5);
    assert(varD >= 0.0);

}

// Calculate the interaction variance across the whole population
void Locus::calcVarI(const size_t &n) {

    // Interaction variance (across ecotypes only) = V(epsilon)
    varI = gmeans[aa] * gsumgen[tot][aa];
    varI += gmeans[Aa] * gsumgen[tot][Aa];
    varI += gmeans[AA] * gsumgen[tot][AA];
    varI *= -2.0;
    varI += gcounts[tot][aa] * utl::sqr(gmeans[aa]);
    varI += gcounts[tot][Aa] * utl::sqr(gmeans[Aa]);
    varI += gcounts[tot][AA] * utl::sqr(gmeans[AA]);
    varI += gssqgen[tot][all];
    varI /= n;
    utl::correct(varI, 0.0, 1.0E-5);
    assert(varI >= 0.0);

}

// Calculate within-ecotype variance in heterozygozity
void Locus::calcVarS(const size_t &n0, const size_t &n1, const size_t &n) {

    // Variance in within-ecotype heterozygosity
    varS = n0 * utl::sqr(freqs[0u]);
    varS += n1 * utl::sqr(freqs[1u]);
    varS /= n;
    varS -= utl::sqr(freqs[tot]);

}

// Calculate average heterozygosity within ecotypes
void Locus::calcHWithin(const size_t &n0, const size_t &n1, const size_t &n) {

    // Average heterozygosity within ecotypes
    h = n0 * freqs[0u] * (1.0 - freqs[0u]);
    h += n1 * freqs[1u] * (1.0 - freqs[1u]);
    h /= n;
    assert(h >= 0.0);
    assert(h <= 1.0);

}

// Calculate heterozygosity across ecotypes
void Locus::calcHAcross() {

    // Heterozygosity across ecotypes
    H = freqs[tot] * (1.0 - freqs[tot]);
    assert(H >= 0.0);
    assert(H <= 1.0);

}

// Measure observed heterozygosity within each ecotype
void Locus::calcHObserved(const size_t &eco) {

    const size_t n = gcounts[eco][3u];

    // Observed heterozygosity is the ecotype frequency of heterozygotes
    hobs[eco] = n > 0u ? gcounts[eco][1u] / n : 0.0;

    assert(hobs[eco] >= 0.0);
    assert(hobs[eco] <= 1.0);

}

// Calculate F-statistics
void Locus::partitionVariance(const std::vector<size_t> &ecounts) {

    // Calculate divergence metrics
    Pst = Xst(varP, ecounts);
    Gst = Xst(varG, ecounts);
    Qst = Xst(varA, ecounts);
    Cst = Xst(varN, ecounts);
    Fst = H > 0.0 ? 1.0 - h / H : 0.0;

    // Check divergence metrics
    assert(Pst <= 1.0);
    assert(Gst <= 1.0);
    assert(Qst <= 1.0);
    assert(Cst <= 1.0);
    assert(Fst <= 1.0);
    assert(Pst >= 0.0);
    assert(Gst >= 0.0);
    assert(Qst >= 0.0);
    assert(Cst >= 0.0);

}

// Calculate the genetic variance of a trait within an ecotype
void Collector::calcVarG(const size_t &trait, const size_t &eco) {

    // Genetic variance
    varG[trait][eco] = essqgen[trait][eco] / ecounts[eco];
    varG[trait][eco] -= utl::sqr(esumgen[trait][eco] / ecounts[eco]);
    utl::correct(varG[trait][eco], 0.0, 1.0E-5);
    assert(varG[trait][eco] >= 0.0);

}

// Calculate the phenotypic variance of a trait within an ecotype
void Collector::calcVarP(const size_t &trait, const size_t &eco) {

    // Phenotypic variance
    varP[trait][eco] = essqphe[trait][eco] / ecounts[eco];
    varP[trait][eco] -= utl::sqr(esumphe[trait][eco] / ecounts[eco]);
    utl::correct(varP[trait][eco], 0.0, 1.0E-5);
    assert(varP[trait][eco] >= 0.0);

}


// Calculate genome-wide F-statistics for a given trait
void Collector::partitionVariance(const size_t &trait) {

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

// Calculate ecological isolation (or divergence)
void Collector::calcEI() {

    EI = Pst[0u];

}

// Calculate spatial isolation
void Collector::calcSI() {

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

}

// Calculate reproductive isolation using mating trials
void Collector::calcRI(const MetaPop &m, const Param &p) {

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
        double norm = crosses[0u][0u] + crosses[0u][1u];
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

}

// Reset the statistics of a given edge in the network
void Connexion::reset() {

    // Reset the edge information
    ggcounts = utl::uzeros(3u, 3u);
    corgen = 0.0;
    corbreed = 0.0;
    corfreq = 0.0;
    avgi = 0.0;
    avgj = 0.0;

}

// Fill relevant matrices by looping through the population
void Connexion::fillMatrices(const MetaPop &m, const size_t &nloci) {

    // Loop through individuals
    for (size_t k = 0u; k < m.getSize(); ++k) {

        const double geni = m.population[k].getLocusValue(i);
        const double genj = m.population[k].getLocusValue(j);
        const size_t zygi = m.population[k].getZygosity(i, nloci);
        const size_t zygj = m.population[k].getZygosity(j, nloci);

        sprgen += geni * genj;
        ++ggcounts[zygi][zygj];

    }

}

// Calculate the correlation in genetic values between two loci
void Connexion::calcCorGen(const size_t &n,
 const std::vector<Locus> &genomescan) {

    double covgen = sprgen / n;
    covgen -= genomescan[i].meanG * genomescan[j].meanG;
    const double norm = sqrt(genomescan[i].varG[tot] * genomescan[j].varG[tot]);
    corgen = norm == 0.0 ? 0.0 : covgen / norm;

}

// Calculate the correlation in breeding values between two loci
void Connexion::calcCorBreed(const size_t &n,
 const std::vector<Locus> &genomescan) {

    double covbreed = 0.0;
    for (size_t zygi : { aa, Aa, AA }) {
        for (size_t zygj : { bb, Bb, BB }) {
            covbreed += ggcounts[zygi][zygj] * genomescan[i].gbeta[zygi] *
             genomescan[j].gbeta[zygj];
        }
    }
    covbreed /= n;

    const double norm = sqrt(genomescan[i].varA[tot] * genomescan[j].varA[tot]);
    corbreed = norm == 0.0 ? 0.0 : covbreed / norm;

}

// Calculate the correlation in allele frequencies between two loci
void Connexion::calcCorFreq(const size_t &n,
 const std::vector<Locus> &genomescan) {

    double covzyg = 4.0 * ggcounts[AA][BB];
    covzyg += 2.0 * ggcounts[Aa][BB];
    covzyg += 2.0 * ggcounts[AA][Bb];
    covzyg += ggcounts[Aa][Bb]; // aa and bb do not count
    covzyg /= n;
    covzyg -= 4.0 * genomescan[i].freqs[tot] * genomescan[j].freqs[tot];

    const double norm = sqrt(genomescan[i].varQ * genomescan[j].varQ);
    corfreq = norm == 0.0 ? 0.0 : covzyg / norm;

}

// Calculate the expected effect of genetic variation at locus i
// on the variation in the additive effect of allele substitutions at locus j
// This is mostly for plotting purposes, to detect genes that are
// expected to modify the additive effects of their interacting partners
double Connexion::calcBkgdEffect(const Locus &locusi, const Locus &locusj,
 const double &effectj, const double &weightij, const double &scaleA,
  const double &scaleI) const {

    const double alphaj = locusj.alpha;
    const double additj = scaleA * effectj;
    const double ratioj = additj == 0.0 ? 0.0 : scaleI * weightij / additj;
    const double varexpi = locusi.varX;
    return sqrt(utl::sqr(alphaj * ratioj) * varexpi);

}

// Use this calculation in a setter for both partner loci i and j
void Connexion::calcBkgdIJ(const std::vector<Locus> &genomescan,
 const GenArch &a, const Param &p) {

    avgi = calcBkgdEffect(genomescan[i], genomescan[j], a.effects[j],
     a.networks[trait].weights[id], p.scaleA[trait], p.scaleI[trait]);

    avgj = calcBkgdEffect(genomescan[j], genomescan[i], a.effects[i],
     a.networks[trait].weights[id], p.scaleA[trait], p.scaleI[trait]);

}

void Collector::analyzeLocus(const size_t &l, const MetaPop &m,
 const GenArch &a, const Param &p) {

    // Reset locus
    genomescan[l].reset();

    const double locusvarE = p.locusE[genomescan[l].trait];

    genomescan[l].fillMatrices(m, p.nloci);
    genomescan[l].calcAlleleFreqs(ecounts);

    // Regress genetic values on allele counts
    genomescan[l].calcMeanQ();
    genomescan[l].calcMeanG(ecounts[tot]);
    genomescan[l].calcVarQ(ecounts[tot]);
    genomescan[l].calcCovQG(ecounts[tot]);
    genomescan[l].calcAvgMutEffect();

    // Variance in gene expression
    genomescan[l].calcVarX(p.scaleD[genomescan[l].trait], a.dominances[l],
     ecounts[tot]);

    // Calculate genotype-specific statistics
    genomescan[l].calcGenotypeStats();

    // Calculate variance components within and across ecotypes
    for (size_t eco = 0u; eco < 3u; ++eco) {

        if (!ecounts[eco]) continue;

        genomescan[l].calcVarG(eco, ecounts[eco]);
        genomescan[l].calcVarP(eco, locusvarE);
        genomescan[l].calcVarA(eco, ecounts[eco]);
        genomescan[l].calcVarN(eco, ecounts[eco]);

        // Contribute to genome-wide statistics
        varA[genomescan[l].trait][eco] += genomescan[l].varA[eco];
        varN[genomescan[l].trait][eco] += genomescan[l].varN[eco];

        // Measure within-ecotype observed heterozygosity
        if (eco < 2u) genomescan[l].calcHObserved(eco);

    }

    // Calculate dominance and interaction variance across the population
    genomescan[l].calcVarD(ecounts[tot]);
    genomescan[l].calcVarI(ecounts[tot]);

    // Calculate expected heterozygosity
    genomescan[l].calcVarS(ecounts[0u], ecounts[1u], ecounts[tot]);
    genomescan[l].calcHWithin(ecounts[0u], ecounts[1u], ecounts[tot]);
    genomescan[l].calcHAcross();

    // Calculate F-statistics
    genomescan[l].partitionVariance(ecounts);

    // Contribute to genome-wide statistics
    varD[genomescan[l].trait] += genomescan[l].varD;
    varI[genomescan[l].trait] += genomescan[l].varI;
    varS[genomescan[l].trait] += genomescan[l].varS;
    varT[genomescan[l].trait] += genomescan[l].H;

}

void Collector::analyzeEdge(const size_t &e, const MetaPop &m, const GenArch &a,
 const Param &p) {

    networkscan[e].reset();
    networkscan[e].fillMatrices(m, p.nloci);
    networkscan[e].calcCorGen(ecounts[tot], genomescan);
    networkscan[e].calcCorBreed(ecounts[tot], genomescan);
    networkscan[e].calcCorFreq(ecounts[tot], genomescan);
    networkscan[e].calcBkgdIJ(genomescan, a, p);

}

void Collector::analyzeTrait(const size_t &trait) {

    // Within and across ecotypes
    for (size_t eco = 0u; eco < 3u; ++eco) {

        if (!ecounts[eco]) continue;

        // Calculate genetic and phenotypic variance
        calcVarG(trait, eco);
        calcVarP(trait, eco);

    }

    // Calculate F-statistics
    partitionVariance(trait);

}

// Full analysis of the simulation
void Collector::analyze(const MetaPop &m, const Param &p, const GenArch &a)
{

    // Reset
    reset();

    // Calculate sums within ecotypes and habitats
    fillMatrices(m);

    // Scan the genome and compute locus-specific variances
    for (size_t l = 0u; l < genomescan.size(); ++l)
        analyzeLocus(l, m, a, p);

    // Calculate genome-wide variances for each trait
    for (size_t trait = 0u; trait < 3u; ++trait)
        analyzeTrait(trait);

    // Speciation metrics
    calcEI();
    calcSI();
    calcRI(m, p);

    // Network scan
    for (size_t e = 0u; e < a.getNetworkSize(); ++e)
        analyzeEdge(e, m, a, p);

}

// Various public getters used in tests
//-------------------------------------

double Collector::getEI() const
{
    return EI;
}
double Collector::getSI() const
{
    return SI;
}
double Collector::getRI() const
{
    return RI;
}
double Collector::getVarP(const size_t &t) const
{
    return varP[t][2u];
}
double Collector::getVarN(const size_t &t) const
{
    return varN[t][2u];
}
double Collector::getVarI(const size_t &t) const
{
    return varI[t];
}

// Variance in locus-specific genetic value within a given genotype z
double Collector::calcLocusGenotypeVarG(const size_t &l, const size_t &z,
 const MetaPop &m, const size_t &nloci) const
{
    double ssq = 0.0;
    double sum = 0.0;
    size_t n = 0u;
    for (size_t i = 0u; i < m.getSize(); ++i) {
        const size_t zyg = m.population[i].getZygosity(l, nloci);
        if (zyg != z) continue;
        const double gen = m.population[i].getLocusValue(l);
        ssq += utl::sqr(gen);
        sum += gen;
        ++n;
    }
    return n > 0u ? ssq / n - utl::sqr(sum / n) : 0.0;
}


// Functions by Thijs for the GUI
//--------------------------------

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
