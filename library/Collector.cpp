#include "Collector.h"

vecStrings Collector::whattosave() const
{
    return { "time", "ecotype0", "ecotype1", "popsize0", "popsize1",
     "nfemales0", "nfemales1", "resource00", "resource01", "resource10",
      "resource11", "mean_eco0", "mean_eco1", "mean_eco", "varP_eco",
       "varG_eco", "varA_eco", "varD_eco", "varI_eco", "varN_eco", "Pst_eco",
        "Gst_eco", "Qst_eco", "Cst_eco", "Fst_eco", "mean_mat0", "mean_mat1",
         "mean_mat", "varP_mat", "varG_mat", "varA_mat", "varD_mat", "varI_mat",
          "varN_mat", "Pst_mat", "Gst_mat", "Qst_mat", "Cst_mat", "Fst_mat",
           "mean_ntr0", "mean_ntr1", "mean_ntr", "varP_ntr", "varG_ntr",
            "varA_ntr", "varD_ntr", "varI_ntr", "varN_ntr", "Pst_ntr",
             "Gst_ntr", "Qst_ntr", "Cst_ntr", "Fst_ntr", "ecological_iso",
              "spatial_iso", "mating_iso", "varP_scan", "varG_scan",
               "varA_scan", "varD_scan", "varI_scan", "varN_scan", "Pst_scan",
                "Gst_scan", "Qst_scan", "Cst_scan", "Fst_scan" };
}

vecLoci Collector::emptyloci(const GenArch &arch) const
{
    std::vector<Locus> loci;
    loci.reserve(arch.traits.size());
    for (size_t locus = 0u; locus < arch.traits.size(); ++locus) {
        loci.push_back(Locus(locus, arch.traits[locus]));
    }
    return loci;
}

double Xst(const vecDbl &v, const vecUns &n)
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

    xst = utl::round(xst, 10u); // avoid very small negatives by imprecision

    assert(xst >= 0.0);
    assert(xst <= 1.0);    

    return xst;
}

void Collector::analyze(const MetaPop &m, const Param &p)
{

    // Indices for readability
    const size_t tot = 2u; // across ecotypes
    const size_t all = 3u; // across genotypes
    const size_t aa = 0u;
    const size_t Aa = 1u;
    const size_t AA = 2u;

    // Reset
    counts = utl::uzeros(3u, 3u); // per habitat per ecotype
    means = utl::zeros(3u, 3u, 3u); // per trait per habitat per ecotype
    varG = utl::zeros(3u, 3u); // per trait per ecotype
    varP = utl::zeros(3u, 3u); // per trait per ecotype
    varA = utl::zeros(3u, 3u); // per trait per ecotype
    varN = utl::zeros(3u, 3u); // per trait per ecotype
    varD = utl::zeros(3u); // per trait
    varI = utl::zeros(3u); // per trait
    Pst = utl::zeros(3u); // per trait
    Gst = utl::zeros(3u); // per trait
    Qst = utl::zeros(3u); // per trait
    Cst = utl::zeros(3u); // per trait
    Fst = utl::zeros(3u); // per trait
    EI = 0.0;
    SI = 0.0;
    RI = 0.0;

    // Create tables for counts
    Matx3d sumgen = utl::zeros(3u, 3u, 3u); // per trait per habitat per ecotype
    Matx3d sumphe = utl::zeros(3u, 3u, 3u); // per trait per habitat per ecotype
    Matx3d ssqgen = utl::zeros(3u, 3u, 3u); // per trait per habitat per ecotype
    Matx3d ssqphe = utl::zeros(3u, 3u, 3u); // per trait per habitat per ecotype

    // Calculate sums within ecotypes and habitats
    for (size_t i = 0u; i < m.population.size(); ++i) {

        // Count densities
        const size_t eco = m.population[i].getEcotype();
        const size_t hab = m.population[i].getHabitat();
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
    vecUns ecounts = counts[2u]; // per ecotype
    assert(ecounts[tot] = m.population.size());
    assert(ecounts[tot] == ecounts[0u] + ecounts[1u]);
    Matrix esumgen = utl::zeros(3u, 3u); // per trait per ecotype
    Matrix esumphe = utl::zeros(3u, 3u); // per trait per ecotype
    Matrix essqgen = utl::zeros(3u, 3u); // per trait per ecotype
    Matrix essqphe = utl::zeros(3u, 3u); // per trait per ecotype

    for (size_t trait = 0u; trait < 3u; ++trait) {
        esumgen[trait] = sumgen[trait][2u];
        esumphe[trait] = sumphe[trait][2u];
        essqgen[trait] = ssqgen[trait][2u];
        essqphe[trait] = ssqphe[trait][2u];
    }

    // Initialize genome-wide variances in heterozygosity
    vecDbl varS = utl::zeros(3u); // per trait
    vecDbl varT = utl::zeros(3u); // per trait

    // Scan the genome and compute locus-specific variances
    for (size_t l = 0u; l < genomescan.size(); ++l) {

        // Reset locus
        genomescan[l].varG = utl::zeros(3u); // per ecotype
        genomescan[l].varP = utl::zeros(3u); // per ecotype
        genomescan[l].varA = utl::zeros(3u); // per ecotype
        genomescan[l].varN = utl::zeros(3u); // per ecotype
        genomescan[l].varD = 0.0;
        genomescan[l].varI = 0.0;
        genomescan[l].Pst = 0.0;
        genomescan[l].Gst = 0.0;
        genomescan[l].Qst = 0.0;
        genomescan[l].Cst = 0.0;
        genomescan[l].Fst = 0.0;

        const double locusvarE = p.locusE[genomescan[l].trait];

        // Create count tables
        MatUns gcounts = utl::uzeros(3u, 4u); // per ecotype per genotype
        Matrix gsumgen = utl::zeros(3u, 4u); // per ecotype per genotype
        Matrix gssqgen = utl::zeros(3u, 4u); // per ecotype per genotype

        // Calculate sums within ecotypes and genotypes
        for (size_t i = 0u; i < m.population.size(); ++i) {

            const size_t eco = m.population[i].getEcotype();
            const double gen = m.population[i].getLocusValue(genomescan[l].id);
            const double zyg = m.population[i].getZygosity(genomescan[l].id);

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
        vecDbl allfreq = utl::zeros(3u);
        for (size_t eco = 0u; eco < 3u; ++eco) {
            if (ecounts[eco]) {
                allfreq[eco] = gcounts[eco][AA] + 0.5 * gcounts[eco][Aa];
                allfreq[eco] /= ecounts[eco];
            }
            assert(allfreq[eco] >= 0.0);
            assert(allfreq[eco] <= 1.0);
        }

        // Regress genetic values on allele counts

        // meanq = mean allele count

        const double meanq = 2.0 * allfreq[tot];
        assert(meanq >= 0.0);
        assert(meanq <= 2.0);

        // meang = mean genetic value

        const double meang = gsumgen[tot][all] / ecounts[tot];

        // varq = variance in allele count

        // varq = E(q2) - E(q)2
        // using ssq(q) = 4 nAA + nAa

        double varq;
        varq = (4.0 * gcounts[tot][AA] + gcounts[tot][Aa]) / ecounts[tot];
        varq -= utl::sqr(meanq);
        assert(varq >= 0.0);

        // covqg = covariance between allele count and genetic value

        // covqg = E(qg) - E(q)E(g)
        // using sum(qg) = sum_Aa(g) + 2 sum_AA(g)

        double covqg;
        covqg = (2.0 * gsumgen[tot][AA] + gsumgen[tot][Aa]) / ecounts[tot];
        covqg -= meanq * meang;

        // alpha = average effect of a mutation

        double alpha = 0.0;
        if (varq != 0.0) alpha = covqg / varq;

        // Calculate genotype-specific statistics

        // delta = deviation in genetic value due to dominance
        // expec = expected genetic value under additivity
        // beta = breeding value
        // = deviation from population mean genetic value due to additivity

        // delta(AA) = meang(AA) - expec(AA)
        // expec(AA) = meang(all) + beta(AA)
        // beta(AA) = alpha (AA - meanq)

        vecDbl gbeta = utl::zeros(3u); // per genotype
        vecDbl gexpec = utl::zeros(3u); // per genotype
        vecDbl gmeans = utl::zeros(3u); // per genotype
        vecDbl gdelta = utl::zeros(3u); // per genotype
        for (size_t zyg : { aa, Aa, AA }) {
            if (gcounts[tot][zyg]) {
                gbeta[zyg] = alpha * (zyg - meanq);
                gexpec[zyg] = meang + gbeta[zyg];
                gmeans[zyg] = gsumgen[tot][zyg] / gcounts[tot][zyg];
                gdelta[zyg] = gmeans[zyg] - gexpec[zyg];
            }
        }

        // Calculate variance components within and across ecotypes
        for (size_t eco = 0u; eco < 3u; ++eco) {

            if (!ecounts[eco]) continue;

            // Genetic variance = V(g)
            genomescan[l].varG[eco] = gssqgen[eco][all] / ecounts[eco];
            const double emean = gsumgen[eco][all] / ecounts[eco];
            genomescan[l].varG[eco] -= utl::sqr(emean);
            utl::correct(genomescan[l].varG[eco]);
            assert(genomescan[l].varG[eco] >= 0.0);

            // Phenotypic variance
            genomescan[l].varP[eco] = genomescan[l].varG[eco] + locusvarE;
            utl::correct(genomescan[l].varP[eco]);
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
            utl::correct(genomescan[l].varA[eco]);
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
            utl::correct(genomescan[l].varN[eco]);
            assert(genomescan[l].varN[eco] >= 0.0);

            // Contribute to genome-wide non-additive variance
            varN[genomescan[l].trait][eco] += genomescan[l].varN[eco];

        }

        // Dominance variance (across ecotypes only) = V(delta)
        genomescan[l].varD = gcounts[tot][aa] * utl::sqr(gdelta[aa]);
        genomescan[l].varD += gcounts[tot][Aa] * utl::sqr(gdelta[Aa]);
        genomescan[l].varD += gcounts[tot][AA] * utl::sqr(gdelta[AA]);
        genomescan[l].varD /= ecounts[tot];
        utl::correct(genomescan[l].varI);
        assert(genomescan[l].varD >= 0.0);

        // Contribute to genome-wide dominance variance
        varD[genomescan[l].trait] += genomescan[l].varD;

        // Interaction variance (across ecotypes only) = V(epsilon)
        genomescan[l].varI = gmeans[aa] * gsumgen[tot][aa];
        genomescan[l].varI += gmeans[Aa] * gsumgen[tot][Aa];
        genomescan[l].varI += gmeans[AA] * gsumgen[tot][AA];
        genomescan[l].varI *= 2.0;
        genomescan[l].varI += gcounts[tot][aa] * utl::sqr(gmeans[aa]);
        genomescan[l].varI += gcounts[tot][Aa] * utl::sqr(gmeans[Aa]);
        genomescan[l].varI += gcounts[tot][AA] * utl::sqr(gmeans[AA]);
        genomescan[l].varI += gssqgen[tot][all];
        genomescan[l].varI /= ecounts[tot];
        utl::correct(genomescan[l].varI);
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
            utl::correct(varG[trait][eco]);
            assert(varG[trait][eco] >= 0.0);

            // Phenotypic variance
            varP[trait][eco] = essqphe[trait][eco] / ecounts[eco];
            varP[trait][eco] -= utl::sqr(esumphe[trait][eco] / ecounts[eco]);
            utl::correct(varP[trait][eco]);
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
    if (norm == 0.0) {
        SI = 0.0; // if an ecotype or a habitat is empty
    }
    else {
        SI = counts[0u][0u] * counts[1u][1u] - counts[0u][1u] * counts[1u][0u];
        SI /= sqrt(norm);
    }
    assert(SI >= 0.0);
    assert(SI <= 1.0);

    // Mating isolation

    // Make a vector of IDs for males of both ecotypes
    std::vector<vecUns> males(2u);
    for (size_t eco = 0u; eco < 2u; ++eco)
        males[eco].reserve(m.getSize());

    // Females no matter what ecotype
    vecUns females;
    females.reserve(m.population.size());

    // Count the sexes in each ecotype
    MatUns esexes = utl::uzeros(2u, 2u); // per ecotype per sex

    for (size_t i = 0u; i < m.population.size(); ++i) {
        const size_t sex = m.population[i].getGender();
        const size_t eco = m.population[i].getEcotype();
        ++esexes[eco][sex];
        if (!sex)
            males[eco].push_back(i);
        else
            females.push_back(i);
    }

    for (size_t eco = 0u; eco < 2u; ++eco)
        males[eco].shrink_to_fit();

    females.shrink_to_fit();

    // RI = 0.0 if a sex is missing from an ecotype

    // Perform mating trials
    if (esexes[0u][0u] && esexes[0u][1u] && esexes[1u][0u] && esexes[1u][1u]) {

        // Loop through the females
        for (size_t f = 0u; f < females.size(); ++f) {

            const size_t fem = females[f];

            // Record ecotype
            const size_t eco = m.population[fem].getEcotype();
            const size_t alt = eco == 0u ? 1u : 0u;

            size_t ntrials = p.ntrials;

            // Repeat x times and perform the average
            while (ntrials) {

                // Sample a male from each ecotype
                const size_t i = males[eco][rnd::random(males[eco].size())];
                const size_t j = males[alt][rnd::random(males[alt].size())];

                assert(m.population[i].getEcotype() == eco);
                assert(m.population[j].getEcotype() == alt);

                // Record male ecological trait values
                const double xi = m.population[i].getEcoTrait();
                const double xj = m.population[j].getEcoTrait();

                // Calculate probabilities of mating with each male
                const double hom = m.population[fem].mate(xi, p); // homogamy
                const double het = m.population[fem].mate(xj, p); // heterogamy

                assert(hom >= 0.0);
                assert(het >= 0.0);
                assert(hom <= 1.0);
                assert(het <= 1.0);

                // Assortment score (from -1 to +1)
                double assort = 0.0;
                if (hom + het) assort = (hom - het) / (hom + het);

                assert(assort >= -1.0);
                assert(assort <= 1.0);

                RI += assort;

                --ntrials;
            }
       }

       RI /= p.ntrials * females.size();

       assert(RI >= -1.0);
       assert(RI <= 1.0);

    }

    // This way involves giving a choice between the two ecotypes to the female
    // Another way would be to assess homogamy and heterogamy by sampling
    // random males from the whole populations and count the number of crossings
    // The two should give similar results, except if an ecotype is rare
    // (across the whole population)
    // If an ecotype is rare, it will be rarely encountered in the second
    // algorithm and mating isolation will be high
    // In the first algorithm, mating isolation will depend only on difference
    // in trait and mating preference, not on densities
    // Good to keep in mind

}

namespace stf // save to file
{
    void write(const double &x, std::shared_ptr<std::ofstream> &out)
    {
        out->write((char *) &x, sizeof(x));
    }

    void write(const vecDbl &vec, std::shared_ptr<std::ofstream> &out)
    {
        if (vec.size() > 0.0)
            for (size_t i = 0u; i < vec.size(); ++i)
                stf::write(vec[i], out);
    }
}

void Collector::print(const size_t &t, const MetaPop &m)
{
    size_t f = 0u; // file id

    // Time
    stf::write(utl::size2dbl(t), files[0u]); ++f;

    // Census
    stf::write(utl::size2dbl(counts[2u][2u]), files[f]); ++f; // total
    stf::write(utl::size2dbl(counts[0u][2u]), files[f]); ++f; // eco 0
    stf::write(utl::size2dbl(counts[1u][2u]), files[f]); ++f; // eco 1
    stf::write(utl::size2dbl(counts[2u][0u]), files[f]); ++f; // hab 0
    stf::write(utl::size2dbl(counts[2u][1u]), files[f]); ++f; // hab 1
    stf::write(utl::size2dbl(counts[0u][0u]), files[f]); ++f; // eco 0 hab 0
    stf::write(utl::size2dbl(counts[0u][1u]), files[f]); ++f; // eco 0 hab 1
    stf::write(utl::size2dbl(counts[1u][0u]), files[f]); ++f; // eco 1 hab 0
    stf::write(utl::size2dbl(counts[1u][1u]), files[f]); ++f; // eco 1 hab 1
    stf::write(utl::size2dbl(m.sexcounts[0u][0u]), files[f]); ++f; // fem hab 0
    stf::write(utl::size2dbl(m.sexcounts[0u][1u]), files[f]); ++f; // fem hab 1

    // Resources in each habitat
    stf::write(m.resources[0u][0u], files[f]); ++f; // hab 0 res 0
    stf::write(m.resources[0u][1u], files[f]); ++f; // hab 0 res 1
    stf::write(m.resources[1u][0u], files[f]); ++f; // hab 1 res 0
    stf::write(m.resources[1u][1u], files[f]); ++f; // hab 1 res 1

    // Quantitative genetics
    for (size_t trait = 0u; trait < 3u; ++trait) {
        stf::write(means[trait][0u][0u], files[f]); ++f; // eco 0 hab 0
        stf::write(means[trait][0u][1u], files[f]); ++f; // eco 0 hab 1
        stf::write(means[trait][1u][0u], files[f]); ++f; // eco 1 hab 0
        stf::write(means[trait][1u][1u], files[f]); ++f; // eco 1 hab 1
        stf::write(means[trait][0u][2u], files[f]); ++f; // eco 0
        stf::write(means[trait][1u][2u], files[f]); ++f; // eco 1
        stf::write(means[trait][2u][0u], files[f]); ++f; // hab 0
        stf::write(means[trait][2u][1u], files[f]); ++f; // hab 1
        stf::write(varP[trait], files[f]); ++f;
        stf::write(varG[trait], files[f]); ++f;
        stf::write(varA[trait], files[f]); ++f;
        stf::write(varD[trait], files[f]); ++f;
        stf::write(varI[trait], files[f]); ++f;
        stf::write(varN[trait], files[f]); ++f;
        stf::write(Pst[trait], files[f]); ++f;
        stf::write(Gst[trait], files[f]); ++f;
        stf::write(Qst[trait], files[f]); ++f;
        stf::write(Cst[trait], files[f]); ++f;
        stf::write(Fst[trait], files[f]); ++f;
    }

    // Speciation metrics
    stf::write(EI, files[f]); ++f;
    stf::write(SI, files[f]); ++f;
    stf::write(RI, files[f]); ++f;

    // Genome scans
    for (size_t l = 0u; l < genomescan.size(); ++l) {
        size_t off = 0u;
        stf::write(genomescan[l].varP, files[f + off]); ++off;
        stf::write(genomescan[l].varG, files[f + off]); ++off;
        stf::write(genomescan[l].varA, files[f + off]); ++off;
        stf::write(genomescan[l].varD, files[f + off]); ++off;
        stf::write(genomescan[l].varI, files[f + off]); ++off;
        stf::write(genomescan[l].varN, files[f + off]); ++off;
        stf::write(genomescan[l].Pst, files[f + off]); ++off;
        stf::write(genomescan[l].Gst, files[f + off]); ++off;
        stf::write(genomescan[l].Qst, files[f + off]); ++off;
        stf::write(genomescan[l].Cst, files[f + off]); ++off;
        stf::write(genomescan[l].Fst, files[f + off]); ++off;
    }
}


