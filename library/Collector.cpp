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
    if (v[2u] == 0.0) return 0.0;

    const double xst = 1.0 - (n[0u] * v[0u] + n[1u] * v[1u]) / (n[2u] * v[2u]);

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

        // Calculate sums across genotypes
        for (size_t eco = 0u; eco < 2u; ++eco) {
            for (size_t zyg : { aa, Aa, AA }) {
                gcounts[eco][all] += gcounts[eco][zyg];
                gsumgen[eco][all] += gsumgen[eco][zyg];
                gssqgen[eco][all] += gssqgen[eco][zyg];
            }
        }

        // Calculate sums across ecotypes
        for (size_t zyg : { aa, Aa, AA }) {
            gcounts[tot][zyg] = gcounts[0u][zyg] + gcounts[1u][zyg];
            gsumgen[tot][zyg] = gsumgen[0u][zyg] + gsumgen[1u][zyg];
            gssqgen[tot][zyg] = gssqgen[0u][zyg] + gssqgen[1u][zyg];
        }

        // Check sums of squares
        for (size_t eco = 0u; eco < 3u; ++eco) {
            for (size_t zyg : { aa, Aa, AA } ) {
                assert(gssqgen[eco][zyg] >= 0.0);
            }
        }

        // Calculate allele frequencies within and across ecotypes
        vecDbl allfreq = utl::zeros(3u);
        for (size_t eco = 0u; eco < 3u; ++eco) {
            if (ecounts[0u]) {
                allfreq[eco] = gcounts[eco][AA] + 0.5 * gcounts[eco][Aa];
                allfreq[eco] /= ecounts[eco];
            }
            assert(allfreq[eco] >= 0.0);
            assert(allfreq[eco] <= 1.0);
        }

        // Regress genetic values on allele counts
        const double meanq = 2.0 * allfreq[tot];
        const double meang = gsumgen[tot][all] / ecounts[tot];
        double varq;
        varq = (4.0 * gcounts[tot][AA] + gcounts[tot][Aa]) / ecounts[tot];
        varq -= utl::sqr(meanq);
        assert(varq >= 0.0);
        double covqg;
        covqg = (2.0 * gsumgen[tot][AA] + gsumgen[tot][Aa]) / ecounts[tot];
        covqg -= meanq * meang;
        double alpha = 0.0;
        if (varq != 0.0) alpha = covqg / varq;

        // Calculate genotype-specific statistics
        vecDbl gbeta = utl::zeros(3u); // per genotype
        vecDbl gexpec = utl::zeros(3u); // per genotype
        vecDbl gdelta = utl::zeros(3u); // per genotype
        for (size_t zyg : { aa, Aa, AA }) {
            if (gcounts[tot][zyg] == 0u) {
                gbeta[zyg] = alpha * (zyg - meanq);
                gexpec[zyg] = meang - gbeta[zyg];
                gdelta[zyg] = gsumgen[all][zyg] / gcounts[all][zyg];
                gdelta[zyg] -= gexpec[zyg];
            }
        }

        // Calculate variance components within ecotypes
        for (size_t eco = 0u; eco < 2u; ++eco) {

            // Genetic variance
            genomescan[l].varG[eco] = gssqgen[eco][all] / ecounts[eco];
            genomescan[l].varG[eco] -= utl::sqr(gsumgen[eco][all] / ecounts[eco]);
            assert(genomescan[l].varG[eco] >= 0.0);

            // Phenotypic variance
            genomescan[l].varP[eco] = genomescan[l].varG[eco] + locusvarE;
            assert(genomescan[l].varP[eco] >= 0.0);

            // Additive variance
            genomescan[l].varA[eco] = gcounts[eco][aa] * utl::sqr(gbeta[aa]);
            genomescan[l].varA[eco] += gcounts[eco][Aa] * utl::sqr(gbeta[Aa]);
            genomescan[l].varA[eco] += gcounts[eco][AA] * utl::sqr(gbeta[AA]);
            genomescan[l].varA[eco] /= ecounts[eco];
            double meanb = gcounts[eco][aa] * gbeta[aa];
            meanb += gcounts[eco][Aa] * gbeta[Aa];
            meanb += gcounts[eco][AA] * gbeta[AA];
            meanb /= ecounts[eco];
            genomescan[l].varA[eco] -= utl::sqr(meanb);
            assert(genomescan[l].varA[eco] >= 0.0);

            // Contribute to genome-wide additive variance
            varA[genomescan[l].trait][eco] += genomescan[l].varA[eco];

            // Non-additive variance
            for (size_t zyg : { aa, Aa, AA }) {
                if (gcounts[eco][zyg]) {
                    double x = gsumgen[tot][zyg] * gsumgen[eco][zyg];
                    x /= gcounts[eco][zyg];
                    genomescan[l].varN[eco] += x;
                }
            }
            genomescan[l].varN[eco] *= -2.0;
            genomescan[l].varN[eco] += gssqgen[eco][all];
            for (size_t zyg : { aa, Aa, AA }) {
                if (gcounts[tot][zyg]) {
                    double x = gcounts[eco][zyg] * utl::sqr(gsumgen[tot][zyg]);
                    x /= gcounts[tot][zyg];
                    genomescan[l].varN[eco] += x;
                }
            }
            genomescan[l].varN[eco] /= ecounts[eco];
            double meandev = gsumgen[eco][all];
            for (size_t zyg : { aa, Aa, AA }) {
                if (gcounts[tot][zyg]) {
                    double x = gcounts[eco][zyg] * gsumgen[tot][zyg];
                    x /= gcounts[tot][zyg];
                    meandev -= x;
                }
            }
            meandev /= ecounts[eco];
            genomescan[l].varN[eco] -= utl::sqr(meandev);
            assert(genomescan[l].varN[eco] >= 0.0);

            // Contribute to genome-wide non-additive variance
            varN[genomescan[l].trait][eco] += genomescan[l].varN[eco];

        }

        // Calculate variance components across ecotypes

        // Genetic variance
        genomescan[l].varG[tot] = gssqgen[tot][all] / ecounts[tot] - utl::sqr(meang);
        assert(genomescan[l].varG[tot] >= 0.0);

        // Phenotypic variance
        genomescan[l].varP[tot] = genomescan[l].varG[tot] + locusvarE;
        assert(genomescan[l].varP[tot] >= 0.0);

        // Additive variance
        genomescan[l].varA[tot] = utl::sqr(alpha) * varq;
        assert(genomescan[l].varA[tot] >= 0.0);

        // Contribute to genome-wide additive variance
        varA[genomescan[l].trait][tot] += genomescan[l].varA[tot];

        // Dominance variance
        genomescan[l].varD = gcounts[tot][aa] * utl::sqr(gdelta[aa]);
        genomescan[l].varD += gcounts[tot][Aa] * utl::sqr(gdelta[Aa]);
        genomescan[l].varD += gcounts[tot][AA] * utl::sqr(gdelta[AA]);
        genomescan[l].varD /= ecounts[tot];
        assert(genomescan[l].varD >= 0.0);

        // Contribute to genome-wide dominance variance
        varD[genomescan[l].trait] += genomescan[l].varD;

        // Interaction variance
        genomescan[l].varI = 0.0;
        for (size_t zyg : { aa, Aa, AA }) {
            if (gcounts[tot][zyg]) {
                double x = gssqgen[tot][zyg];
                x -= utl::sqr(gsumgen[tot][zyg]) / gcounts[tot][zyg];
                genomescan[l].varI += x;
            }
        }
        genomescan[l].varI /= ecounts[tot];
        assert(genomescan[l].varI >= 0.0);

        // Contribute to genome-wide interaction variance
        varI[genomescan[l].trait] += genomescan[l].varI;

        // Non-additive variance
        genomescan[l].varN[tot] = genomescan[l].varD + 0.5 * genomescan[l].varI;
        assert(genomescan[l].varN[tot] >= 0.0);

        // Contribute to genome-wide non-additive variance
        varN[genomescan[l].trait][tot] += genomescan[l].varN[tot];

        // Variance in within-ecotype heterozygosity
        double lvarS = ecounts[0u] * utl::sqr(allfreq[0u]);
        lvarS += ecounts[1u] * utl::sqr(allfreq[1u]);
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
        genomescan[l].Fst = 1.0 - h / H;

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

            // Genetic variance
            varG[trait][eco] = essqgen[trait][eco] / ecounts[eco];
            varG[trait][eco] -= utl::sqr(esumgen[trait][eco] / ecounts[eco]);
            assert(varG[trait][eco] >= 0.0);

            // Phenotypic variance
            varP[trait][eco] = essqphe[trait][eco] / ecounts[eco];
            varP[trait][eco] -= utl::sqr(esumphe[trait][eco] / ecounts[eco]);
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
    SI = counts[0u][0u] * counts[1u][1u] - counts[0u][1u] * counts[1u][0u];
    double norm = counts[0u][0u] + counts[0u][1u];
    norm *= counts[1u][0u] + counts[1u][1u];
    norm *= counts[0u][0u] + counts[1u][0u];
    norm *= counts[0u][1u] + counts[1u][1u];
    if (norm == 0.0) SI /= sqrt(norm);
    assert(SI >= 0.0);
    assert(SI <= 1.0);

    // Perform mating trials
    mating = m.matingtrials(p);

    // Mating isolation
    RI = mating[0u][0u] * mating[1u][1u] - mating[0u][1u] * mating[1u][0u];
    norm = mating[0u][0u] + mating[0u][1u];
    norm *= mating[1u][0u] + mating[1u][1u];
    norm *= mating[0u][0u] + mating[1u][0u];
    norm *= mating[0u][1u] + mating[1u][1u];
    if (norm == 0.0) RI /= sqrt(norm);
    assert(RI >= 0.0);
    assert(RI <= 1.0);
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


