#include "MetaPop.h"

void MetaPop::resetEcoTraits(const size_t &p, const double &x)
{
    pops[p]->resetEcoTraits(x, ecosel, maxfeed);
}

void MetaPop::resetMatePrefs(const size_t &p, const double &y)
{
    pops[p]->resetMatePrefs(y);
}

void MetaPop::resetEcotypes(const size_t &p, const size_t &e)
{
    pops[p]->resetEcotypes(e);
}

void MetaPop::resetGenders(const size_t &p, const bool &sex)
{
    pops[p]->resetGenders(sex);
    pops[p]->sortSexes();
}

double Xst(const vecDbl &v, const vecUns &n, const double &tiny = 1E-15)
{
    if (v[2u] < tiny) return 0.0;
    const double xst = 1.0 - (n[0u] * v[0u] + n[1u] * v[1u]) / (n[2u] * v[2u]);
    if (xst < tiny) return 0.0;
    if (xst > 1.0 - tiny) return 1.0;
    return xst;
}

void MetaPop::analyze(const GenArch &arch)
{

    // Reset statistics

    // Note:
    // In this function, vectors with values "per ecotype" have three elements
    // The two first are for each of the two ecotypes
    // The third is for the value across the whole metapopulation

    const double tiny = 1E-15;
    const size_t nloci = arch.nLoci;

    // Mean and variance components
    // One value per trait per ecotype
    // The within-ecotype variances are used in computing divergence statistics
    meanPhenotypes = utl::matzeros(3u, 3u);
    pheVariances = utl::matzeros(3u, 3u);
    genVariances = utl::matzeros(3u, 3u);
    addVariances = utl::matzeros(3u, 3u);
    nadVariances = utl::matzeros(3u, 3u);

    // Dominance and interaction variance
    // One value per trait
    // These are not used for computing divergence statistics
    domVariances = utl::zeros(3u);
    intVariances = utl::zeros(3u);

    // Divergence between ecotypes
    // One value per trait
    Pst = utl::zeros(3u); // phenotypic
    Gst = utl::zeros(3u); // genetic (values)
    Qst = utl::zeros(3u); // additive
    Cst = utl::zeros(3u); // non-additive
    Fst = utl::zeros(3u); // in allele frequency (true genetic divergence)

    // Genome scans of divergence between ecotypes
    // One value per locus (no matter what trait the locus is coding for)
    varPScan = utl::zeros(nloci);
    varGScan = utl::zeros(nloci);
    varAScan = utl::zeros(nloci);
    varNScan = utl::zeros(nloci);
    PstScan = utl::zeros(nloci);
    GstScan = utl::zeros(nloci);
    QstScan = utl::zeros(nloci);
    CstScan = utl::zeros(nloci);
    FstScan = utl::zeros(nloci);

    // Variances in heterozygosity
    // Used for genome-wide Fst calculations
    // One value per trait
    vecDbl varS = utl::zeros(3u);
    vecDbl varT = utl::zeros(3u);

    // Mean genetic values
    // One value per trait per ecotype
    Matrix meanGenValues = utl::matzeros(3u, 3u);

    // Reset ecotypes
    for (size_t eco = 0u; eco < 2u; ++eco) {
        ecotypes[eco].clear();
        assert(ecotypes[eco].size() == 0u);
    }

    // Metapopulation size, sum and SS phenotypes and genetic values
    size_t metapopsize = 0u;
    for (size_t p = 0u; p < 2u; ++p) {
        metapopsize += pops[p]->getPopSize();
        for (auto ind : pops[p]->individuals) {
            const vecDbl traitvalues = ind->getTraits();
            const vecDbl genvalues = ind->getGenValues();
            for (size_t trait = 0u; trait < 3u; ++trait) {
                meanPhenotypes[trait][2u] += traitvalues[trait];
                meanGenValues[trait][2u] += traitvalues[trait];
                pheVariances[trait][2u] += utl::sqr(traitvalues[trait]);
                genVariances[trait][2u] += utl::sqr(traitvalues[trait]);
            }
        }
    }

    // Mean and variance in phenotypes and genetic values
    // The mean ecological trait will serve to classify ecotypes
    for (size_t trait = 0u; trait < 3u; ++trait) {
        meanPhenotypes[trait][2u] /= metapopsize;
        meanGenValues[trait][2u] /= metapopsize;
        pheVariances[trait][2u] /= metapopsize;
        genVariances[trait][2u] /= metapopsize;
        pheVariances[trait][2u] -= utl::sqr(meanPhenotypes[trait][2u]);
        genVariances[trait][2u] -= utl::sqr(meanGenValues[trait][2u]);

        if (pheVariances[trait][2u] < tiny) pheVariances[trait][2u] = 0.0;
        if (genVariances[trait][2u] < tiny) genVariances[trait][2u] = 0.0;

        assert(pheVariances[trait][2u] >= 0.0);
        assert(genVariances[trait][2u] >= 0.0);
    }

    // Census per ecotype and whole population
    vecUns census = { 0u, 0u, metapopsize };

    // Within-ecotype sum and SS in phenotype and genetic values
    for (size_t p = 0u; p < 2u; ++p) {
        for (auto ind : pops[p]->individuals) {

            vecDbl traitvalues = ind->getTraits();
            vecDbl genvalues = ind->getGenValues();

            size_t eco = ind->getEcoTrait() < meanPhenotypes[0u][2u];
            assert(eco == 0u || eco == 1u );
            ind->setEcotype(eco);
            ecotypes[eco].push_back(ind);
            ++census[eco];

            for (size_t trait = 0u; trait < 3u; ++trait) {
                meanPhenotypes[trait][eco] += traitvalues[trait];
                meanGenValues[trait][eco] += genvalues[trait];
                pheVariances[trait][eco] += utl::sqr(traitvalues[trait]);
                genVariances[trait][eco] += utl::sqr(genvalues[trait]);
            }
        }
    }

    assert(census[0u] + census[1u] == metapopsize);

    // Test: var or mean should be zero if no-one in ecotype

    // Within-ecotype mean and variance in phenotypes and genetic values
    for (size_t trait = 0u; trait < 3u; ++trait) {
        for (size_t eco = 0u; eco < 2u; ++eco) {
            if (census[eco]) {
                meanPhenotypes[trait][eco] /= census[eco];
                meanGenValues[trait][eco] /= census[eco];
                pheVariances[trait][eco] /= census[eco];
                genVariances[trait][eco] /= census[eco];
            }

            pheVariances[trait][eco] -= utl::sqr(meanPhenotypes[trait][eco]);
            genVariances[trait][eco] -= utl::sqr(meanGenValues[trait][eco]);

            if (pheVariances[trait][eco] < tiny) pheVariances[trait][eco] = 0.0;
            if (genVariances[trait][eco] < tiny) genVariances[trait][eco] = 0.0;

            assert(pheVariances[trait][eco] >= 0.0);
            assert(genVariances[trait][eco] >= 0.0);
        }
    }

    // Locus-specific environmental variance
    vecDbl locusVarE = utl::zeros(3u);
    for (size_t trait = 0u; trait < 3u; ++trait) {
        locusVarE[trait] = utl::sqr(arch.scaleE[trait]) / nloci;
        if (locusVarE[trait] < tiny) locusVarE[trait] = 0.0;
        assert(locusVarE[trait] >= 0.0);
    }

    // Locus-specific variance decomposition
    for (size_t locus = 0u; locus < nloci; ++locus) {

        // Trait encoded
        const size_t trait = arch.traits[locus];

        // Genetic values and allele counts
        double meanAlleleCount = 0.0;
        double varAlleleCount = 0.0;
        vecDbl locusMeanG = utl::zeros(3u);
        vecDbl locusVarG = utl::zeros(3u);
        double covGenValueAlleleCount = 0.0;
        MatUns genotypeCounts = utl::matuzeros(3u, 3u);
        vecDbl meanGenotypeGenValues = utl::zeros(3u);

        for (size_t eco = 0u; eco < 2u; ++eco) {
            for (auto ind : ecotypes[eco]) {

                size_t zyg = ind->getZygosity(locus);
                double genvalue = ind->getLocusValue(locus);

                meanAlleleCount += zyg;
                varAlleleCount += utl::sqr(zyg);

                ++genotypeCounts[2u][zyg];
                ++genotypeCounts[eco][zyg];
                meanGenotypeGenValues[zyg] += genvalue;

                locusMeanG[eco] += genvalue;
                locusMeanG[2u] += genvalue;
                locusVarG[eco] += utl::sqr(genvalue);
                locusVarG[2u] += utl::sqr(genvalue);
                covGenValueAlleleCount += zyg * genvalue;

            }
        }

        // Means, variances and covariances in genetic values and allele
        // counts
        meanAlleleCount /= metapopsize;
        varAlleleCount /= metapopsize;
        varAlleleCount -= utl::sqr(meanAlleleCount);

        if (varAlleleCount < tiny) varAlleleCount = 0.0;
        assert(varAlleleCount >= 0.0);

        for (size_t eco = 0u; eco < 3u; ++eco) {
            const size_t n = census[eco];
            if (n) {
                locusMeanG[eco] /= n;
                locusVarG[eco] /= n;
            }
            locusVarG[eco] -= utl::sqr(locusMeanG[eco]);
            if (locusVarG[eco] < tiny) locusVarG[eco] = 0.0;

            assert(locusVarG[eco] >= 0.0);
        }
        covGenValueAlleleCount /= metapopsize;
        covGenValueAlleleCount -= meanAlleleCount * locusMeanG[2u];
        for (size_t zyg = 0u; zyg < 3u; ++zyg) {
            const size_t n = genotypeCounts[2u][zyg];
            if (n) meanGenotypeGenValues[zyg] /= n;
        }

        // Population-wide additive variance
        vecDbl locusVarA = utl::zeros(3u);
        double avgMutEffect;
        if (!varAlleleCount) {
            avgMutEffect = 0.0;
            locusVarA[2u] = 0.0;
        } else {
            avgMutEffect = covGenValueAlleleCount / varAlleleCount;
            locusVarA[2u] = utl::sqr(avgMutEffect) * varAlleleCount;
            if (locusVarA[2u] < tiny) locusVarA[2u] = 0.0;
        }

        assert(locusVarA[2u] >= 0.0);

        addVariances[trait][2u] += locusVarA[2u];

        // Population-wide breeding values and dominance variance
        vecDbl domDeviations = utl::zeros(3u);
        vecDbl breedingValues = utl::zeros(3u);
        vecDbl addExpectations = utl::zeros(3u);
        double locusVarD = 0.0;
        for (size_t zyg = 0u; zyg < 3u; ++zyg) {
            breedingValues[zyg] = avgMutEffect;
            breedingValues[zyg] *= (zyg - meanAlleleCount);
            addExpectations[zyg] = locusMeanG[2u];
            addExpectations[zyg] -= breedingValues[zyg];
            domDeviations[zyg] = meanGenotypeGenValues[zyg];
            domDeviations[zyg] -= addExpectations[zyg];
            double genotypeSSDeviation = genotypeCounts[2u][zyg];
            genotypeSSDeviation += utl::sqr(domDeviations[zyg]);
            locusVarD += genotypeSSDeviation;
        }
        locusVarD /= metapopsize;
        if (locusVarD < tiny) locusVarD = 0.0;
        assert(locusVarD >= 0.0);
        domVariances[trait] += locusVarD;

        // Deviations from additivity, and within-ecotype additive and
        // non-additive variance
        double locusVarI = 0.0;
        vecDbl locusMeanDev = utl::zeros(2u); // within-ecotype
        vecDbl locusVarN = utl::zeros(3u);
        for (size_t eco = 0u; eco < 2u; ++eco) {

            // Interaction deviations and within-ecotype
            // non-additive variance
            for (auto ind : ecotypes[eco]) {

                const size_t zyg = ind->getZygosity(locus);
                const double genvalue = ind->getLocusValue(locus);

                double intDeviation = genvalue;
                intDeviation -= meanGenotypeGenValues[zyg];
                locusVarI += utl::sqr(intDeviation);

                double totDeviation = genvalue - addExpectations[zyg];
                locusMeanDev[eco] += totDeviation;
                locusVarN[eco] += utl::sqr(totDeviation);

            }
            if (census[eco]) {
                locusMeanDev[eco] /= census[eco];
                locusVarN[eco] /= census[eco];
            }
            locusVarN[eco] -= utl::sqr(locusMeanDev[eco]);
            if (locusVarN[eco] < tiny) locusVarN[eco] = 0.0;

            assert(locusVarN[eco] >= 0.0);
            nadVariances[trait][eco] += locusVarN[eco];

            // Within-ecotype additive variance
            double meanBreed = 0.0;
            double meanBreedSq = 0.0;
            for (size_t zyg = 0u; zyg < 3u; ++zyg) {
                size_t n = genotypeCounts[eco][zyg];
                double brv = breedingValues[zyg];
                meanBreed += n * brv;
                meanBreedSq += n * utl::sqr(brv);
            }
            if (census[eco]) {
                meanBreed /= census[eco];
                meanBreedSq /= census[eco];
            }
            locusVarA[eco] = meanBreedSq - utl::sqr(meanBreed);
            if (locusVarA[eco] < tiny) locusVarA[eco] = 0.0;
            assert(locusVarA[eco] >= 0.0);
            addVariances[trait][eco] += locusVarA[eco];
        }

        // Population-wide interaction variance
        locusVarI /= metapopsize;
        locusVarI *= 2.0;
        if (locusVarI < tiny) locusVarI = 0.0;
        assert(locusVarI >= 0.0);
        intVariances[trait] += locusVarI;

        // Population-wide non-additive variance
        locusVarN[2u] = locusVarD + 0.5 * locusVarI;
        if (locusVarN[2u] < tiny) locusVarN[2u] = 0.0;
        assert(locusVarN[2u] >= 0.0);
        nadVariances[trait][2u] += locusVarN[2u];

        // Phenotypic variance
        vecDbl locusVarP = utl::zeros(3u);
        for (size_t eco = 0u; eco < 3u; ++eco) {
            locusVarP[eco] = locusVarG[eco] + locusVarE[trait];
            if (locusVarP[eco] < tiny) locusVarP[eco] = 0.0;
            assert(locusVarP[eco] >= 0.0);
        }

        // Genome scans
        varPScan[locus] = locusVarP[2u];
        varGScan[locus] = locusVarG[2u];
        varAScan[locus] = locusVarA[2u];
        varNScan[locus] = locusVarN[2u];

        // Divergence statistics
        PstScan[locus] = Xst(locusVarP, census);
        GstScan[locus] = Xst(locusVarG, census);
        QstScan[locus] = Xst(locusVarA, census);
        CstScan[locus] = Xst(locusVarN, census);

        assert(PstScan[locus] >= 0.0);
        assert(GstScan[locus] >= 0.0);
        assert(QstScan[locus] >= 0.0);
        assert(CstScan[locus] >= 0.0);

        assert(PstScan[locus] <= 1.0);
        assert(GstScan[locus] <= 1.0);
        assert(QstScan[locus] <= 1.0);
        assert(CstScan[locus] <= 1.0);

        // Within-ecotype heterozygosity (locus-specific and genome-wide)
        double Hwithin = 0.0;
        for (size_t eco = 0u; eco < 2u; ++eco) {

            // Within-ecotype allele frequencies
            double alleleFreq = genotypeCounts[eco][0u]; // AA?
            alleleFreq += 0.5 * genotypeCounts[eco][1u]; // Aa
            if (census[eco]) alleleFreq /= census[eco];
            if (alleleFreq > tiny) alleleFreq = 0.0;
            if (alleleFreq > 1.0 - tiny) alleleFreq = 1.0;
            assert(alleleFreq >= 0.0);
            assert(alleleFreq <= 1.0);
            varS[trait] += census[eco] * utl::sqr(alleleFreq);
            double heterozFreq = 2.0 * alleleFreq * (1.0 - alleleFreq);
            if (heterozFreq < tiny) heterozFreq = 0.0;
            if (heterozFreq > 1.0 - tiny) heterozFreq = 1.0;
            assert(heterozFreq >= 0.0);
            assert(heterozFreq <= 1.0);
            Hwithin += census[eco] * heterozFreq;
        }
        varS[trait] /= metapopsize;
        Hwithin /= metapopsize;
        if (Hwithin < tiny) Hwithin = 0.0;
        if (Hwithin > 1.0 - tiny) Hwithin = 1.0;
        assert(Hwithin >= 0.0);
        assert(Hwithin <= 1.0);

        // Population-wide heterozygosity (locus-specific and genome-wide)
        double alleleFreq = genotypeCounts[2u][0u];
        alleleFreq += 0.5 * genotypeCounts[2u][1u];
        alleleFreq /= metapopsize;
        if (alleleFreq < tiny) alleleFreq = 0.0;
        if (alleleFreq > 1.0 - tiny) alleleFreq = 1.0;
        assert(alleleFreq >= 0.0);
        assert(alleleFreq <= 1.0);
        varS[trait] -= utl::sqr(alleleFreq);
        varT[trait] += alleleFreq * (1.0 - alleleFreq);
        double Htotal = meanAlleleCount * (1.0 - 0.5 * meanAlleleCount);
        if (Htotal < tiny) Htotal = 0.0;
        if (Htotal < 1.0 - tiny) Htotal = 1.0;
        assert(Htotal >= 0.0);
        assert(Htotal <= 1.0);

        // Locus-specific genetic divergence
        FstScan[locus] = 1.0 - Hwithin / Htotal;
        if (FstScan[locus] < tiny) FstScan[locus] = 0.0;
        if (FstScan[locus] > 1.0 - tiny) FstScan[locus] = 1.0;
        assert(FstScan[locus] >= 0.0);
        assert(FstScan[locus] <= 1.0);
    }

    // Genome-wide divergence
    for (size_t trait = 0u; trait < 3u; ++trait) {

        Pst[trait] = Xst(pheVariances[trait], census);
        Gst[trait] = Xst(genVariances[trait], census);
        Qst[trait] = Xst(addVariances[trait], census);
        Cst[trait] = Xst(nadVariances[trait], census);

        assert(Pst[trait] >= 0.0);
        assert(Gst[trait] >= 0.0);
        assert(Qst[trait] >= 0.0);
        assert(Cst[trait] >= 0.0);

        assert(Pst[trait] <= 1.0);
        assert(Gst[trait] <= 1.0);
        assert(Qst[trait] <= 1.0);
        assert(Cst[trait] <= 1.0);

        if (varS[trait] < tiny) varS[trait] = 0.0;
        if (varT[trait] < tiny) varT[trait] = 0.0;
        assert(varS[trait] >= 0.0);
        assert(varT[trait] >= 0.0);
        Fst[trait] = varS[trait] / varT[trait];
        if (Fst[trait] < tiny) Fst[trait] = 0.0;
        if (Fst[trait] > 1.0 - tiny) Fst[trait] = 1.0;
        assert(Fst[trait] >= 0.0);
        assert(Fst[trait] <= 1.0);

    }

    assert(meanPhenotypes.size() == 3u);
    assert(pheVariances.size() == 3u);
    assert(genVariances.size() == 3u);
    assert(addVariances.size() == 3u);
    assert(nadVariances.size() == 3u);
    assert(domVariances.size() == 3u);
    assert(intVariances.size() == 3u);
    assert(Pst.size() == 3u);
    assert(Gst.size() == 3u);
    assert(Qst.size() == 3u);
    assert(Cst.size() == 3u);
    assert(Fst.size() == 3u);
    assert(varPScan.size() == nloci);
    assert(varGScan.size() == nloci);
    assert(varAScan.size() == nloci);
    assert(varNScan.size() == nloci);
    assert(PstScan.size() == nloci);
    assert(GstScan.size() == nloci);
    assert(QstScan.size() == nloci);
    assert(CstScan.size() == nloci);
    assert(FstScan.size() == nloci);

}

int MetaPop::evolve(const GenArch &arch)
{
    bool isBurnin = true;
    t = - tburnin;

    Output out;
    if (record) out.openAll();    

    for (; t < tmax; ++t) {

        std::clog << "t = " << t << '\n';

        // Sort out the sexes
        pops[0u]->sortSexes();
        pops[1u]->sortSexes();

        if (t > 0 && isBurnin) {
            isBurnin = false;
            for (auto pop : pops) pop->exitBurnIn();
        }

        // Analyze and record
        if (record && t % tsave == 0u && t > 0) {
            analyze(arch);
            save(out);
        }

        // Dispersal (only if not burnin)
        if (t > 0) {
            Crowd migrants1 = pops[0u]->emigrate(dispersal);
            Crowd migrants2 = pops[1u]->emigrate(dispersal);
            pops[0u]->immigrate(migrants2);
            pops[1u]->immigrate(migrants1);
        }

        size_t isExtant = 0u;

        for (auto pop : pops) {
            pop->consume();
            pop->reproduce(birth, sexsel, matingcost, ecosel, maxfeed, arch);
            isExtant += pop->survive(survival);
        }

        if (isExtant == 0u) {
            std::clog << "The population went extinct at t = " << t << '\n';
            break;
        }
    }

    if (record) out.closeAll();

    return t;
}

void MetaPop::save(Output &out)
{
    // Time
    write(utl::size2dbl(t), out.files[0u]);

    // Census
    write(utl::size2dbl(getEcotypeSize(0u)), out.files[1u]);
    write(utl::size2dbl(getEcotypeSize(1u)), out.files[2u]);
    write(utl::size2dbl(getPopSize(0u)), out.files[3u]);
    write(utl::size2dbl(getPopSize(1u)), out.files[4u]);
    write(utl::size2dbl(getNFemales(0u)), out.files[5u]);
    write(utl::size2dbl(getNFemales(1u)), out.files[6u]);

    // Resources in each habitat
    write(getResource(0u, 0u), out.files[7u]);
    write(getResource(1u, 0u), out.files[8u]);
    write(getResource(0u, 1u), out.files[9u]);
    write(getResource(1u, 1u), out.files[10u]);

    // Quantitative genetics
    for (size_t trait = 0u; trait < 3u; ++trait) {
        const size_t off = trait * 13u;
        write(meanPhenotypes[trait][0u], out.files[11u + off]); // ecotype 0
        write(meanPhenotypes[trait][1u], out.files[12u + off]); // ecotype 1
        write(meanPhenotypes[trait][2u], out.files[13u + off]); // whole pop
        write(pheVariances[trait][2u], out.files[14u + off]); // whole pop
        write(genVariances[trait][2u], out.files[15u + off]); // whole pop
        write(addVariances[trait][2u], out.files[16u + off]); // whole pop
        write(domVariances[trait], out.files[17u + off]); // whole pop
        write(intVariances[trait], out.files[18u + off]); // whole pop
        write(Pst[trait], out.files[19u + off]);
        write(Gst[trait], out.files[20u + off]);
        write(Qst[trait], out.files[21u + off]);
        write(Cst[trait], out.files[22u + off]);
        write(Fst[trait], out.files[23u + off]);
    }

    // Speciation metrics
    write(getEcoIsolation(), out.files[49u]);
    write(getSpatialIsolation(), out.files[50u]);
    write(getMatingIsolation(), out.files[51u]);

    // Genome scans
    write(varPScan, out.files[52u]);
    write(varGScan, out.files[53u]);
    write(varAScan, out.files[54u]);
    write(varNScan, out.files[55u]);
    write(PstScan, out.files[56u]);
    write(GstScan, out.files[57u]);
    write(QstScan, out.files[58u]);
    write(CstScan, out.files[59u]);
    write(FstScan, out.files[60u]);
}

void MetaPop::write(const double &x, std::ofstream * &out)
{
    out->write((char *) &x, sizeof(x));
}

void MetaPop::write(const vecDbl &vec, std::ofstream * &out)
{
    if (vec.size() > 0.0)
        for (auto x : vec)
            write(x, out);
}

double MetaPop::getEcoIsolation()
{
    // Ecological isolation is the standard deviation in ecological trait value
    // Ripa et al. use the standard deviation of the ecological trait
    // We use Pst

    return Pst[0u];

}

double MetaPop::getSpatialIsolation()
{

    const double tiny = 1E-15;

    // Ecotype-by-habitat table
    MatUns n = { utl::uzeros(2u), utl::uzeros(2u) };

    for (size_t p = 0u; p < 2u; ++p) {
        auto pop = pops[p];
        for (size_t i = 0u; i < pop->getPopSize(); ++i) {
            auto ind = pop->individuals[i];
            size_t ecotype = ind->getEcotype();
            ++n[p][ecotype];
        }
    }

    double si = n[0u][0u] * n[1u][1u] - n[0u][1u] * n[1u][0u];

    // Normalizing product
    double norm = n[0u][0u] + n[0u][1u];
    norm *= n[1u][0u] + n[1u][1u];
    norm *= n[0u][0u] + n[1u][0u];
    norm *= n[0u][1u] + n[1u][1u];
    if (norm == 0.0) return 0.0;

    si /= sqrt(norm);

    if (si < tiny) si = 0.0;
    if (si > 1.0 - tiny) si = 1.0;
    assert(si >= 0.0);
    assert(si <= 1.0);

    return si;
}


double MetaPop::getMatingIsolation()
{

    const double tiny = 1E-15;

    // Count homogamic and heterogamic crossings

    // For each of them sample a number of males to encounter, from the metapop
    // Evaluate each male by a yes or no
    // Update the table of matings accordingly by looking at ecotypes
    // RI = 0 if mating tests are not possible

    size_t nfemales = 0u;
    size_t nmales = 0u;

    for (size_t p = 0u; p < 2u; ++p) {
        pops[p]->sortSexes();
        nfemales += pops[p]->getNFemales();
        nmales += pops[p]->getNMales();
    }

    if (nfemales == 0u || nmales == 0u)
        return 0.0;

    // Table of crossings
    MatUns m = { utl::uzeros(2u), utl::uzeros(2u) };

    // Make a vector with all males of the metapop
    Crowd allMales;
    for (size_t p = 0u; p < 2u; ++p)
        for (size_t i = 0u; i < pops[p]->getNMales(); ++i)
            allMales.push_back(pops[p]->individuals[i]);

    // Loop through the females of the metapop and test their preference
    for (size_t p = 0u; p < 2u; ++p) {
        for (size_t i = 0u; i < pops[p]->getNFemales(); ++i) {
            auto fem = pops[p]->females[i];
            size_t nencounters = rnd::poisson(1.0 / matingcost);
            while (nencounters) {
                auto candidate = allMales[rnd::random(allMales.size())];
                if (fem->acceptMate(candidate->getEcoTrait(), sexsel))
                    ++m[fem->getEcotype()][candidate->getEcotype()];
                --nencounters;
            }
        }
    }

    double ri = m[0u][0u] * m[1u][1u] - m[0u][1u] * m[1u][0u];
    double norm = m[0u][0u] + m[0u][1u];
    norm *= m[1u][0u] + m[1u][1u];
    norm *= m[0u][0u] + m[1u][0u];
    norm *= m[0u][1u] + m[1u][1u];
    if (norm == 0.0) return 0.0;

    ri /= sqrt(norm);

    if (ri < tiny) ri = 0.0;
    if (ri > 1.0 - tiny) ri = 1.0;
    assert(ri >= 0.0);
    assert(ri <= 1.0);

    return ri;

}

size_t MetaPop::getNOffspring(const size_t &p) const
{
    return pops[p]->getNOffspring();
}

double MetaPop::getResource(const size_t &p, const size_t &r) const
{
    return pops[p]->getResource(r);
}
