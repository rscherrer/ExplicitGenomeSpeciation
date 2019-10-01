#include "Stats.h"

size_t Stats::Locus::totcount = 0u;
vecUns Stats::Locus::ecocounts = { 0u, 0u };



//-///////////////////////////////////
// Module for single-locus statistics
//-///////////////////////////////////


// Perform statistical analyses

void Stats::Locus::accumulate(const vecPop& pops)
{
    for (auto pop : pops) {
        for (size_t ind = 0u; ind < pop->getPopSize(); ++ind) {

            const size_t eco = pop->getEcotype(ind);
            const size_t gen = pop->getLocusValue(ind, id);
            const size_t zyg = pop->getZygosity(ind, id);

            ++ecostats[eco].genocounts[zyg];
            ecostats[eco].genosumgens[zyg] += gen;
            ecostats[eco].genossqgens[zyg] += utl::sqr(gen);

        }
    }

    for (size_t eco : { 0u, 1u }) {

        ecostats[eco].setup();

        for (size_t zyg : { aa, Aa, AA }) {

            genocounts[zyg] += ecostats[eco].genocounts[zyg];
            genosumgens[zyg] += ecostats[eco].genosumgens[zyg];
            genossqgens[zyg] += ecostats[eco].genossqgens[zyg];

            assert(ecostats[eco].genossqgens[zyg] >= 0.0);
            assert(genossqgens[zyg] >= 0.0);
        }

        sumgen += ecostats[eco].sumgen;
        ssqgen += ecostats[eco].ssqgen;
        assert(ssqgen >= 0.0);
    }
}

void Stats::Locus::regress()
{
    // Regression

    p = (genocounts[AA] + 0.5 * genocounts[Aa]) / totcount;
    assert(p >= 0.0);
    assert(p <= 1.0);

    meanq = 2.0 * p;
    assert(meanq >= 0.0);
    assert(meanq <= 2.0);

    varq = (4.0 * genocounts[AA] + genocounts[Aa]) / totcount - utl::sqr(meanq);
    assert(varq >= 0.0);

    covqg = (2.0 * genosumgens[AA] + genosumgens[Aa]) / totcount;
    covqg -= meanq * sumgen / totcount;

    if (varq != 0.0) alpha = covqg / varq;

    for (size_t zyg : { aa, Aa, AA }) {
        if (genocounts[zyg] > 0u) { // if there is someone with that genotype
            genobetas[zyg] = alpha * (zyg - meanq);
            genoexpecs[zyg] = sumgen / totcount - genobetas[zyg];
            genodeltas[zyg] = genosumgens[zyg] / genocounts[zyg];
            genodeltas[zyg] -= genoexpecs[zyg];
        }
    }
}

void Stats::Locus::setEcotypes()
{
    // Analyze ecotypes
    for (size_t eco : { 0u, 1u }) {
        if (ecocounts[eco] > 0u) {
            ecostats[eco].setVarG();
            ecostats[eco].setVarP(varE);
            ecostats[eco].setVarA(genobetas);
            ecostats[eco].setVarN(genocounts, genosumgens);
        }
    }
}

void Stats::Locus::setVarG()
{
    varG = ssqgen / totcount - utl::sqr(sumgen / totcount);
    assert(varG >= 0.0);
}

void Stats::Locus::setVarP()
{
    varP = varG + varE;
    assert(varP >= 0.0);
}

void Stats::Locus::setVarA()
{
    varA = utl::sqr(alpha) * varq;
    assert(varA >= 0.0);
}

void Stats::Locus::setVarD()
{
    varD += genocounts[aa] * utl::sqr(genodeltas[aa]);
    varD += genocounts[Aa] * utl::sqr(genodeltas[Aa]);
    varD += genocounts[AA] * utl::sqr(genodeltas[AA]);
    varD /= totcount;
    assert(varD >= 0.0);
}

void Stats::Locus::setVarI()
{
    if (genocounts[aa])
        varI += genossqgens[aa] - utl::sqr(genosumgens[aa]) / genocounts[aa];

    if (genocounts[Aa])
        varI += genossqgens[Aa] - utl::sqr(genosumgens[Aa]) / genocounts[Aa];

    if (genocounts[AA])
        varI += genossqgens[AA] - utl::sqr(genosumgens[AA]) / genocounts[AA];

    varI /= totcount;
    assert(varI >= 0.0);
}

void Stats::Locus::setVarN()
{
    varN = varD + 0.5 * varI;
    assert(varN >= 0.0);
}

void Stats::Locus::setVarS()
{
    varS += ecocounts[0u] * utl::sqr(ecostats[0u].p);
    varS += ecocounts[1u] * utl::sqr(ecostats[1u].p);
    varS /= totcount;
    varS -= utl::sqr(p); // I think varS can be negative
}

void Stats::Locus::setVarT()
{
    varT += p * (1.0 - p);
    assert(varT >= 0.0);
}

void Stats::Locus::setPst()
{
    const double v0 = getVarP(0u);
    const double v1 = getVarP(1u);
    const double v = varP;
    const size_t n0 = ecocounts[0u];
    const size_t n1 = ecocounts[1u];
    const size_t n = totcount;

    Pst = utl::Xst(v0, v1, v, n0, n1, n);
    assert(Pst >= 0.0);
    assert(Pst <= 1.0);
}

void Stats::Locus::setGst()
{
    const double v0 = getVarG(0u);
    const double v1 = getVarG(1u);
    const double v = varG;
    const size_t n0 = ecocounts[0u];
    const size_t n1 = ecocounts[1u];
    const size_t n = totcount;

    Gst = utl::Xst(v0, v1, v, n0, n1, n);
    assert(Gst >= 0.0);
    assert(Gst <= 1.0);
}

void Stats::Locus::setQst()
{
    const double v0 = getVarA(0u);
    const double v1 = getVarA(1u);
    const double v = varA;
    const size_t n0 = ecocounts[0u];
    const size_t n1 = ecocounts[1u];
    const size_t n = totcount;

    Qst = utl::Xst(v0, v1, v, n0, n1, n);
    assert(Qst >= 0.0);
    assert(Qst <= 1.0);
}

void Stats::Locus::setCst()
{
    const double v0 = getVarN(0u);
    const double v1 = getVarN(1u);
    const double v = varN;
    const size_t n0 = ecocounts[0u];
    const size_t n1 = ecocounts[1u];
    const size_t n = totcount;

    Cst = utl::Xst(v0, v1, v, n0, n1, n);
    assert(Cst >= 0.0);
    assert(Cst <= 1.0);
}

void Stats::Locus::setFst()
{
    const double p0 = getP(0u);
    const double p1 = getP(1u);
    const size_t n0 = ecocounts[0u];
    const size_t n1 = ecocounts[1u];
    const double h = (n0 * p0 * (1 - p0) + n1 * p1 * (1 - p1)) / totcount;
    const double H = p * (1.0 - p);
    assert(h >= 0.0);
    assert(h <= 1.0);
    assert(H >= 0.0);
    assert(H <= 1.0);

    if (H == 0u) return;

    Fst = 1.0 - h / H;
    assert(Fst <= 1.0); // Fst could be negative
}





//-//////////////////////////////////////////////////
// Module for within-ecotype single-locus statistics
//-//////////////////////////////////////////////////


// Initialization

std::vector<Stats::Locus::Ecotype> Stats::Locus::makeEmptyEcotypes()
{
    std::vector<Stats::Locus::Ecotype> vec;
    for (size_t i = 0u; i < 2u; ++i)
        vec.push_back(Stats::Locus::Ecotype(i));
    return vec;
}


// Variance calculation

void Stats::Locus::Ecotype::setup()
{
    for (size_t zyg : { aa, Aa, AA }) {
        ecocounts[eco] += genocounts[zyg];
        sumgen += genosumgens[zyg];
        ssqgen += genossqgens[zyg];
    }

    assert(ssqgen >= 0.0);

    if (ecocounts[eco] == 0u) return;
    p = (genocounts[AA] + 0.5 * genocounts[Aa]) / ecocounts[eco];
    assert(p >= 0.0);
    assert(p <= 1.0);
}

void Stats::Locus::Ecotype::setVarG()
{
    varG = ssqgen / ecocounts[eco] - utl::sqr(sumgen / ecocounts[eco]);
    assert(varG >= 0.0);
}

void Stats::Locus::Ecotype::setVarP(const double &varE)
{
    varP = varG + varE;
    assert(varP >= 0.0);
}

void Stats::Locus::Ecotype::setVarA(const vecDbl &brvs)
{

    // Mean squared breeding value...
    varA = genocounts[aa] * utl::sqr(brvs[aa]);
    varA += genocounts[Aa] * utl::sqr(brvs[Aa]);
    varA += genocounts[AA] * utl::sqr(brvs[AA]);
    varA /= ecocounts[eco];

    // ... minus squared mean breeding value
    double meanb = genocounts[aa] * brvs[aa];
    meanb += genocounts[Aa] * brvs[Aa];
    meanb += genocounts[AA] * brvs[AA];
    meanb /= ecocounts[eco];
    varA -= utl::sqr(meanb);
    assert(varA >= 0.0);
}

void Stats::Locus::Ecotype::setVarN(const vecUns &ntot, const vecDbl &sumgentot)
{

    // Mean squared non-additive deviation...
    if (genocounts[aa])
        varN += sumgentot[aa] * genosumgens[aa] / genocounts[aa];
    if (genocounts[Aa])
        varN += sumgentot[Aa] * genosumgens[Aa] / genocounts[Aa];
    if (genocounts[AA])
        varN += sumgentot[AA] * genosumgens[AA] / genocounts[AA];
    varN *= -2.0;
    varN += ssqgen;
    if (ntot[aa])
        varN += genocounts[aa] / ntot[aa] * utl::sqr(sumgentot[aa]);
    if (ntot[Aa])
        varN += genocounts[Aa] / ntot[Aa] * utl::sqr(sumgentot[Aa]);
    if (ntot[AA])
        varN += genocounts[AA] / ntot[AA] * utl::sqr(sumgentot[AA]);
    varN /= ecocounts[eco];

    // ... minus squared mean non-additive deviation
    double meandev = sumgen;
    if (ntot[aa])
        meandev -= genocounts[aa] / ntot[aa] * sumgentot[aa];
    if (ntot[Aa])
        meandev -= genocounts[Aa] / ntot[Aa] * sumgentot[Aa];
    if (ntot[AA])
        meandev -= genocounts[AA] / ntot[AA] * sumgentot[AA];
    meandev /= ecocounts[eco];
    varN -= utl::sqr(meandev);
    assert(varN >= 0.0);
}
