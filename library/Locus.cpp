#include "Stats.h"

size_t Stats::Locus::totcount = 0u;
vecUns Stats::Locus::ecocounts = { 0u, 0u };



//-///////////////////////////////////
// Module for single-locus statistics
//-///////////////////////////////////


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
        }

        sumgen += ecostats[eco].sumgen;
        ssqgen += ecostats[eco].ssqgen;
    }
}

void Stats::Locus::regress()
{
    // Regression
    p = (genocounts[AA] + 0.5 * genocounts[Aa]) / totcount;
    meanq = 2.0 * p;
    varq = (4.0 * genocounts[AA] + genocounts[Aa]) / totcount - utl::sqr(meanq);
    covqg = (2.0 * genosumgens[AA] + genosumgens[Aa]) / totcount;
    covqg -= meanq * sumgen / totcount;
    alpha = covqg / varq;

    for (size_t zyg : { aa, Aa, AA }) {
        genobetas[zyg] = alpha * (zyg - meanq);
        genoexpecs[zyg] = sumgen / totcount - genobetas[zyg];
        genodeltas[zyg] = genosumgens[zyg] / genocounts[zyg] * genoexpecs[zyg];
    }
}

void Stats::Locus::setEcotypes()
{
    // Analyze ecotypes
    for (size_t eco : { 0u, 1u }) {
        ecostats[eco].setVarG();
        ecostats[eco].setVarP(varE);
        ecostats[eco].setVarA(genobetas);
        ecostats[eco].setVarN(genocounts, genosumgens);
    }
}

void Stats::Locus::setVarG()
{
    varG = ssqgen / totcount - utl::sqr(sumgen / totcount);
}

void Stats::Locus::setVarP()
{
    varP = varG + varE;
}

void Stats::Locus::setVarA()
{
    varA = utl::sqr(alpha) * varq;
}

void Stats::Locus::setVarD()
{
    varD += genocounts[aa] * utl::sqr(genodeltas[aa]);
    varD += genocounts[Aa] * utl::sqr(genodeltas[Aa]);
    varD += genocounts[AA] * utl::sqr(genodeltas[AA]);
    varD /= totcount;
}

void Stats::Locus::setVarI()
{
    varI += genossqgens[AA] - utl::sqr(genosumgens[AA]) / genocounts[AA];
    varI += genossqgens[Aa] - utl::sqr(genosumgens[Aa]) / genocounts[Aa];
    varI += genossqgens[aa] - utl::sqr(genosumgens[aa]) / genocounts[aa];
    varI /= totcount;
}

void Stats::Locus::setVarN()
{
    varN = varD + 0.5 * varI;
}

void Stats::Locus::setVarS()
{
    varS += ecocounts[0u] * utl::sqr(ecostats[0u].p);
    varS += ecocounts[1u] * utl::sqr(ecostats[1u].p);
    varS /= totcount;
    varS -= utl::sqr(p);
}

void Stats::Locus::setVarT()
{
    varT += p * (1.0 - p);
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
}

void Stats::Locus::setFst()
{
    const double p0 = getP(0u);
    const double p1 = getP(1u);
    const size_t n0 = ecocounts[0u];
    const size_t n1 = ecocounts[1u];
    const double h = (n0 * p0 * (1 - p0) + n1 * p1 * (1 - p1)) / totcount;
    const double H = p * (1.0 - p);
    Fst = 1.0 - h / H;
}

std::vector<Stats::Locus::Ecotype> Stats::Locus::makeEmptyEcotypes()
{
    std::vector<Stats::Locus::Ecotype> vec;
    for (size_t i = 0u; i < 2u; ++i)
        vec.push_back(Stats::Locus::Ecotype(i));
    return vec;
}

void Stats::Locus::setTotCount(const double &x)
{
    totcount = x;
}

void Stats::Locus::setEcoCounts(const vecUns &v)
{
    ecocounts = v;
}





//-//////////////////////////////////////////////////
// Module for within-ecotype single-locus statistics
//-//////////////////////////////////////////////////


void Stats::Locus::Ecotype::setup()
{
    for (size_t zyg : { aa, Aa, AA }) {
        ecocounts[eco] += genocounts[zyg];
        sumgen += genosumgens[zyg];
        ssqgen += genossqgens[zyg];
    }

    p = (genocounts[AA] + 0.5 * genocounts[Aa]) / ecocounts[eco];
}

void Stats::Locus::Ecotype::setVarG()
{
    varG = ssqgen / ecocounts[eco] - utl::sqr(sumgen / ecocounts[eco]);
}

void Stats::Locus::Ecotype::setVarP(const double &varE)
{
    varP = varG + varE;
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
}

void Stats::Locus::Ecotype::setVarN(const vecUns &ntot, const vecDbl &sumgentot)
{
    // Mean squared non-additive deviation...
    varN += sumgentot[aa] * genosumgens[aa] / genocounts[aa];
    varN += sumgentot[Aa] * genosumgens[Aa] / genocounts[Aa];
    varN += sumgentot[AA] * genosumgens[AA] / genocounts[AA];
    varN *= -2.0;
    varN += ssqgen;
    varN += genocounts[aa] / ntot[aa] * utl::sqr(sumgentot[aa]);
    varN += genocounts[Aa] / ntot[Aa] * utl::sqr(sumgentot[Aa]);
    varN += genocounts[AA] / ntot[AA] * utl::sqr(sumgentot[AA]);
    varN /= ecocounts[eco];

    // ... minus squared mean non-additive deviation
    double meandev = sumgen;
    meandev -= genocounts[aa] / ntot[aa] * sumgentot[aa];
    meandev -= genocounts[Aa] / ntot[Aa] * sumgentot[Aa];
    meandev -= genocounts[AA] / ntot[AA] * sumgentot[AA];
    meandev /= ecocounts[eco];
    varN -= utl::sqr(meandev);
}
