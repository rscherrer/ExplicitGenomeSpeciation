#include "Stats.h"









//-/////////////////////////////
// Main module for statistics
//-/////////////////////////////


void Stats::reset(const size_t &t, const GenArch &arch)
{
    time = t;
    genomescan = makeEmptyGenomeScan(arch);
    traitstats = makeEmptyTraits();
    resources = utl::matzeros(2u, 2u);
    popcounts = utl::uzeros(2u);
    nfemales = utl::uzeros(2u);
    ecocounts = utl::uzeros(2u);
    totcount = 0u;
    EI = 0.0;
    SI = 0.0;
    RI = 0.0;
}

void Stats::analyze(const vecPop &pops, const GenArch &arch)
{
    for (size_t p = 0u; p < 2u; ++p) {
        for (size_t r = 0u; r < 2u; ++p)
            resources[p][r] = pops[p]->getResource(r);
        popcounts[p] = pops[p]->getPopSize();
        nfemales[p] = pops[p]->getNFemales();
    }
    totcount = popcounts[0u] + popcounts[1u];
    Stats::Trait::totcount = totcount;
    Stats::Locus::totcount = totcount;
    Stats::Trait::ecocounts = ecocounts;
    Stats::Locus::ecocounts = ecocounts;
    accumulate(pops);
    scanGenome(pops, arch.traits);
    setVariance();
}

void Stats::setEcoIsolation()
{
    EI = traitstats[0u].Pst;
}

void Stats::setSpatialIsolation(const vecPop &pops)
{
    const double tiny = 1E-15;

    // Ecotype-by-habitat table
    MatUns n = { utl::uzeros(2u), utl::uzeros(2u) };

    for (size_t p = 0u; p < 2u; ++p) {
        auto pop = pops[p];
        for (size_t ind = 0u; ind < pop->getPopSize(); ++ind) {
            size_t ecotype = pop->getEcotype(ind);
            ++n[p][ecotype];
        }
    }

    SI = n[0u][0u] * n[1u][1u] - n[0u][1u] * n[1u][0u];

    // Normalizing product
    double norm = n[0u][0u] + n[0u][1u];
    norm *= n[1u][0u] + n[1u][1u];
    norm *= n[0u][0u] + n[1u][0u];
    norm *= n[0u][1u] + n[1u][1u];
    if (norm == 0.0) return;

    SI /= sqrt(norm);

    if (SI < tiny) SI = 0.0;
    if (SI > 1.0 - tiny) SI = 1.0;
    assert(SI >= 0.0);
    assert(SI <= 1.0);
}

void Stats::setMatingIsolation(const vecPop &pops, const double &matingcost,
 const double &sexsel)
{
    const double tiny = 1E-15;

    // Count homogamic and heterogamic crossings

    // For each of them sample a number of males to encounter, from the metapop
    // Evaluate each male by a yes or no
    // Update the table of matings accordingly by looking at ecotypes
    // RI = 0 if mating tests are not possible

    size_t nfem = 0u;
    size_t nmal = 0u;

    for (size_t p = 0u; p < 2u; ++p) {
        pops[p]->sortSexes();
        nfem += pops[p]->getNFemales();
        nmal += pops[p]->getNMales();
    }

    if (nfem == 0u || nmal == 0u)
        return;

    // Table of crossings
    MatUns m = { utl::uzeros(2u), utl::uzeros(2u) };

    // Make a vector with all males of the metapop
    Crowd allMales;
    for (size_t p = 0u; p < 2u; ++p)
        for (size_t i = 0u; i < pops[p]->getNMales(); ++i)
            allMales.push_back(pops[p]->getMale(i));

    // Loop through the females of the metapop and test their preference
    for (size_t p = 0u; p < 2u; ++p) {
        for (size_t i = 0u; i < pops[p]->getNFemales(); ++i) {
            auto fem = pops[p]->getFemale(i);
            size_t nencounters = rnd::poisson(1.0 / matingcost);
            while (nencounters) {
                auto candidate = allMales[rnd::random(allMales.size())];
                if (fem->acceptMate(candidate->getEcoTrait(), sexsel))
                    ++m[fem->getEcotype()][candidate->getEcotype()];
                --nencounters;
            }
        }
    }

    RI = m[0u][0u] * m[1u][1u] - m[0u][1u] * m[1u][0u];
    double norm = m[0u][0u] + m[0u][1u];
    norm *= m[1u][0u] + m[1u][1u];
    norm *= m[0u][0u] + m[1u][0u];
    norm *= m[0u][1u] + m[1u][1u];
    if (norm == 0.0) return;

    RI /= sqrt(norm);

    if (RI < tiny) RI = 0.0;
    if (RI > 1.0 - tiny) RI = 1.0;
    assert(RI >= 0.0);
    assert(RI <= 1.0);
}


void Stats::save(Output &out)
{

    size_t f = 0u; // file id

    // Time
    write(utl::size2dbl(time), out.files[0u]); ++f;

    // Census
    write(utl::size2dbl(ecocounts[0u]), out.files[f]); ++f;
    write(utl::size2dbl(ecocounts[1u]), out.files[f]); ++f;
    write(utl::size2dbl(popcounts[0u]), out.files[f]); ++f;
    write(utl::size2dbl(popcounts[1u]), out.files[f]); ++f;
    write(utl::size2dbl(nfemales[0u]), out.files[f]); ++f;
    write(utl::size2dbl(nfemales[1u]), out.files[f]); ++f;

    // Resources in each habitat
    write(resources[0u][0u], out.files[f]); ++f;
    write(resources[0u][1u], out.files[f]); ++f;
    write(resources[1u][0u], out.files[f]); ++f;
    write(resources[1u][1u], out.files[f]); ++f;

    // Quantitative genetics
    for (size_t trait = 0u; trait < 3u; ++trait) {
        write(traitstats[trait].getMeanP(0u), out.files[f]); ++f;
        write(traitstats[trait].getMeanP(1u), out.files[f]); ++f;
        write(traitstats[trait].getMeanP(), out.files[f]); ++f;
        write(traitstats[trait].varP, out.files[f]); ++f;
        write(traitstats[trait].varG, out.files[f]); ++f;
        write(traitstats[trait].varA, out.files[f]); ++f;
        write(traitstats[trait].varD, out.files[f]); ++f;
        write(traitstats[trait].varI, out.files[f]); ++f;
        write(traitstats[trait].varN, out.files[f]); ++f;
        write(traitstats[trait].Pst, out.files[f]); ++f;
        write(traitstats[trait].Gst, out.files[f]); ++f;
        write(traitstats[trait].Qst, out.files[f]); ++f;
        write(traitstats[trait].Cst, out.files[f]); ++f;
        write(traitstats[trait].Fst, out.files[f]); ++f;
    }

    // Speciation metrics
    write(EI, out.files[f]); ++f;
    write(SI, out.files[f]); ++f;
    write(RI, out.files[f]); ++f;

    // Genome scans
    for (size_t locus = 0u; locus < genomescan.size(); ++locus) {
        size_t off = 0u;
        write(genomescan[locus]->varP, out.files[f + off]); ++off;
        write(genomescan[locus]->varG, out.files[f + off]); ++off;
        write(genomescan[locus]->varA, out.files[f + off]); ++off;
        write(genomescan[locus]->varD, out.files[f + off]); ++off;
        write(genomescan[locus]->varI, out.files[f + off]); ++off;
        write(genomescan[locus]->varN, out.files[f + off]); ++off;
        write(genomescan[locus]->Pst, out.files[f + off]); ++off;
        write(genomescan[locus]->Gst, out.files[f + off]); ++off;
        write(genomescan[locus]->Qst, out.files[f + off]); ++off;
        write(genomescan[locus]->Cst, out.files[f + off]); ++off;
        write(genomescan[locus]->Fst, out.files[f + off]); ++off;
    }

}

void Stats::write(const double &x, std::ofstream * &out)
{
    out->write((char *) &x, sizeof(x));
}

void Stats::write(const vecDbl &vec, std::ofstream * &out)
{
    if (vec.size() > 0.0)
        for (auto x : vec)
            write(x, out);
}


void Stats::accumulate(const vecPop &pops)
{
    for (auto pop : pops) {
        for (size_t ind = 0u; ind < pop->getPopSize(); ++ind) {

            const size_t eco = pop->getEcotype(ind);
            ++ecocounts[eco];

            for (size_t trait : { x, y, z }) {
                const double gen = pop->getGenValue(ind, trait);
                const double phe = pop->getTraitValue(ind, trait);
                traitstats[trait].ecostats[eco].sumgen += gen;
                traitstats[trait].ecostats[eco].sumphe += phe;
                traitstats[trait].ecostats[eco].ssqgen += utl::sqr(gen);
                traitstats[trait].ecostats[eco].ssqphe += utl::sqr(phe);
            }
        }
    }

    for (size_t trait : { x, y, z }) {
        for (size_t eco : { 0u, 1u }) {
            traitstats[trait].sumgen += traitstats[trait].ecostats[eco].sumgen;
            traitstats[trait].sumphe += traitstats[trait].ecostats[eco].sumphe;
            traitstats[trait].ssqgen += traitstats[trait].ecostats[eco].ssqgen;
            traitstats[trait].ssqphe += traitstats[trait].ecostats[eco].ssqphe;
        }
    }
}

void Stats::scanGenome(const vecPop &pops, const vecUns &traits)
{

    for (size_t locus = 0u; locus < genomescan.size(); ++locus) {

        genomescan[locus]->accumulate(pops);
        genomescan[locus]->regress();
        genomescan[locus]->setEcotypes();
        genomescan[locus]->setVarG();
        genomescan[locus]->setVarP();
        genomescan[locus]->setVarA();
        genomescan[locus]->setVarD();
        genomescan[locus]->setVarI();
        genomescan[locus]->setVarN();
        genomescan[locus]->setPst();
        genomescan[locus]->setGst();
        genomescan[locus]->setQst();
        genomescan[locus]->setCst();
        genomescan[locus]->setFst();

        contributeVarA(genomescan[locus]->varA, traits[locus]);
        contributeVarD(genomescan[locus]->varD, traits[locus]);
        contributeVarI(genomescan[locus]->varI, traits[locus]);
        contributeVarN(genomescan[locus]->varN, traits[locus]);

        for (size_t e : { 0u, 1u }) {
            contributeEcoVarA(genomescan[locus]->getVarA(e), e, traits[locus]);
            contributeEcoVarN(genomescan[locus]->getVarN(e), e, traits[locus]);
        }

    }
}

void Stats::setVariance()
{
    for (size_t zyg : { x, y, z }) {

        for (size_t eco : { 0u, 1u }) {
            traitstats[zyg].ecostats[eco].setVarP(ecocounts[eco]);
            traitstats[zyg].ecostats[eco].setVarG(ecocounts[eco]);
        }

        traitstats[zyg].setVarP();
        traitstats[zyg].setVarG();
    }
}


void Stats::contributeVarA(const double &x, const size_t &trait)
{
    traitstats[trait].varA += x;
}

void Stats::contributeVarD(const double &x, const size_t &trait)
{
    traitstats[trait].varD += x;
}

void Stats::contributeVarI(const double &x, const size_t &trait)
{
    traitstats[trait].varI += x;
}

void Stats::contributeVarN(const double &x, const size_t &trait)
{
    traitstats[trait].varN += x;
}

void Stats::contributeEcoVarA(const double &x, const size_t &eco,
 const size_t &trait)
{
    traitstats[trait].ecostats[eco].varA = x;
}

void Stats::contributeEcoVarN(const double &x, const size_t &eco,
 const size_t &trait)
{
    traitstats[trait].ecostats[eco].varN = x;
}

std::vector<Stats::Trait> Stats::makeEmptyTraits()
{
    std::vector<Stats::Trait> traits(totcount);
    for (size_t i = 0u; i < 3u; ++i)
        traits[i] = Stats::Trait();
    return traits;
}

std::vector<Stats::Locus *> Stats::makeEmptyGenomeScan(const GenArch &arch)
{
    std::vector<Stats::Locus *> scan(arch.nLoci);
    for (size_t locus = 0u; locus < arch.nLoci; ++locus) {
        const size_t trait = arch.traits[locus];
        const double envar = arch.locusVarE[trait];
        scan[locus] = new Stats::Locus(locus, envar);
    }
    return scan;
}










