#include "Stats.h"

size_t Stats::Trait::totcount = 0u;
vecUns Stats::Trait::ecocounts = { 0u, 0u };



//-//////////////////////////////////////
// Module for trait-specific statistics
//-//////////////////////////////////////


// Variance calculations

void Stats::Trait::setVarG()
{
    varG = ssqgen / totcount - utl::sqr(sumgen / totcount);
    assert(varG >= 0.0);
}

void Stats::Trait::setVarP()
{
    varP = ssqphe / totcount - utl::sqr(sumphe / totcount);
    assert(varP >= 0.0);
}

void Stats::Trait::setVarG(const size_t &eco)
{
    ecostats[eco].setVarG(ecocounts[eco]);
}

void Stats::Trait::setVarP(const size_t &eco)
{
    ecostats[eco].setVarP(ecocounts[eco]);
}

void Stats::Trait::setPst()
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

void Stats::Trait::setGst()
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

void Stats::Trait::setQst()
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

void Stats::Trait::setCst()
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

void Stats::Trait::setFst()
{
    if (varT == 0.0) return;
    Fst = varS / varT;
    assert(Fst <= 1.0); // Fst can be negative
}



// Getters

double Stats::Trait::getMeanP() const
{
    return sumphe / totcount;
}

double Stats::Trait::getMeanP(const size_t &eco) const
{
    return ecostats[eco].sumphe / ecocounts[eco];
}

double Stats::Trait::getSumGen(const size_t &eco) const
{
    return ecostats[eco].sumgen;
}

double Stats::Trait::getSumPhe(const size_t &eco) const
{
    return ecostats[eco].sumphe;
}

double Stats::Trait::getSsqGen(const size_t &eco) const
{
    return ecostats[eco].ssqgen;
}

double Stats::Trait::getSsqPhe(const size_t &eco) const
{
    return ecostats[eco].ssqphe;
}

double Stats::Trait::getVarP(const size_t &eco) const
{
    return ecostats[eco].varP;
}

double Stats::Trait::getVarG(const size_t &eco) const
{
    return ecostats[eco].varG;
}

double Stats::Trait::getVarA(const size_t &eco) const
{
    return ecostats[eco].varA;
}

double Stats::Trait::getVarN(const size_t &eco) const
{
    return ecostats[eco].varN;
}





//-////////////////////////////////////////////////////
// Module for within-ecotype trait-specific statistics
//-////////////////////////////////////////////////////


// Initialization

std::vector<Stats::Trait::Ecotype> Stats::Trait::makeEmptyEcotypes()
{
    std::vector<Stats::Trait::Ecotype> vec;
    for (size_t i = 0u; i < 3u; ++i)
        vec.push_back(Stats::Trait::Ecotype());
    return vec;
}


// Variance calculation

void Stats::Trait::Ecotype::setVarG(const size_t &n)
{
   if (n == 0u) return;
   varG = ssqgen / n - utl::sqr(sumgen / n);
   assert(varG >= 0.0);
}

void Stats::Trait::Ecotype::setVarP(const size_t &n)
{
   if (n == 0u) return;
   varP = ssqphe / n - utl::sqr(sumphe / n);
   assert(varP >= 0.0);
}


