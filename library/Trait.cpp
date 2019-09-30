#include "Stats.h"


size_t Stats::Trait::totcount = 0u;
vecUns Stats::Trait::ecocounts = { 0u, 0u };


//-//////////////////////////////////////
// Module for trait-specific statistics
//-//////////////////////////////////////

double Stats::Trait::getMeanP() const
{
    return sumphe / totcount;
}

double Stats::Trait::getMeanP(const size_t &eco) const
{
    return ecostats[eco].sumphe / ecocounts[eco];
}

void Stats::Trait::setVarG()
{
    varG = ssqgen / totcount - utl::sqr(sumgen / totcount);
}

void Stats::Trait::setVarP()
{
    varP = ssqphe / totcount - utl::sqr(sumphe / totcount);
}

std::vector<Stats::Trait::Ecotype> Stats::Trait::makeEmptyEcotypes()
{
    std::vector<Stats::Trait::Ecotype> vec;
    for (size_t i = 0u; i < 3u; ++i)
        vec.push_back(Stats::Trait::Ecotype());
    return vec;
}

void Stats::Trait::setTotCount(const double &x)
{
    totcount = x;
}

void Stats::Trait::setEcoCounts(const vecUns &v)
{
    ecocounts = v;
}








//-////////////////////////////////////////////////////
// Module for within-ecotype trait-specific statistics
//-////////////////////////////////////////////////////


void Stats::Trait::Ecotype::setVarG(const size_t &n)
{
   varG = ssqgen / n - utl::sqr(sumgen / n);
}

void Stats::Trait::Ecotype::setVarP(const size_t &n)
{
   varP = ssqphe / n - utl::sqr(sumphe / n);
}


