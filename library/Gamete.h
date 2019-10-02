#ifndef EXPLICITGENOMESPECIATION_GAMETE_H
#define EXPLICITGENOMESPECIATION_GAMETE_H

#include "Random.h"
#include "Types.h"
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<> Haplotype;

class Gamete {

    friend class Deme;

public:

    Gamete(const Haplotype &hap, const double &mutrate) :
        seq(hap)
    {
        mutate(mutrate);
    }

    size_t getAlleleSum() const { return seq.count(); }
    size_t getAllele(const size_t &l) const { return seq.test(l); }

private:

    void mutate(const double&);

    Haplotype seq;

};

#endif
