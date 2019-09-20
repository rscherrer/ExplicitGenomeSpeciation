#ifndef EXPLICITGENOMESPECIATION_UTILITIES_H
#define EXPLICITGENOMESPECIATION_UTILITIES_H

#include "types.h"
#include <vector>
#include <cstddef>
#include <string>
#include <stddef.h>

namespace utl
{

    double sqr(const double&);
    vecDbl ones(const size_t &n);
    vecDbl zeros(const size_t&);
    vecUns uzeros(const size_t&); // unsigned zeros
    Matrix matzeros(const size_t&, const size_t&);
    MatUns matuzeros(const size_t&, const size_t&);
    vecBool falses(const size_t&);
    vecDbl rep(const double&, const size_t&);
    vecUns repUns(const size_t&, const size_t&);
    double sum(vecDbl&);
    size_t argmin(vecDbl&);
    size_t sumbool(vecBool);
    size_t sumu(vecUns&);
    double size2dbl(const size_t&);

}

#endif
