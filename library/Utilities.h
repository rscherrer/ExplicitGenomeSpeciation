#ifndef EXPLICITGENOMESPECIATION_UTILITIES_H
#define EXPLICITGENOMESPECIATION_UTILITIES_H

#include "Types.h"
#include <cstddef>
#include <stddef.h>
#include <numeric>
#include <cassert>
#include <algorithm>

namespace utl
{

    double sqr(const double&);
    vecDbl ones(const size_t &n);
    vecDbl zeros(const size_t&);    
    Matrix zeros(const size_t&, const size_t&);
    Matx3d zeros(const size_t&, const size_t&, const size_t&);
    vecUns uzeros(const size_t&); // unsigned zeros
    MatUns uzeros(const size_t&, const size_t&);
    vecDbl rep(const double&, const size_t&);
    vecUns repUns(const size_t&, const size_t&);
    double sum(vecDbl&);
    double sum(Matrix&);
    size_t argmin(vecDbl&);
    size_t sumu(const vecUns&);
    void marginalize(Matrix&);
    void marginalize(MatUns&);
    Matrix dividemat(const Matrix&, const MatUns&);
    double size2dbl(const size_t&);

}

#endif
