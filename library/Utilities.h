#ifndef EXPLICITGENOMESPECIATION_UTILITIES_H
#define EXPLICITGENOMESPECIATION_UTILITIES_H

#include "Types.h"
#include <cstddef>
#include <stddef.h>
#include <numeric>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>

namespace utl
{

    double sqr(const double&);
    vecDbl ones(const size_t&);
    Matrix ones(const size_t&, const size_t&);
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
    size_t dbl2size(const double&);
    double round(const double&, const size_t&);
    void correct(double&, const double&, const double&);

}

namespace stf // save to file
{
    void write(const double&, std::shared_ptr<std::ofstream>&);
    void write(const vecDbl&, std::shared_ptr<std::ofstream>&);
}


#endif
