#ifndef EXPLICITGENOMESPECIATION_UTILITIES_H
#define EXPLICITGENOMESPECIATION_UTILITIES_H

#include <vector>
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
    std::vector<std::vector<double> > ones(const size_t&, const size_t&);
    std::vector<std::vector<double> > zeros(const size_t&, const size_t&);
    std::vector<std::vector<std::vector<double> > > zeros(const size_t&,
     const size_t&, const size_t&);
    std::vector<std::vector<size_t> > uzeros(const size_t&, const size_t&);
    double sum(std::vector<double>&);
    double sum(std::vector<std::vector<double> >&);
    size_t sum(const std::vector<size_t>&);
    size_t argmin(std::vector<double>&);
    void marginalize(std::vector<std::vector<double> >&);
    void marginalize(std::vector<std::vector<size_t> >&);
    std::vector<std::vector<double> > dividemat(
     const std::vector<std::vector<double> >&,
      const std::vector<std::vector<size_t> >&);
    double size2dbl(const size_t&);
    size_t dbl2size(const double&);
    double round(const double&, const size_t&);
    void correct(double&, const double&, const double&);

}

namespace stf // save to file
{
    void write2(const size_t&, std::shared_ptr<std::ofstream>&);
    void write(const double&, std::shared_ptr<std::ofstream>&);
    void write(const std::vector<double>&, std::shared_ptr<std::ofstream>&);
}


#endif
