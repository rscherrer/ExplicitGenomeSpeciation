#include "utils.h"

double sqr(const double &number)
{
    return number * number;
}


std::vector<double> zeros(const size_t &n)
{
    std::vector<double> zeros;
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0.0);
    return zeros;
}
