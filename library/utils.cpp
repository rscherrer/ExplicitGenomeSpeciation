#include "utils.h"
#include <numeric>
#include <algorithm>

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


std::vector<size_t> uzeros(const size_t &n)
{
    std::vector<size_t> zeros;
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0u);
    return zeros;
}


double sum(std::vector<double> &v)
{
    return std::accumulate(v.begin(), v.end(), 0);
}


size_t argmin(std::vector<double> &v)
{
    return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}
