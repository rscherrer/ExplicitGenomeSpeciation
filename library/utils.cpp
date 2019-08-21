#include "utils.h"
#include <numeric>
#include <algorithm>
#include <string>
#include <iostream>

double sqr(const double &number)
{
    return number * number;
}


dVector zeros(const size_t &n)
{
    dVector zeros;
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0.0);
    return zeros;
}


uVector uzeros(const size_t &n)
{
    uVector zeros;
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0u);
    return zeros;
}


double sum(dVector &v)
{
    return std::accumulate(v.begin(), v.end(), 0);
}


size_t argmin(dVector &v)
{
    return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}

size_t sumbool(Haplotype &v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i)
        sum += v[i];
    return sum;
}

size_t sumu(uVector &v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i)
        sum += v[i];
    return sum;
}
