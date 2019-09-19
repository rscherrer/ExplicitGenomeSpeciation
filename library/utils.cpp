#include "utils.h"
#include <numeric>
#include <algorithm>
#include <string>
#include <iostream>

double sqr(const double &number)
{
    return number * number;
}


vecDbl ones(const size_t &n)
{
    vecDbl ones;
    for (size_t i = 0u; i < n; ++i)
        ones.push_back(1.0);
    return ones;
}


vecDbl zeros(const size_t &n)
{
    vecDbl zeros;
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0.0);
    return zeros;
}


vecUns uzeros(const size_t &n)
{
    vecUns zeros;
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0u);
    return zeros;
}


vecBool falses(const size_t &n)
{
    vecBool falses;
    for (size_t i = 0u; i < n; ++i)
        falses.push_back(false);
    return falses;
}

vecDbl rep(const double &x, const size_t &n)
{
    vecDbl reps;
    for (size_t i = 0u; i < n; ++i)
        reps.push_back(x);
    return reps;
}


double sum(vecDbl &v)
{
    return std::accumulate(v.begin(), v.end(), 0);
}


size_t argmin(vecDbl &v)
{
    return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}

size_t sumbool(vecBool v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i)
        sum += v[i];
    return sum;
}

size_t sumu(vecUns &v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i)
        sum += v[i];
    return sum;
}

double size2dbl(const size_t &x)
{
    return static_cast<double>(x);
}
