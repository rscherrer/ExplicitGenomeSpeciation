#include "Utilities.h"
#include <numeric>
#include <algorithm>
#include <string>
#include <iostream>

double utl::sqr(const double &number)
{
    return number * number;
}


vecDbl utl::ones(const size_t &n)
{
    vecDbl ones;
    for (size_t i = 0u; i < n; ++i)
        ones.push_back(1.0);
    return ones;
}


vecDbl utl::zeros(const size_t &n)
{
    vecDbl zeros;
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0.0);
    return zeros;
}


vecUns utl::uzeros(const size_t &n)
{
    vecUns zeros;
    for (size_t i = 0u; i < n; ++i)
        zeros.push_back(0u);
    return zeros;
}


Matrix utl::matzeros(const size_t &ncol, const size_t &nrow)
{
    Matrix mat;
    for (size_t i = 0u; i < ncol; ++i)
        mat.push_back(zeros(nrow));
    return mat;
}

MatUns utl::matuzeros(const size_t &ncol, const size_t &nrow)
{
    MatUns mat;
    for (size_t i = 0u; i < ncol; ++i)
        mat.push_back(uzeros(nrow));
    return mat;
}

vecBool utl::falses(const size_t &n)
{
    vecBool falses;
    for (size_t i = 0u; i < n; ++i)
        falses.push_back(false);
    return falses;
}

vecDbl utl::rep(const double &x, const size_t &n)
{
    vecDbl reps;
    for (size_t i = 0u; i < n; ++i)
        reps.push_back(x);
    return reps;
}

vecUns utl::repUns(const size_t &x, const size_t &n)
{
    vecUns reps;
    for (size_t i = 0u; i < n; ++i)
        reps.push_back(x);
    return reps;
}


double utl::sum(vecDbl &v)
{
    return std::accumulate(v.begin(), v.end(), 0);
}


size_t utl::argmin(vecDbl &v)
{
    return std::distance(v.begin(), std::min_element(v.begin(), v.end()));
}

size_t utl::sumbool(vecBool v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i)
        sum += v[i];
    return sum;
}

size_t utl::sumu(vecUns &v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i)
        sum += v[i];
    return sum;
}

double utl::size2dbl(const size_t &x)
{
    return static_cast<double>(x);
}
