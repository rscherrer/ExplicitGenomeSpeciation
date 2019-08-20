#include "utils.h"
#include <numeric>
#include <algorithm>
#include <string>
#include <iostream>

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

size_t sumbool(std::vector<bool> &v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i)
        sum += v[i];
    return sum;
}

size_t sumu(Uector &v) {
    size_t sum = 0u;
    for (size_t i = 0u; i < v.size(); ++i)
        sum += v[i];
    return sum;
}

void printvec(Vector v, std::string sep)
{
    for (size_t i = 0u; i < v.size(); ++i)
        std::cout << v[i] << sep;
}

void printuec(Uector v, std::string sep)
{
    for (size_t i = 0u; i < v.size(); ++i)
        std::cout << v[i] << sep;
}

void printbvec(Haplotype v, std::string sep)
{
    for (size_t i = 0u; i < v.size(); ++i)
        std::cout << v[i] << sep;
}
