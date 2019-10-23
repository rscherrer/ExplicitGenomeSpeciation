#ifndef EXPLICITGENOMESPECIATION_RANDOM_H
#define EXPLICITGENOMESPECIATION_RANDOM_H

#include "Types.h"
#include <random>
#include <stddef.h>

namespace rnd
{

    typedef std::discrete_distribution<size_t> discrete;
    typedef std::uniform_int_distribution<size_t> random;
    typedef std::exponential_distribution<double> exponential;
    typedef std::binomial_distribution<size_t> binomial;

    /// Random number generator
    extern std::mt19937_64 rng;

    /// Functions to generate random numbers
    bool bernoulli(const double&);
    size_t poisson(const double&);
    size_t geometric(const double&);
    double uniform(const double&);
    double normal(const double&, const double&);
    double hnormal(const double&);
    size_t sample(const vecDbl&);
    double flip(const double&, const double&);
    double gamma(const double&, const double&);
    double bigamma(const double&, const double&);

}

#endif
