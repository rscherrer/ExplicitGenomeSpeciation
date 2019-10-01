#ifndef EXPLICITGENOMESPECIATION_RANDOM_H
#define EXPLICITGENOMESPECIATION_RANDOM_H

#include "Types.h"
#include <random>
//#include <vector>
#include <stddef.h>

namespace rnd
{

    /// Random number generator
    extern std::mt19937_64 rng;

    /// Functions to generate random numbers
    bool bernoulli(const double&);
    size_t binomial(const size_t&, const double&);
    size_t poisson(const double&);
    size_t geometric(const double&);
    size_t random(const size_t&);
    double uniform(const double&);
    double normal(const double&, const double&);
    double hnormal(const double&);
    double exponential(const double&);
    size_t sample(const vecDbl&);
    double flip(const double&, const double&);
    double gamma(const double&, const double&);
    double bigamma(const double&, const double&);

}

#endif
