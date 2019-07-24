#ifndef EXPLICITGENOMESPECIATION_RANDOM_H
#define EXPLICITGENOMESPECIATION_RANDOM_H

#include <random>
#include <vector>

struct Random
{

    /// Random number generator
    std::mt19937_64 rng;

    /// Functions to generate random numbers
    bool bernoulli(const double&);
    size_t binomial(const size_t&, const double&);
    size_t poisson(const double&);
    size_t geometric(const double&);
    size_t random(const size_t&);
    double uniform(const double&);
    double normal(const double&, const double&);
    double exponential(const double&);
    size_t sample(const double[], const size_t&);
};

#endif
