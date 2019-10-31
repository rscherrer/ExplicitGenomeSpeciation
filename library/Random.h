#ifndef EXPLICITGENOMESPECIATION_RANDOM_H
#define EXPLICITGENOMESPECIATION_RANDOM_H

#include "Types.h"
#include "rndutils.hpp"
#include <random>
#include <stddef.h>

namespace rnd
{

    // Probability distributions
    typedef std::discrete_distribution<size_t> discrete;
    typedef rndutils::mutable_discrete_distribution<size_t> mdiscrete;
    typedef std::uniform_int_distribution<size_t> random;
    typedef std::exponential_distribution<double> exponential;
    typedef std::binomial_distribution<size_t> binomial;
    typedef std::poisson_distribution<size_t> poisson;
    typedef std::geometric_distribution<size_t> geometric;
    typedef std::uniform_real_distribution<double> uniform;
    typedef std::normal_distribution<double> normal;
    typedef std::gamma_distribution<double> gamma;
    typedef std::bernoulli_distribution bernoulli;

    // Random number generator
    extern std::mt19937_64 rng;

}

#endif
