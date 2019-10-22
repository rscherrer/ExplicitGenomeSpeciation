#include "Random.h"

/// Random number generator
std::mt19937_64 rnd::rng;

/// Function to sample a binary event
bool rnd::bernoulli(const double &p)
{
    return std::bernoulli_distribution(p)(rng);
}


/// Function to sample a number of successes
size_t rnd::binomial(const size_t &n, const double &p)
{
    return std::binomial_distribution<size_t>(n, p)(rng);
}


/// Function to sample from a Poisson distribution
size_t rnd::poisson(const double &lambda)
{
    return std::poisson_distribution<size_t>(lambda)(rng);
}


/// Function to sample from a geometric distribution
size_t rnd::geometric(const double &p)
{
    return std::geometric_distribution<size_t>(p)(rng);
}


/// Function to sample from a uniform integer distribution
size_t rnd::random(const size_t &n)
{
    return std::uniform_int_distribution<size_t>(0, n - 1)(rng);
}


/// Function to sample from a uniform distribution
double rnd::uniform(const double &max)
{
    return std::uniform_real_distribution<double>(0.0, max)(rng);
}


/// Function to sample from a normal distribution
double rnd::normal(const double &mean, const double &sdev)
{
    return sdev == 0.0 ? 0.0 : std::normal_distribution<double>(mean, sdev)(
     rng);
}


/// Half-normal distribution
double rnd::hnormal(const double &sdev)
{
    double x = rnd::normal(0.0, sdev);
    return x > 0.0 ? x : -1.0 * x;
}


/// Function to sample from an exponential distribution
double rnd::exponential(const double &lambda)
{
    return std::exponential_distribution<double>(lambda)(rng);
}


/// Discrete distribution from vector of probabilities
size_t rnd::sample(const vecDbl &probs)
{
    return std::discrete_distribution<size_t>(probs.begin(), probs.end())(rng);
}


/// Flip a number with a certain probability
double rnd::flip(const double &x, const double &p)
{
    return rnd::bernoulli(p) ? x : -1.0 * x;
}


/// Gamma distribution
double rnd::gamma(const double &shape, const double &scale)
{
    return std::gamma_distribution<double>(shape, scale)(rng);
}


/// Two-sided gamma
double rnd::bigamma(const double &shape, const double &scale)
{
    return rnd::flip(rnd::gamma(shape, scale), 0.5);
}


