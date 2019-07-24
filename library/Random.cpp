#include "Random.h"


/// Function to make a random number generator
std::mt19937_64 Random::makeRandomGenerator(const size_t &seed)
{
    std::mt19937_64 generator;

    generator.seed(seed);

    return generator;
}


/// Function to sample a binary event
bool Random::bernoulli(const double &p)
{
    return std::bernoulli_distribution(p)(rng);
}


/// Function to sample a number of successes
size_t Random::binomial(const size_t &n, const double &p)
{
    return std::binomial_distribution<size_t>(n, p)(rng);
}


/// Function to sample from a Poisson distribution
size_t Random::poisson(const double &lambda)
{
    return std::poisson_distribution<size_t>(lambda)(rng);
}


/// Function to sample from a geometric distribution
size_t Random::geometric(const double &p)
{
    return std::geometric_distribution<size_t>(p)(rng);
}


/// Function to sample from a uniform integer distribution
size_t Random::random(const size_t &n)
{
    return std::uniform_int_distribution<size_t>(0, n - 1)(rng);
}


/// Function to sample from a uniform distribution
double Random::uniform(const double &max)
{
    return std::uniform_real_distribution<double>(0.0, max)(rng);
}


/// Function to sample from a normal distribution
double Random::normal(const double &mean, const double &stdev)
{
    return stdev == 0.0 ? 0.0 : std::normal_distribution<double>(mean, stdev)(rng);
}


/// Function to sample from an exponential distribution
double Random::exponential(const double &lambda)
{
    return std::exponential_distribution<double>(lambda)(rng);
}


/// Function to sample from a discrete distribution
size_t Random::sample(const double pdf[], const size_t &J)
{
    return std::discrete_distribution<size_t>(pdf, pdf + J)(rng);
}
