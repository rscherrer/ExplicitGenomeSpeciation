#include "random.h"
#include <chrono>

namespace rnd
{
	std::mt19937_64 rng;

	bool bernoulli(const double p)
	{
		return std::bernoulli_distribution(p)(rng);
	}

	int binomial(const int n, const double p)
	{
		return std::binomial_distribution<int>(n, p)(rng);
	}

	size_t binomial(const size_t n, const double p)
	{
		return std::binomial_distribution<size_t>(n, p)(rng);
	}

	size_t poisson(const double lambda)
	{
		return std::poisson_distribution<int>(lambda)(rng);
	}
    
    size_t geometric(const double p)
    {
        return std::geometric_distribution<size_t>(p)(rng);
    }
    
    int random_int(const int n)
    {
        return std::uniform_int_distribution<int>(0, n - 1)(rng);
    }
    
    size_t random_int(const size_t n)
    {
        return std::uniform_int_distribution<size_t>(0u, n - 1u)(rng);
    }

	double uniform(const double max)
	{
		return std::uniform_real_distribution<double>(0.0, max)(rng);
	}
    
    double normal(const double mean, const double stdev)
	{
		return stdev == 0.0 ? 0.0 : std::normal_distribution<double>(mean, stdev)(rng);
	}
    
    double exponential(const double lambda)
	{
		return std::exponential_distribution<double>(lambda)(rng);
	}

	int sample(const double pdf[], const int J)
    {
        return std::discrete_distribution<int>(pdf, pdf + J)(rng);
    }
}