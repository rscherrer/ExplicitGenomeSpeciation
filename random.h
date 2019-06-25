#ifndef __genomic_signatures_of_speciation__random__
#define __genomic_signatures_of_speciation__random__

#include <random>

namespace rnd
{
	extern std::mt19937_64 rng;
	int random_int(const int);
    size_t random_int(const size_t);
	bool bernoulli(const double = 0.5);
	int binomial(const int, const double = 0.5);
	size_t binomial(const size_t, const double = 0.5);
	size_t poisson(const double = 1.0);
    size_t geometric(const double);
	double uniform(const double = 1.0);
    double normal(const double = 0.0, const double = 1.0);
    double exponential(const double = 1.0);
    int sample(const double[], const int);
}

#endif