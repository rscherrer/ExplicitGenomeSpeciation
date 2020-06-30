/*==================================================================================================================================
                                                     random.cpp
====================================================================================================================================

C++-code accompanying:	
		 
		(ms. in prep).

Written by:
        G. Sander van Doorn
       	Centre for Ecological and Evolutionary Studies - Theoretical Biology Group
        University of Groningen
        the Netherlands

Program version
		xx/xx/2018	:	

Instructions for compiling and running the program
		
	Versions of this program were compiled and run on Windows and Mac, using Microsoft Visual C++
	2010 and XCode, respectively. The code is written in standard C++ and should be compatible with 
	other compilers. 

=================================================================================================================================*/

#include "random.h"
#include <chrono>
#include <sstream>

namespace rnd
{
	std::mt19937_64 rng;
    
    unsigned int set_seed()
	{
        unsigned int seed = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
		return set_seed(seed);
	}

	unsigned int set_seed(const unsigned int seed)
	{
		rng.seed(seed);
        return seed;
	}

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

	int poisson(const double lambda) 
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