/*==================================================================================================================================
                                                     random.h
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

#ifndef __genomic_signatures_of_speciation__random__
#define __genomic_signatures_of_speciation__random__

#include <random>

namespace rnd
{
	extern std::mt19937_64 rng;
	unsigned int set_seed();
	unsigned int set_seed(const unsigned int);
	int random_int(const int);
    size_t random_int(const size_t);
	bool bernoulli(const double = 0.5);
	int binomial(const int, const double = 0.5);
	size_t binomial(const size_t, const double = 0.5);
	int poisson(const double = 1.0);
    size_t geometric(const double);
	double uniform(const double = 1.0);
    double normal(const double = 0.0, const double = 1.0);
    double exponential(const double = 1.0);
    int sample(const double[], const int);
}

#endif