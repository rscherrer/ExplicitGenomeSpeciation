/*==================================================================================================================================
                                                     analysis.h
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


#ifndef __genomic_signatures_of_speciation__analysis__
#define __genomic_signatures_of_speciation__analysis__

#include "Individual.h"
#include "BufferBox.h"
#include "Genome.h"


void decomposeVariance(int,
        const ParameterSet&,
        BufferBox&,
        std::ofstream&,
        std::ofstream&,
        Population&,
        Genome&);

void analyseNetwork(int, const ParameterSet&, const Population&, const Genome&);

#endif
