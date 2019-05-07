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

#include "individual.h"
#include "BufferBox.h"

class Buffer {
public:
    Buffer(const std::string&, const ParameterSet&);
    ~Buffer() { ofs.close(); }
    double &operator[](size_t j) {return data[i][j];}
    void flush(const ParameterSet&);
private:
    size_t i, k;
    int t;
    const size_t n;
    const std::string label;
    std::ofstream ofs;
    std::vector< std::vector<double> > data;
    static const char sep;
};

void decomposeVariance(int,
        const ParameterSet&,
        BufferBox&,
        const Individual::TradeOffPt&,
        const std::array<std::pair<double, double>, nHabitat>&,
        const std::array<std::pair<double, double>, nHabitat>&,
        std::ofstream&,
        std::ofstream&,
        std::array<std::pair<size_t, size_t>, nHabitat>&,
        const std::list<PInd>&);
void analyseNetwork(int, const ParameterSet&, const std::list<PInd>&);

#endif
