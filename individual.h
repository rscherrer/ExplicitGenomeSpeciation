/*==================================================================================================================================
                                                     individual.h
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


#ifndef __genomic_signatures_of_speciation__individual__
#define __genomic_signatures_of_speciation__individual__

#include <bitset>
#include <array>
#include <set>
#include <list>
#include <vector>
#include "ParameterSet.h"
#include "Population.h"
#include "Genome.h"

size_t divide_rounding_up( std::size_t dividend, std::size_t divisor );

std::string to_string( std::vector< bool > const & bitvector );

class Buffer;

inline double sqr(double x) { return x * x;}

class Individual {

    /*

    friend void decomposeVariance(int);
    friend void recordData(int, const std::array<size_t, 7u>&,
            const ParameterSet&,
            const std::vector<std::pair<double, double> >&,
            const std::vector<std::pair<double, double> >&,
            std::ofstream&, std::ofstream&, std::vector<std::pair<size_t, size_t> >&
            const std::list<PInd>&,
            const Genome& genome);
    friend void analyseNetwork(int);
    friend double computeMatingIsolation();
    friend double computePostIsolation();

     */

    friend class Buffer;

public:



    typedef std::pair<double, double> TradeOffPt;

    // A Trait object is locus- and individual-specific
    struct Trait
    {
        size_t alleleCount;
        double expression, geneticValue;
    };

    // Constructors
    Individual(const ParameterSet&, const Genome&);
    Individual(const std::vector<bool>&, const ParameterSet&, const Genome&);
    Individual(Individual const * const,
            Individual const * const,
            const ParameterSet&,
            const Population&,
            const Genome&);

    // Getters
    bool isFemale(const bool& isFemaleHeteroGamety) const {return isHeteroGamous == isFemaleHeteroGamety;}
    TradeOffPt getAttackRate() const { return attackRate; }
    double getBurnInRpSc(double) const;
    double getViability() const {return viability; }
    size_t getHabitat() const { return habitat; }
    size_t getEcotype() const { return ecotype; }
    std::vector<bool> getGenome() const { return genomeSequence; }
    std::vector<double> getTraitP() const { return traitP; } // size nCharacter
    std::vector<double> getTraitG() const { return traitG; } // size nCharacter
    std::vector<double> getTraitE() const { return traitE; } // size nCharacter
    std::vector<Trait> getTraitLocus() const { return traitLocus; } // size nLoci

    // Setters
    void disperse(const size_t& nHabitat) const { habitat = (habitat + 1u) % nHabitat; }
    size_t setEcotype(const Individual::TradeOffPt &threshold) const;
    void prepareChoice() const;
    bool acceptMate(Individual const * const, const ParameterSet&) const;


private:

    // Private setters
    void mutate(const ParameterSet&);
    void develop(const ParameterSet&, const Genome&);

    // Fields
    bool isHeteroGamous;
    mutable size_t habitat, ecotype;
    mutable std::list<double> obs;
    mutable double xsum, xxsum;
    std::vector<double> traitP, traitG, traitE; // size nCharacter
    double viability;
    TradeOffPt attackRate;
    std::vector<bool> genomeSequence;
    std::vector<Trait> traitLocus; // size nLoci

};

typedef Individual const * PInd;

bool tradeOffCompare (const Individual::TradeOffPt&, const Individual::TradeOffPt&);

#endif
