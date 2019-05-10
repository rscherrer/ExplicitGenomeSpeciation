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
#include "GeneticArchitecture.h"

class Buffer;
class Genome;
typedef std::pair<double, double> TradeOffPt;

inline double sqr(double x) { return x * x;}

class Individual {

    friend class Buffer;

public:

    struct Locus
    {
        size_t alleleCount;
        double expression;
        double locusGeneticValue;
    };

    // Constructors
    Individual(const ParameterSet&, const GeneticArchitecture&);
    Individual(const std::vector<bool>&, const ParameterSet&, const GeneticArchitecture&);
    Individual(Individual const * const, Individual const * const, const ParameterSet&, const GeneticArchitecture&);

    void recombineFreely(size_t&, size_t&, const size_t&, const double&, double&);
    void crossOver(size_t&, const double&, const double&, double&);
    void inheritLocus(Individual const * const, const bool&, const size_t&, const size_t&);
    void determineSex(const bool&, const bool&, const size_t&);
    void inheritGamete(Individual const * const, const ParameterSet&, const GeneticArchitecture&);

   


    bool isHeteroGamous;
    size_t habitat
    size_t ecotype;
    std::list<double> obs;
    double xsum;
    double xxsum;
    std::vector<double> phenotypes;
    std::vector<double> geneticValues;
    std::vector<double> envirValues;
    double viability;
    TradeOffPt attackRates;
    std::vector<bool> genomeSequence;
    std::vector<Locus> genotypes;

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
    size_t setEcotype(const TradeOffPt &threshold) const;
    void prepareChoice() const;
    bool acceptMate(Individual const * const, const ParameterSet&) const;

    void expressGene(const size_t&, const size_t&, const double&, const double&);
    void setAdditiveValue(const size_t&, const double&, const double&);
    void setEpistaticValue(const size_t&, const double&, const std::list<std::pair<size_t, double> >&)
    void setAttackRates(const double&);
    void setViability(const double&, const GeneticArchitecture&);
    void setPhenotype(const size_t&);
    void setEnvirValue(const size_t&, const double&);
    void setGeneticValue(const size_t&, const GeneticArchitecture&);
    void setLocusGeneticValue(const size_t&, const GeneticArchitecture&, const ParameterSet&);

    void setGenomeSequence(const size_t&, const double&);

private:

    // Private setters
    void mutate(const ParameterSet&);
    void develop(const ParameterSet&, const GeneticArchitecture&);


};

typedef Individual const * PInd;

// This functions compares two individuals on a trade-off line
// Should it be a member of class Individual then?
bool tradeOffCompare (const TradeOffPt &x, const TradeOffPt &y) {
    bool yOnLeft = y.first < y.second;
    if(x.first < x.second) {
        if(yOnLeft) return (x.first < y.first);
        else return true;
    }
    else {
        if(yOnLeft) return false;
        else return (x.second < y.second);
    }
}



#endif
