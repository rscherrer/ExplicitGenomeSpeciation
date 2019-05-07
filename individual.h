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
#include "ParameterSet.h"

/*
const size_t nEcoLoci         = 400u;
const size_t nMatLoci         = 200u;
const size_t nNtrLoci         = 400u;
const size_t nEcoInteractions = 1000u;
const size_t nMatInteractions = 500u;
const size_t nNtrInteractions = 0u;
const size_t nChromosomes     = 3u;
const size_t nHabitat         = 2u;
const size_t nCharacter       = 3u;
const double tiny             = 1.0e-12;    // for clipping towards zero
const size_t nLoci = nEcoLoci + nMatLoci + nNtrLoci;
const size_t nBits = 2u * nLoci;



extern bool isFemaleHeteroGamety;

*/

class Buffer;

inline double sqr(double x) { return x * x;}

class Individual {
    friend void decomposeVariance(int);
    friend void recordData(int, const std::array<size_t, 7u>&);
    friend void analyseNetwork(int);
    friend double computeMatingIsolation();
    friend double computePostIsolation();
    friend class Buffer;
public:
    typedef std::pair<double, double> TradeOffPt;
    struct Character
    {
        size_t character, linkageGroup;
        double location, effectSize, dominanceCoeff, avgEffectOfSubstitution,
            varD, F_it, F_is, F_st, P_st, G_st, Q_st, C_st;
        std::array<double, 3u> alleleFrequency, meanEffect, varP, varG, varA, varI;
        std::list<std::pair<size_t, double> > edges;
    };
    struct Trait
    {
        size_t alleleCount;
        double expression, geneticValue;
    };
    Individual(const ParameterSet&);
    Individual(const std::string&, const ParameterSet&);
    Individual(Individual const * const, Individual const * const, const ParameterSet&);
    void disperse() const { habitat = (habitat + 1u) % nHabitat; }
    bool isFemale(const bool& isFemaleHeteroGamety) const {return isHeteroGamous == isFemaleHeteroGamety;}
    std::string getSequence() const { return genome.to_string();}
    TradeOffPt getAttackRate() const { return attackRate; }
    double getBurnInRpSc(double) const;
    double getViability() const {return viability; }
    void prepareChoice() const;
    bool acceptMate(Individual const * const, const ParameterSet&) const;
    size_t getHabitat() const { return habitat; }
    size_t getEcotype() const { return ecotype; }

    static void setNBits(const size_t& nbits) { nBits = nbits; }

    std::bitset<nBits> getGenome() { return genome; }

    std::array<double, nCharacter> getTraitP() const { return traitP; }
    std::array<double, nCharacter> getTraitG() const { return traitG; }
    std::array<double, nCharacter> getTraitE() const { return traitE; }
    std::array<Trait, nLoci> getTraitLocus() const { return traitLocus; }
    size_t setEcotype(const Individual::TradeOffPt &threshold) const;
    static void generateGeneticArchitecture(const ParameterSet&);
    static void storeGeneticArchitecture(const std::string&);
    static void loadGeneticArchitecture(const std::string&);
    static std::array<std::array<double, 3u>, nCharacter>
            avgG, varP, varG, varA, varI;
    static std::array<double, nCharacter> varD, F_st, P_st, G_st, Q_st, C_st;
    static std::array<double, nChromosomes - 1u> chromosomeSize;
    static std::array<std::set<size_t>, nCharacter> vertices;
    static std::array<Character, nLoci> characterLocus;
private:
    static size_t nBits;
    void mutate(const ParameterSet&);
    void develop(const ParameterSet&);
    bool isHeteroGamous;
    mutable size_t habitat, ecotype;
    mutable std::list<double> obs;
    mutable double xsum, xxsum;
    std::array<double, nCharacter> traitP, traitG, traitE;
    double viability;
    TradeOffPt attackRate;
    std::bitset<nBits> genome;
    std::array<Trait, nLoci> traitLocus;

};

typedef Individual const * PInd;
bool tradeOffCompare (const Individual::TradeOffPt&, const Individual::TradeOffPt&);

#endif
