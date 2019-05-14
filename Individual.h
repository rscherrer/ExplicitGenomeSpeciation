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

#include <list>
#include <vector>
#include "ParameterSet.h"
#include "Population.h"
#include "GeneticArchitecture.h"
#include "square.h"

class Buffer;
class Genome;

class Individual {

    friend class Buffer;
    friend class Population;

public:

    typedef Individual const * PInd;

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

    // Getters
    bool isFemale(const bool&) const;
    double getFitness() const { return fitness;}
    size_t getHabitat() const { return habitat; }
    size_t getEcotype() const { return ecotype; }
    std::pair<double, double> getAttackRates() const { return attackRates; }
    std::vector<size_t> getMates() const { return mates; }
    std::vector<bool> getGenomeSequence() const { return genomeSequence; }
    std::vector<double> getPhenotypes() const { return phenotypes; }
    std::vector<double> getGeneticValues() const { return geneticValues; }
    std::vector<double> getEnvirValues() const { return envirValues; }
    std::vector<Locus> getLoci() const { return genotypes; }

private:

    // Ecological attributes
    mutable size_t habitat;
    size_t ecotype;
    double fitness;
    double nOffspring;
    std::vector<size_t> mates;
    double matePreference;
    std::pair<double, double> attackRates;

    // Genetic attributes
    bool isHeterogamous;
    std::vector<bool> genomeSequence;
    std::vector<Locus> genotypes;
    std::vector<double> phenotypes;
    std::vector<double> geneticValues;
    std::vector<double> envirValues;

    // Ecology
    void disperse(const size_t& nHabitat) const;
    void setFitness(const std::pair<double, double>&);
    void setBurninFitness(const std::pair<double, double>&, const double&);
    void setAttackRates(const double&);
    void setMatePreference(const double&);
    void chooseMates(const double&, std::discrete_distribution<size_t>&, const std::vector<PInd>&, const ParameterSet&);
    double assessMatingProb(const double&, const double&) const;
    bool acceptMate(Individual const * const, const ParameterSet&) const;
    size_t sampleClutchSize(const double&) const;
    bool survive(const double&) const;

    // Genetics
    void mutate(const ParameterSet&);
    void develop(const ParameterSet&, const GeneticArchitecture&);
    void expressGene(const size_t&, const size_t&, const double&, const double&);
    void setAdditiveValue(const size_t&, const double&, const double&);
    void setEpistaticValue(const size_t&, const double&, const std::list<std::pair<size_t, double> >&);
    void setPhenotype(const size_t&);
    void setEnvirValue(const size_t&, const double&);
    void setGeneticValue(const size_t&, const GeneticArchitecture&);
    void setLocusGeneticValue(const size_t&, const GeneticArchitecture&, const ParameterSet&);
    void setGenomeSequence(const size_t&, const double&);
    void recombineFreely(size_t&, size_t&, const size_t&, const double&, double&) const;
    void crossOver(size_t&, const double&, const double&, double&) const;
    void inheritLocus(Individual const * const, const bool&, const size_t&, const size_t&);
    void determineSex(const bool&, const bool&, const size_t&);
    void inheritGamete(Individual const * const, const ParameterSet&, const GeneticArchitecture&);

    // To be taken care of
    size_t setEcotype(const std::pair<double, double> &threshold);

};

#endif
