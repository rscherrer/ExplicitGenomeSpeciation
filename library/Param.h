#ifndef EXPLICITGENOMESPECIATION_PARAM_H
#define EXPLICITGENOMESPECIATION_PARAM_H


#include "Utilities.h"
#include "Random.h"
#include <fstream>
#include <iostream>
#include <chrono>
#include <cassert>
#include <cstdint>

// Parameter set. Contains values of the parameters of the simulation.
// All parameters have default values that can be modified by calling
// the program with a parameter file name as unique argument.

struct Param {

    Param();

    void read(const std::string&);
    void update();
    void import(std::ifstream&);

    void write(std::ofstream&) const;
    void save() const;
    void check() const;

    size_t makeDefaultSeed();

    // Ecological parameters    
    size_t rdynamics;
    double capacity;
    double replenish;
    double inflow;
    double outflow;
    double hsymmetry;
    double ecosel;
    double dispersal;
    double birth;
    double survival;
    double sexsel;
    double matingcost;
    double ecoscale;
    std::vector<size_t> demesizes;

    // Genetic parameters
    mutable size_t nloci;
    mutable std::vector<size_t> nvertices;
    mutable std::vector<size_t> nedges;
    mutable size_t nchrom;
    double  mutation;
    double  recombination;
    double  allfreq;

    // Genotype-phenotype map
    std::vector<double> scaleA;
    std::vector<double> scaleD;
    std::vector<double> scaleI;
    std::vector<double> scaleE;
    std::vector<double> locusE;
    std::vector<double> skews;
    double effectshape;
    double effectscale;
    double interactionshape;
    double interactionscale;
    double dominancevar;

    // Simulation parameters
    int tburnin;
    int tend;
    int tsave;
    int tfreeze;
    int tpedigree;
    bool talkative;
    bool record;
    bool datsave;
    bool choosewhattosave;
    bool gensave;
    bool archsave;
    bool archload;
    bool parsave;
    bool pedigreesave;
    std::string archfile;
    std::string parfile;
    std::string orderfile;
    std::string logfile;
    std::string freezerfile;
    std::string locifile;
    std::string pedigreefile;
    size_t seed;
    size_t ntrials;
    size_t pedigreetrials;
    size_t pedigreeoffspring;

};

#endif
