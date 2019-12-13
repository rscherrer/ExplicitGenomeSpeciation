#ifndef EXPLICITGENOMESPECIATION_PARAM_H
#define EXPLICITGENOMESPECIATION_PARAM_H

#include "Types.h"
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

    Param() :
        rdynamics(1u),
        trenewal(0.001),
        capacity(100.0),
        replenish(1.0),
        hsymmetry(0.0),
        ecosel(1.8),
        dispersal(1.0E-2),
        birth(4.0),
        survival(0.8),
        sexsel(10.0),
        matingcost(0.01),
        maxfeed(4.0E-4),
        demesizes({ 100u, 0u }),
        nloci(90u), // cannot be provided
        nvertices({ 30u, 30u, 30u }),
        nedges({ 0u, 0u, 0u }),
        nchrom(3u),
        mutation(1.0E-3),
        recombination(3.0),
        allfreq(0.2),
        scaleA({ 1.0, 1.0, 1.0 }),
        scaleD({ 0.0, 0.0, 0.0 }),
        scaleI({ 0.0, 0.0, 0.0 }),
        scaleE({ 0.0, 0.0, 0.0 }),
        locusE({ 0.0, 0.0, 0.0 }), // cannot be provided
        skews({ 1.0, 1.0, 1.0 }),
        effectshape(2.0),
        effectscale(1.0),
        interactionshape(5.0),
        interactionscale(1.0),
        dominancevar(1.0),
        tburnin(0),
        tend(10),
        tsave(10),
        talkative(true),
        record(true),
        archsave(false),
        archload(false),
        parsave(true),
        archfile("architecture.txt"),        
        parfile("paramlog.txt"),
        seed(makeDefaultSeed()),
        ntrials(100u)
    {

        // Make sure parameter values make sense
        checkParams();

        // Seed the random number generator
        rnd::rng.seed(seed);
    }

    void read(const std::string&);
    void update();
    void import(std::ifstream&);

    void write(std::ofstream&) const;
    void save() const;
    void checkParams() const;

    size_t makeDefaultSeed();

    // Ecological parameters    
    size_t rdynamics;
    double trenewal;
    double capacity;
    double replenish;
    double hsymmetry;
    double ecosel;
    double dispersal;
    double birth;
    double survival;
    double sexsel;
    double matingcost;
    double maxfeed;
    vecUns demesizes;

    // Genetic parameters
    mutable size_t nloci;
    mutable vecUns nvertices;
    mutable vecUns nedges;
    mutable size_t nchrom;
    double  mutation;
    double  recombination;
    double  allfreq;

    // Genotype-phenotype map
    vecDbl scaleA;
    vecDbl scaleD;
    vecDbl scaleI;
    vecDbl scaleE;
    vecDbl locusE;
    vecDbl skews;
    double effectshape;
    double effectscale;
    double interactionshape;
    double interactionscale;
    double dominancevar;

    // Simulation parameters
    int  tburnin;
    int  tend;
    int  tsave;
    bool talkative;
    bool record;
    bool archsave;
    bool archload;
    bool parsave;
    std::string archfile;
    std::string parfile;
    size_t seed;
    size_t ntrials;

};

#endif
