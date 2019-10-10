#ifndef EXPLICITGENOMESPECIATION_PARAM_H
#define EXPLICITGENOMESPECIATION_PARAM_H

#include "Types.h"
#include "Utilities.h"
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
        capacity(100.0),
        replenish(10.0),
        hsymmetry(1.0),
        ecosel(1.0),
        dispersal(1.0E-3),
        birth(4.0),
        survival(0.6),
        sexsel(10.0),
        matingcost(0.01),
        maxfeed(1.0),
        demesizes({ 10u, 0u }),        
        nloci(300u), // cannot be provided
        nvertices({ 100u, 100u, 100u }),
        nedges({ 0u, 0u, 0u }),
        nchrom(3u),
        mutation(1.0E-6),
        recombination(0.01),
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
        tburnin(10),
        tend(10),
        tsave(10),
        record(true),
        seed(makeDefaultSeed()),
        ntrials(10u)
    {
        // Make sure there are no more edges than feasible
        capEdges();

        // Make sure parameter values make sense
        checkParams();
    }

    void read(const std::string&);
    void update();

    void import(std::ifstream&);
    void capEdges();
    void checkParams();
    size_t makeDefaultSeed();

    // Ecological parameters    
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
    size_t nloci;
    vecUns nvertices;
    vecUns nedges;
    size_t nchrom;
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
    bool record;
    size_t seed;
    size_t ntrials;

};

#endif
