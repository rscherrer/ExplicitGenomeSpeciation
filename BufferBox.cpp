//
// Created by p278834 on 7-5-2019.
//

#include "BufferBox.h"

// Constructor and initializer
BufferBox::BufferBox(const ParameterSet& parameters, const Population& population, const Genome& genome)
{

    bufferFreq = new Buffer("freq", parameters, population, genome);
    bufferF_it = new Buffer("Fit", parameters, population, genome);
    bufferF_is = new Buffer("Fis", parameters, population, genome);
    bufferF_st = new Buffer("Fst", parameters, population, genome);
    bufferP_st = new Buffer("Pst", parameters, population, genome);
    bufferG_st = new Buffer("Gst", parameters, population, genome);
    bufferQ_st = new Buffer("Qst", parameters, population, genome);
    bufferC_st = new Buffer("Cst", parameters, population, genome);
    bufferVarP = new Buffer("varP", parameters, population, genome);
    bufferVarG = new Buffer("varG", parameters, population, genome);
    bufferVarA = new Buffer("varA", parameters, population, genome);
    bufferVarD = new Buffer("varD", parameters, population, genome);
    bufferVarI = new Buffer("varI", parameters, population, genome);

}
