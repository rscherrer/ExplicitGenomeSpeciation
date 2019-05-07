//
// Created by p278834 on 7-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_BUFFERBOX_H
#define EXPLICITGENOMESPECIATION_BUFFERBOX_H

#include "individual.h"

class BufferBox {

public:

    Buffer *bufferFreq;
    Buffer *bufferF_it;
    Buffer *bufferF_is;
    Buffer *bufferF_st;
    Buffer *bufferP_st;
    Buffer *bufferG_st;
    Buffer *bufferQ_st;
    Buffer *bufferC_st;
    Buffer *bufferVarP;
    Buffer *bufferVarG;
    Buffer *bufferVarA;
    Buffer *bufferVarD;
    Buffer *bufferVarI;

};


#endif //EXPLICITGENOMESPECIATION_BUFFERBOX_H
