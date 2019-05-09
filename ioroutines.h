//
// Created by p278834 on 9-5-2019.
//

#ifndef EXPLICITGENOMESPECIATION_IOROUTINES_H
#define EXPLICITGENOMESPECIATION_IOROUTINES_H

#include <string>
#include <iostream>
#include <fstream>
#include "ParameterSet.h"
#include "random.h"

namespace io {


    void writeParameters(std::ofstream&, const ParameterSet& , const char = ' ');

}



#endif //EXPLICITGENOMESPECIATION_IOROUTINES_H
