#ifndef EXPLICITGENOMESPECIATION_BUFFER_H
#define EXPLICITGENOMESPECIATION_BUFFER_H

#include <vector>
#include <iostream>
#include <fstream>

typedef std::vector<double> vecDbl;

class Buffer
{

    friend class MetaPop;

private:

    vecDbl fields;

public:

    Buffer() : fields({ }) {}
    ~Buffer() {}

    void flush();
    void add(const double&);
    void write(std::ofstream *&, const double&);
};

#endif
