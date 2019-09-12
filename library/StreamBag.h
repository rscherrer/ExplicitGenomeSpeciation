#ifndef EXPLICITGENOMESPECIATION_STREAMBAG_H
#define EXPLICITGENOMESPECIATION_STREAMBAG_H

#include "MetaPop.h"
#include "Buffer.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

typedef std::vector<std::ofstream *> vecStreams;
typedef std::vector<std::string> vecStrings;

class StreamBag
{

    friend class MetaPop;
    friend class Buffer;

private:

    vecStreams files;
    vecStrings names;

public:

    StreamBag() :
        files({ }),
        names({
              "time",
              "ecotype0",
              "ecotype1",
              "popsize0",
              "popsize1",
              "nfemales0",
              "nfemales1",
              "resource00",
              "resource01",
              "resource10",
              "resource11",
              "meanx0",
              "meanx1",
              "meany0",
              "meany1",
              "meanz0",
              "meanz1",
              "ecological_iso",
              "spatial_iso",
              "mating_iso"
        })
    {}
    ~StreamBag() {}

    void openAll();
    void closeAll();

};

#endif
