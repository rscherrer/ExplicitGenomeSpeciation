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
              "mean_eco0",
              "mean_eco1",
              "mean_eco",
              "varP_eco",
              "varG_eco",
              "varA_eco",
              "mean_mat0",
              "mean_mat1",
              "mean_mat",
              "varP_mat",
              "varG_mat",
              "varA_mat",
              "mean_ntr0",
              "mean_ntr1",
              "mean_ntr",
              "varP_ntr",
              "varG_ntr",
              "varA_ntr",
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
