#ifndef EXPLICITGENOMESPECIATION_OUTPUT_H
#define EXPLICITGENOMESPECIATION_OUTPUT_H

#include "Types.h"
//#include <vector>
//#include <iostream>
#include <fstream>
//#include <string>

typedef std::vector<std::ofstream *> vecStreams;

class Output
{

    friend class Stats;

private:

    vecStreams files;
    vecStrings names;

public:

    Output() :
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
              "varD_eco",
              "varI_eco",
              "varN_eco",
              "Pst_eco",
              "Gst_eco",
              "Qst_eco",
              "Cst_eco",
              "Fst_eco",
              "mean_mat0",
              "mean_mat1",
              "mean_mat",
              "varP_mat",
              "varG_mat",
              "varA_mat",
              "varD_mat",
              "varI_mat",
              "varN_mat",
              "Pst_mat",
              "Gst_mat",
              "Qst_mat",
              "Cst_mat",
              "Fst_mat",
              "mean_ntr0",
              "mean_ntr1",
              "mean_ntr",
              "varP_ntr",
              "varG_ntr",
              "varA_ntr",
              "varD_ntr",
              "varI_ntr",
              "varN_ntr",
              "Pst_ntr",
              "Gst_ntr",
              "Qst_ntr",
              "Cst_ntr",
              "Fst_ntr",
              "ecological_iso",
              "spatial_iso",
              "mating_iso",
              "varP_scan",
              "varG_scan",
              "varA_scan",
              "varD_scan",
              "varI_scan",
              "varN_scan",
              "Pst_scan",
              "Gst_scan",
              "Qst_scan",
              "Cst_scan",
              "Fst_scan"
        })
    {}
    ~Output() {}

    void openAll();
    void closeAll();

};

#endif
