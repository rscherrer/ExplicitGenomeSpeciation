#ifndef EXPLICITGENOMESPECIATION_TESTUTILITIES_H
#define EXPLICITGENOMESPECIATION_TESTUTILITIES_H

#include <boost/test/included/unit_test.hpp>
#include <vector>

namespace tst
{

    void makeValidParamFile();
    void makeValidParamFile2();
    void makeInvalidParamName();
    void makeInvalidParamValue();
    void makeInvalidParamValue2();
    void makeParamFileWithArchitecture();
    void makeParamFileWithMissingArchitecture();

}

void tst::makeValidParamFile()
{
    std::ofstream file;
    file.open("validparamfile.txt");
    if (!file.is_open())
        std::cout << "Unable to open valid parameter test file.\n";

    file << "rdynamics" << '\t' << 0 << '\n'
         << "trenewal" << '\t' << 0.001 << '\n'
         << "capacity" << '\t' << 10.0 << '\n'
         << "replenish" << '\t' << 1.0 << '\n'
         << "hsymmetry" << '\t' << 1.0 << '\n'
         << "ecosel" << '\t' << 1.0 << '\n'
         << "dispersal" << '\t' << 1.0E-3 << '\n'
         << "birth" << '\t' << 1.5 << '\n'
         << "survival" << '\t' << 0.8 << '\n'
         << "sexsel" << '\t' << 10.0 << '\n'
         << "matingcost" << '\t' << 0.01 << '\n'
         << "maxfeed" << '\t' << 4.0E-4 << '\n'
         << "demesizes" << '\t' << 100 << '\t' << 100 << '\n'
         << "nvertices" << '\t' << 40 << '\t' << 20 << '\t' << 40 << '\n'
         << "nedges" << '\t' << 1000 << '\t' << 500 << '\t' << 0 << '\n'
         << "nchrom" << '\t' << 3 << '\n'
         << "mutation" << '\t' << 1.0e-5 << '\n'
         << "recombination" << '\t' << 0.01 << '\n'
         << "allfreq" << '\t' << 0.02 << '\n'
         << "scaleA" << '\t' << 1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
         << "scaleD" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
         << "scaleI" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
         << "scaleE" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
         << "skews" << '\t' << 1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
         << "effectshape" << '\t' << 2.0 << '\n'
         << "effectscale" << '\t' << 1.0 << '\n'
         << "interactionshape" << '\t' << 5.0 << '\n'
         << "interactionscale" << '\t' << 1.0 << '\n'
         << "dominancevar" << '\t' << 1.0 << '\n'
         << "tburnin" << '\t' << 5 << '\n'
         << "tend" << '\t' << 2 << '\n'
         << "tsave" << '\t' << 1 << '\n'
         << "record" << '\t' << 1 << '\n'
         << "seed" << '\t' << 42 << '\n'
         << "ntrials" << '\t' << 20 << '\n';

    file.close();
}

void tst::makeValidParamFile2()
{
    std::ofstream file;
    file.open("validparamfile2.txt");
    if (!file.is_open())
        std::cout << "Unable to open valid parameter test file.\n";

    file << "rdynamics" << '\t' << 1 << '\n'
         << "trenewal" << '\t' << 0.001 << '\n'
         << "capacity" << '\t' << 10.0 << '\n'
         << "replenish" << '\t' << 1.0 << '\n'
         << "hsymmetry" << '\t' << 0.0 << '\n'
         << "ecosel" << '\t' << 1.8 << '\n'
         << "dispersal" << '\t' << 1.0E-2 << '\n'
         << "birth" << '\t' << 4.0 << '\n'
         << "survival" << '\t' << 0.8 << '\n'
         << "sexsel" << '\t' << 10.0 << '\n'
         << "matingcost" << '\t' << 0.01 << '\n'
         << "maxfeed" << '\t' << 4.0E-4 << '\n'
         << "demesizes" << '\t' << 100 << '\t' << 0 << '\n'
         << "nvertices" << '\t' << 40 << '\t' << 40 << '\t' << 40 << '\n'
         << "nedges" << '\t' << 0 << '\t' << 0 << '\t' << 0 << '\n'
         << "nchrom" << '\t' << 3 << '\n'
         << "mutation" << '\t' << 1.0e-3 << '\n'
         << "recombination" << '\t' << 0.01 << '\n'
         << "allfreq" << '\t' << 0.2 << '\n'
         << "scaleA" << '\t' << 1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
         << "scaleD" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
         << "scaleI" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
         << "scaleE" << '\t' << 0.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
         << "skews" << '\t' << 1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
         << "effectshape" << '\t' << 2.0 << '\n'
         << "effectscale" << '\t' << 1.0 << '\n'
         << "interactionshape" << '\t' << 5.0 << '\n'
         << "interactionscale" << '\t' << 1.0 << '\n'
         << "dominancevar" << '\t' << 1.0 << '\n'
         << "tburnin" << '\t' << 100 << '\n'
         << "tend" << '\t' << 100 << '\n'
         << "tsave" << '\t' << 50 << '\n'
         << "record" << '\t' << 1 << '\n'
         << "ntrials" << '\t' << 100 << '\n';

    file.close();
}


void tst::makeInvalidParamName()
{
    std::ofstream file;
    file.open("invalidparamname.txt");
    if (!file.is_open())
        std::cout << "Unable to open invalid parameter name test file.\n";
    file << "nonsense" << '\t' << 3.0 << '\n';
    file.close();
}

void tst::makeInvalidParamValue()
{
    std::ofstream file;
    file.open("invalidparamvalue.txt");
    if (!file.is_open())
        std::cout << "Unable to open invalid parameter value test file.\n";

    file << "trenewal" << '\t' << -1.0 << '\n'
         << "capacity" << '\t' << -1.0 << '\n'
         << "replenish" << '\t' << -1.0 << '\n'
         << "hsymmetry" << '\t' << -1.0 << '\n'
         << "ecosel" << '\t' << -1.0 << '\n'
         << "dispersal" << '\t' << 10 << '\n'
         << "birth" << '\t' << -4.0 << '\n'
         << "survival" << '\t' << -0.8 << '\n'
         << "sexsel" << '\t' << -10.0 << '\n'
         << "matingcost" << '\t' << -0.01 << '\n'
         << "maxfeed" << '\t' << -1.0 << '\n'
         << "demesizes" << '\t' << -100 << '\t' << -100 << '\n'
         << "nvertices" << '\t' << 1 << '\t' << 1 << '\t' << 1 << '\n'
         << "nchrom" << '\t' << 0 << '\n'
         << "mutation" << '\t' << -1.0e-5 << '\n'
         << "recombination" << '\t' << -0.01 << '\n'
         << "allfreq" << '\t' << -0.02 << '\n'
         << "scaleA" << '\t' << -1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
         << "scaleD" << '\t' << -1.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
         << "scaleI" << '\t' << -1.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
         << "scaleE" << '\t' << -1.0 << '\t' << 0.0 << '\t' << 0.0 << '\n'
         << "skews" << '\t' << -1.0 << '\t' << 1.0 << '\t' << 1.0 << '\n'
         << "effectshape" << '\t' << -2.0 << '\n'
         << "effectscale" << '\t' << -1.0 << '\n'
         << "interactionshape" << '\t' << -5.0 << '\n'
         << "interactionscale" << '\t' << -1.0 << '\n'
         << "dominancevar" << '\t' << -1.0 << '\n'
         << "ntrials" << '\t' << 0 << '\n';

    file.close();
}

void tst::makeInvalidParamValue2()
{
    std::ofstream file;
    file.open("invalidparamvalue2.txt");
    if (!file.is_open())
        std::cout << "Unable to open invalid parameter value test file.\n";

    file << "dispersal" << '\t' << -1.0e-3 << '\n'
         << "survival" << '\t' << 2.0 << '\n'
         << "mutation" << '\t' << 10 << '\n'
         << "hsymmetry" << '\t' << 2.0 << '\n'
         << "allfreq" << '\t' << 10 << '\n';

    file.close();
}

#endif
