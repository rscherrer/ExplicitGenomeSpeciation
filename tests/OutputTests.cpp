#include "library/Collector.h"
#include "library/Random.h"
#include "library/Utilities.h"
#include "tests/TestUtilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cassert>
#include <fstream>

// Test the output files

// One big test case in which we make sure everything is saved properly

BOOST_AUTO_TEST_CASE(OutputFilesAreCorrectlyWritten)
{
    std::clog << "Testing output files...\n";

    Param pars;
    pars.rdynamics = 1u;
    pars.dispersal = 0.0;
    pars.birth = 0.0;
    pars.survival = 1.0;
    pars.trenewal = 0.001;
    pars.hsymmetry = 1.0;
    pars.demesizes = { 100u, 0u };
    pars.allfreq = 0.5;
    pars.ecosel = 1.0;
    pars.tburnin = 0u;

    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    Collector collector = Collector(arch);

    size_t cumulsize = 0u;
    size_t lastgenfirst = 0u;

    for (int t = 0; t < 10; ++t) {

        metapop.cycle(pars, arch);
        collector.analyze(metapop, pars, arch);
        collector.print(t, metapop);

        if (t != 9) lastgenfirst += metapop.getSize();

        cumulsize += metapop.getSize();

    }
    collector.shutdown();

    // Read output files
    vecDbl time = tst::readfile("time.dat");
    vecDbl count00 = tst::readfile("count00.dat");
    vecDbl count01 = tst::readfile("count01.dat");
    vecDbl count10 = tst::readfile("count10.dat");
    vecDbl count11 = tst::readfile("count11.dat");
    vecDbl fem0 = tst::readfile("fem0.dat");
    vecDbl fem1 = tst::readfile("fem1.dat");
    vecDbl resource00 = tst::readfile("resource00.dat");
    vecDbl resource01 = tst::readfile("resource01.dat");
    vecDbl resource10 = tst::readfile("resource10.dat");
    vecDbl resource11 = tst::readfile("resource11.dat");
    vecDbl mean00x = tst::readfile("mean00_x.dat");
    vecDbl mean01x = tst::readfile("mean01_x.dat");
    vecDbl mean10x = tst::readfile("mean10_x.dat");
    vecDbl mean11x = tst::readfile("mean11_x.dat");
    vecDbl varpx = tst::readfile("varP_x.dat");
    vecDbl vargx = tst::readfile("varG_x.dat");
    vecDbl varax = tst::readfile("varA_x.dat");
    vecDbl vardx = tst::readfile("varD_x.dat");
    vecDbl varix = tst::readfile("varI_x.dat");
    vecDbl varnx = tst::readfile("varN_x.dat");
    vecDbl pstx = tst::readfile("Pst_x.dat");
    vecDbl gstx = tst::readfile("Gst_x.dat");
    vecDbl qstx = tst::readfile("Qst_x.dat");
    vecDbl cstx = tst::readfile("Cst_x.dat");
    vecDbl fstx = tst::readfile("Fst_x.dat");
    vecDbl ei = tst::readfile("EI.dat");
    vecDbl si = tst::readfile("SI.dat");
    vecDbl ri = tst::readfile("RI.dat");

    // Check output files

    BOOST_CHECK_EQUAL(time.size(), 10u);
    BOOST_CHECK_EQUAL(count00.size(), 10u);
    BOOST_CHECK_EQUAL(count01.size(), 10u);
    BOOST_CHECK_EQUAL(count10.size(), 10u);
    BOOST_CHECK_EQUAL(count11.size(), 10u);
    BOOST_CHECK_EQUAL(fem0.size(), 10u);
    BOOST_CHECK_EQUAL(fem1.size(), 10u);
    BOOST_CHECK_EQUAL(resource00.size(), 10u);
    BOOST_CHECK_EQUAL(resource01.size(), 10u);
    BOOST_CHECK_EQUAL(resource10.size(), 10u);
    BOOST_CHECK_EQUAL(resource11.size(), 10u);

    BOOST_CHECK_EQUAL(mean00x.size(), 10u);
    BOOST_CHECK_EQUAL(mean01x.size(), 10u);
    BOOST_CHECK_EQUAL(mean10x.size(), 10u);
    BOOST_CHECK_EQUAL(mean11x.size(), 10u);
    BOOST_CHECK_EQUAL(varpx.size(), 10u);
    BOOST_CHECK_EQUAL(vargx.size(), 10u);
    BOOST_CHECK_EQUAL(varax.size(), 10u);
    BOOST_CHECK_EQUAL(vardx.size(), 10u);
    BOOST_CHECK_EQUAL(varix.size(), 10u);
    BOOST_CHECK_EQUAL(varnx.size(), 10u);
    BOOST_CHECK_EQUAL(pstx.size(), 10u);
    BOOST_CHECK_EQUAL(gstx.size(), 10u);
    BOOST_CHECK_EQUAL(qstx.size(), 10u);
    BOOST_CHECK_EQUAL(cstx.size(), 10u);
    BOOST_CHECK_EQUAL(fstx.size(), 10u);

    BOOST_CHECK_EQUAL(ei.size(), 10u);
    BOOST_CHECK_EQUAL(si.size(), 10u);
    BOOST_CHECK_EQUAL(ri.size(), 10u);

    // Read individual trait values and ecotypes
    vecDbl popx = tst::readfile("population_x.dat");
    vecDbl ecotypes = tst::readfile("population_ecotype.dat");

    BOOST_CHECK_EQUAL(popx.size(), cumulsize);

    // In the last generation
    // Check that ecological trait values are higher in ecotype 1 than 0

    double xmax0 = -1.0;
    double xmin1 = 1.0;

    for (size_t i = lastgenfirst; i < popx.size(); ++i) {

        const size_t e = utl::dbl2size(ecotypes[i]);
        const double x = popx[i];

        if (e == 0u && x > xmax0) xmax0 = x;
        if (e == 1u && x < xmin1) xmin1 = x;

    }

    BOOST_CHECK(xmax0 < xmin1);

}
