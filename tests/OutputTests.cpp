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
    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    Collector collector = Collector(arch);

    vecDbl time(10);
    vecDbl varpx(10);
    size_t cumulsize = 0u;

    for (int t = 0; t < 10; ++t) {

        time[t] = utl::size2dbl(t);
        metapop.cycle(pars, arch);
        collector.analyze(metapop, pars, arch);
        collector.print(t, metapop);

        varpx[t] = collector.getVarP(0u);

        cumulsize += metapop.getSize();

    }
    collector.shutdown();

    // Read output file
    vecDbl rtime = tst::readfile("time.dat");
    vecDbl rvarpx = tst::readfile("varP_x.dat");

    BOOST_CHECK_EQUAL(time.size(), rtime.size());
    BOOST_CHECK_EQUAL(varpx.size(), rvarpx.size());

    // Read individual trait values
    vecDbl rpopx = tst::readfile("population_x.dat");

    BOOST_CHECK_EQUAL(rpopx.size(), cumulsize); // passes just fine

}
