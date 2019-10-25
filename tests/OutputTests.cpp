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

    vecDbl time(10);

    Param pars;
    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    Collector collector = Collector(arch);

    for (int t = 0; t < 10; ++t) {

        time[t] = utl::size2dbl(t);

        metapop.cycle(pars, arch);
        collector.analyze(metapop, pars, arch);
        collector.print(t, metapop);
    }
    collector.shutdown();

    // Read output file
    vecDbl rtime = tst::readfile("time.dat");

    BOOST_CHECK_EQUAL(time.size(), rtime.size());

    // There are 80 bytes in the output file
    // So 10 values of 8 bytes each
    // It's the reading function that reads 11 values, duplicating the last one
    // Why is it doing so?
    // Is the reader in Python doing the same thing?

    // On the C++ side it probably comes from ifstream::eof()
    // It seems that the end of file is not reached at the last element
    // So the loop goes one step further, and for some reason reads the
    // last element again

}
