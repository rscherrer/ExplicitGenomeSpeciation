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
    for (int t = 0; t < 10; ++t) {
        metapop.cycle(pars, arch);
        collector.analyze(metapop, pars, arch);
        collector.print(t, metapop);
        std::clog << "Printing time step " << t << '\n';
    }
    collector.shutdown();

    // Read output file
    vecDbl time = tst::readfile("time.dat");

}
