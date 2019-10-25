#include "library/Collector.h"
#include "library/Random.h"
#include "library/Utilities.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <cassert>

// Test the output files

// One big test case in which we make sure everything is saved properly

BOOST_AUTO_TEST_CASE(OutputFilesAreCorrectlyWritten)
{
    std::clog << "Testing output files...\n";
    Param pars;
    GenArch arch = GenArch(pars);
    MetaPop metapop = MetaPop(pars, arch);
    metapop.cycle(pars, arch);
    Collector collector = Collector(arch);
    collector.analyze(metapop, pars, arch);

}
