#include "library/Collector.h"
#include "library/Printer.h"
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
  pars.hsymmetry = 1.0;
  pars.demesizes = { 100u, 0u };
  pars.allfreq = 0.5;
  pars.ecosel = 1.0;
  pars.tburnin = 0u;

  GenArch arch = GenArch(pars);
  MetaPop metapop = MetaPop(pars, arch);
  Collector collector = Collector(arch);
  Printer printer = Printer();

  size_t cumulsize = 0u;
  size_t lastgenfirst = 0u;
  size_t lastgensize = 0u;

  for (int t = 0; t < 10; ++t) {

    metapop.cycle(pars, arch);
    collector.analyze(metapop, pars, arch);
    printer.print(t, collector, metapop);

    if (t != 9) lastgenfirst += metapop.getSize();
    else lastgensize = metapop.getSize();

    cumulsize += metapop.getSize();

  }

  assert(cumulsize - lastgenfirst == lastgensize);

  printer.shutdown();

    // Read output files
  std::vector<double> time = tst::readfile("time.dat");
  std::vector<double> ei = tst::readfile("EI.dat");
  std::vector<double> si = tst::readfile("SI.dat");
  std::vector<double> ri = tst::readfile("RI.dat");

    // Check output files

  BOOST_CHECK_EQUAL(time.size(), 10u);\
  BOOST_CHECK_EQUAL(ei.size(), 10u);
  BOOST_CHECK_EQUAL(si.size(), 10u);
  BOOST_CHECK_EQUAL(ri.size(), 10u);

    // Read individual trait values and ecotypes
  std::vector<double> popx = tst::readfile("individual_trait.dat");
  std::vector<double> ecotypes = tst::readfile("individual_ecotype.dat");

  BOOST_CHECK_EQUAL(popx.size() / 3u, cumulsize);

}
