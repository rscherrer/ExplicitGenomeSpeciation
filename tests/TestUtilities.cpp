#include "tests/TestUtilities.h"

std::vector<double> tst::readfile(const std::string &filename)
{
    // Open the input file
  std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);

    // Prepare storage for values
  double x;
  std::vector<double> v;

    // If the file is open
  if (file.is_open()) {

      // Loop through the file until we reach the end of the file
    while(file) {

        // Read elements
      file.read((char *) &x, sizeof(double));

        // Exit if reaching the end of the file
      if (!file.gcount()) break;

        // Store elements
      v.push_back(x);

    }
  }

    // Close the file
  file.close();

  return v;

}

std::vector<size_t> tst::readfile2(const std::string &filename)
{
    // Open the input file
  std::ifstream file(filename.c_str(), std::ios::in | std::ios::binary);

    // Prepare storage for values
  size_t x;
  std::vector<size_t> v;

    // If the file is open
  if (file.is_open()) {

      // Loop through the file until we reach the end of the file
    while(file) {

        // Read elements
      file.read((char *) &x, sizeof(size_t));

        // Exit if reaching the end of the file
      if (!file.gcount()) break;

        // Store elements
      v.push_back(x);

    }
  }

    // Close the file
  file.close();

  return v;

}

void tst::makeValidParamFile()
{
  std::ofstream file;
  file.open("validparamfile.txt");
  if (!file.is_open())
    std::cout << "Unable to open valid parameter test file.\n";

  file << "rdynamics" << '\t' << 0 << '\n'
  << "capacity" << '\t' << 10.0 << '\n'
  << "replenish" << '\t' << 1.0 << '\n'
  << "inflow" << '\t' << 400.0 << '\n'
  << "outflow" << '\t' << 100.0 << '\n'
  << "hsymmetry" << '\t' << 1.0 << '\n'
  << "ecosel" << '\t' << 1.0 << '\n'
  << "dispersal" << '\t' << 1.0E-3 << '\n'
  << "birth" << '\t' << 1.5 << '\n'
  << "survival" << '\t' << 0.8 << '\n'
  << "sexsel" << '\t' << 10.0 << '\n'
  << "matingcost" << '\t' << 0.01 << '\n'
  << "demesizes" << '\t' << 100 << '\t' << 100 << '\n'
  << "nvertices" << '\t' << 40 << '\t' << 20 << '\t' << 40 << '\n'
  << "nedges" << '\t' << 50 << '\t' << 50 << '\t' << 50 << '\n'
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

void tst::makeParamFileEmptyFreezer()
{
    std::ofstream file;
    file.open("paramfileemptyfreezer.txt");
    if (!file.is_open())
        std::cout << "Unable to open empty freezer parameter file.\n";

    file << "freezerfile" << '\t' << "empty_freezer.dat" << '\n'
         << "gensave" << '\t' << 1u << '\n'
         << "tend" << '\t' << 10u << '\n'
         << "tsave" << '\t' << 10u << '\n'
         << "tfreeze" << '\t' << 100u << '\n';

    file.close();
}

size_t tst::getFileSize(const std::string &filename) {

    std::ifstream input(filename, std::ios::binary);
    assert(input.is_open());
    input.seekg(0, std::ios::end);
    return static_cast<size_t>(input.tellg());

}

bool tst::isFileEmpty(const std::string &filename) {

    return !static_cast<bool>(getFileSize(filename));

}
