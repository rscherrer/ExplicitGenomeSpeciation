#ifndef EXPLICITGENOMESPECIATION_TESTUTILITIES_H
#define EXPLICITGENOMESPECIATION_TESTUTILITIES_H

#include <fstream>
#include <vector>
#include <iostream>

namespace tst
{

  void makeValidParamFile();
  void makeValidParamFile2();
  void makeInvalidParamName();
  void makeInvalidParamValue();
  void makeInvalidParamValue2();
  void makeParamFileWithArchitecture();
  void makeParamFileWithMissingArchitecture();
  std::vector<double> readfile(const std::string&);
  std::vector<size_t> readfile2(const std::string&);

}

#endif
