# ExplicitGenomeSpeciation

Branch|[![Travis CI logo](ci_setup/pics/TravisCI.png)](https://travis-ci.org)|[![Codecov logo](ci_setup/pics/Codecov.png)](https://www.codecov.io)
---|---|---
master|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=master)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=master)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/master)
develop|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=develop)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=develop)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/develop)
raph|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=raph)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=raph)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/raph)
thijs|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=thijs)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=thijs)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/thijs)
richel|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=richel)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=richel)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/richel)

This repository contains the source code of the simulation program ExplicitGenomeSpeciation.

# Folder structure

The repository has the following folder structure (the content of each folder is explained):

The root folder contains files such as EGS.pro and EGS_test.pro, which are configuration files for building the program under Qt.

## ci
   
Bash files needed by Travis to perform its different tasks, as defined in .travis.yml, during continuous integration to GitHub.

## library

The C++ header and source files needed to build the program in release mode.

## tests

The C++ header and source files needed to build the program in debug mode, with tests.

## build

Where the executables are built. Any user can also build the program themselves to other build foldeors using the IDE of their choice. This one is just made for my own purpose when I develop in QtCreator. I am not sure what platform these executables can be directly run from, and when the user should recompile the program themselves. The instructions for the debug and the release builds are in ../EGS_test.pro and ../EGS.pro, respectively.

### debug
   
Where the debug configuration is built. This version of the program runs all the tests upon execution, in debug mode. Name of the executable: ./EGS_test

### release
 
Where the release configuration is built. This version of the program runs the actual simulation, optimized for speed. Name of the executable: ./EGS This file also contains a parameter file, ./parameters.txt that can be passed to the executable as a command line argument

## gui

This folder contains Thijs's work on the GUI version of the program. I do not know yet what it contains and how to use it.

## cluster

This folder contains scripts that make the program usable on the Peregrine cluster. Those scripts allow to launch large numbers of simulations on the server.

# Notes

For more information on how to run the program once the executable is available and working, please refer to the README file in folder /build.

For more information on how to build the program and use it to run multiple simulations on the Peregrine cluster, please refer to the README in folder /cluster.

# Downloads

 * [Windows executable](http://richelbilderbeek.nl/EGS_gui.zip)

