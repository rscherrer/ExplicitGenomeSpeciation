# ExplicitGenomeSpeciation GUI version

Branch|[![Travis CI logo](ci_setup/pics/TravisCI.png)](https://travis-ci.org)|[![Codecov logo](ci_setup/pics/Codecov.png)](https://www.codecov.io)
---|---|---
master|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=master)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=master)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/master)
develop|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=develop)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=develop)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/develop)
raph|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=raph)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=raph)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/raph)
thijs|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=thijs)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=thijs)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/thijs)
richel|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=richel)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=richel)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/richel)

This repository contains the source code of the simulation program ExplicitGenomeSpeciation.

# GUI

GUI Tabs:

![GUI overview](https://github.com/rscherrer/ExplicitGenomeSpeciation/blob/thijs/gui/Screenshot_2019-10-25_combined.png)

Extended coloration of histograms (Update 29-10-2019):
![GUI overview](https://github.com/rscherrer/ExplicitGenomeSpeciation/blob/thijs/gui/Screenshot_2019-10-29_11.37.21.png)




# Folder structure

The repository has the following folder structure (the content of each folder is explained):

The root folder contains files such as EGS.pro and EGS_test.pro, which are configuration files for building the program under Qt. The main.cpp file is also there.

## ci
   
Bash files needed by Travis to perform its different tasks, as defined in .travis.yml, during continuous integration to GitHub.

## library

The C++ header and source files needed to build the program in release mode.

## tests

The C++ header and source files needed to build the program in debug mode, with tests.

## gui

This folder contains extra code for the GUI.

## cluster

This folder contains scripts that make the program usable on the Peregrine cluster. Those scripts allow to launch large numbers of simulations on the server.

# Download the program

Download the repository by clicking the green button on the top right corner of the main GitHub page. Alternatively, you can download it from the command line using Git:

```{bash}
git clone https://github.com/rscherrer/ExplicitGenomeSpeciation
```

Then navigate to the repository to access the program.

# Build the program

You need to build the program from the source code. The project file `EGS.pro` contains all build instructions in release mode for Qt5, using qmake. The project file `EGS_test.pro` contains the build instructions for the debug version, which runs tests. The project file `EGS_gui` contains the instructions for building the GUI.

# Run the program

To run the program from the command line, use

```{bash}
./EGS
```

For non-default parameters you have to provide them in a parameter file, like this:

```{bash}
./EGS parameters.txt
```

# Default parameters

Refer to the manual of the program for this.

# Use on the cluster

For more information on how to build the program and use it to run multiple simulations on the Peregrine cluster, please refer to the README in folder /cluster.

# Downloads

 * [Windows executable](http://richelbilderbeek.nl/EGS_gui.zip)

