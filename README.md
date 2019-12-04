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

Where the executables are built. Any user can also build the program themselves from source using the IDE of their choice.

### debug
   
Where the debug configuration is built. This version of the program runs all the tests upon execution. Name of the executable: ./EGS_test. Check /EGS\_test.pro for details about compiler flags.

### release
 
Where the release configuration is built. This version of the program runs the actual simulation, optimized for speed. Name of the executable: ./EGS This file also contains a parameter file, ./parameters.txt that can be passed to the executable as a command line argument. Check /EGS.pro for details about compiler flags.

## gui

This folder contains Thijs's work on the GUI version of the program. I do not know yet what it contains and how to use it.

## cluster

This folder contains scripts that make the program usable on the Peregrine cluster. Those scripts allow to launch large numbers of simulations on the server.

# Run the program

From the command line:

```{bash}
./EGS
```

will launch a simulation with default parameters. For custom parameters, please provide a parameter file name as unique argument, for example:

```{bash}
./EGS parameters.txt
```

# Default parameters

Please refer to /MANUAL.md, the manuscript or the comments in the source code for a detailed explanation of what each parameter does. The default values are:

Name | Value
---|---
```rdynamics``` | 1
```trenewal``` | 0.001
```capacity``` | 100
```replenish``` | 1
```hsymmetry``` | 0
```ecosel``` | 1.8
```dispersal``` | 0.01
```birth``` | 4
```survival``` | 0.8
```sexsel``` | 10
```matingcost``` | 0.01
```maxfeed``` | 0.0004
```demesizes``` | 100, 0
```nvertices``` | 30, 30, 30
```nedges``` | 0, 0, 0
```nchrom``` | 3
```mutation``` | 0.001
```recombination``` | 3
```allfreq``` | 0.2
```scaleA``` | 1, 1, 1
```scaleD``` | 0, 0, 0
```scaleI``` | 0, 0, 0
```scaleE``` | 0, 0, 0
```skews``` | 1, 1, 1
```effectshape``` | 2
```effectscale``` | 1
```interactionshape``` | 5
```interactionscale``` | 1
```dominancevar``` | 1
```tburnin``` | 0
```tend``` | 10
```tsave``` | 10
```talkative``` | 1
```record``` | 1
```archsave``` | 0
```archload``` | 0
```archfile``` | architecture.txt
```seed``` | automatically generated
```ntrials``` | 100

# Parameter file

A parameter file should contain parameter names as they appear in the above example with default parameters, followed by their values. Parameters are separated from other parameters by any blank character e.g. new line, tab, or space. For example:

```
ecosel 1.0
hsymmetry 0.2
```

would set the ecological selection coefficient to 1 and habitat symmetry to 0.2, and the other parameters would keep their default values.

Note that each parameter expects a specific type of value, e.g. `ntrials` or `seed` expect integer numbers, while `archfile` expects a file name. 

Also note that some parameter expect multiple values. Parameters `nvertices`, `nedges`, `scaleA`, `scaleD`, `scaleI`, `scaleE` and `skews` all expect three values, one per trait, separated by any blank character. Parameter `demesizes` expects two values, one for each habitat.

# Notes

For more information on how to build the program and use it to run multiple simulations on the Peregrine cluster, please refer to the README in folder /cluster.

# Downloads

 * [Windows executable](http://richelbilderbeek.nl/EGS_gui.zip)

