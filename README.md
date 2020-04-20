# ExplicitGenomeSpeciation

Branch|[![Travis CI logo](ci_setup/pics/TravisCI.png)](https://travis-ci.org)|[![Codecov logo](ci_setup/pics/Codecov.png)](https://www.codecov.io)
---|---|---
master|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=master)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=master)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/master)
develop|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=develop)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=develop)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/develop)
raph|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=raph)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=raph)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/raph)
thijs|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=thijs)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=thijs)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/thijs)
richel|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=richel)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=richel)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/richel)

ExplicitGenomeSpeciation is a program that simulates the evolution of a species using individual-based modelling and incorporating explicit genomes and a flexible genotype-phenotype map.

# Example simulations




# Folder structure

The repository has the following folder structure:

* `ci`: Bash files needed by Travis for continuous integration to GitHub.
* `library`: headers and source code for the simulation
* `tests`: headers and source code for the tests
* `gui`: headers and source code for the graphical user interface (GUI)
* `cluster`: scripts to use the program on the Peregrine cluster

# Download the program

Download the repository by clicking the green button on the top right corner of the main GitHub page. Alternatively, you can download it from the command line using Git:

```{bash}
git clone https://github.com/rscherrer/ExplicitGenomeSpeciation
```

Then navigate to the repository to access the program.

# Build the program

The program was written in standard C++14 on Ubuntu 18.04 using QtCreator. The project file `EGS.pro` contains all build instructions that were used to compile the program. The project file `EGS_test.pro` contains the build instructions for the debug version, which runs tests and assertions. Tests were written using the Boost.Test library. The project file `EGS_gui` contains the instructions for building the GUI. Qt is needed to build either of these versions from the `.pro` files. Alternatively you can build the program from source with the compiler of your choice.

# Run the program

To run the program from the command line, use

```{bash}
./EGS
```

For non-default parameters you have to provide them in a parameter file, like this:

```{bash}
./EGS parameters.txt
```

# Parameters

The following table shows all the parameters of the program and their default values. All of them can be supplied by the user in a parameter file. The parameter file must contain on each line, the name of one parameter followed by a blank (e.g. space or tab) and the value(s) that this parameter must take. Parameter names must match those in the table. Parameters that are not in the parameter file will take default values. Parameters that take multiple values must be supplied as such, with values separated by blanks.

| Parameter | Meaning | Default |
|----|:----:|----:|
| rdynamics | Type of resource dynamics | 1 |
| capacity | Resource carrying capacity | 1.0 |
| replenish | Resource replenishment rate | 2375.0 |
| inflow | Resource inflow rate | 400.0 |
| outflow | Resource outflow rate | 100.0 |
| hsymmetry | Habitat symmetry | 0.0 |
| ecosel | Ecological divergent selection | 1.8 |
| dispersal | Dispersal rate | 1.0E-2 |
| birth | Birth rate | 1.0 |
| survival | Survival rate | 0.8 |
| sexsel | Sexual selection coefficient | 10.0 |
| matingcost | Cost of mate choice | 0.01 |
| demesizes | Initial population sizes | 100 0 |
| nvertices | Numbers of loci | 30 30 30 |
| nedges | Numbers of edges | 30 0 0 |
| nchrom | Number of chromosomes | 3 |
| mutation | Mutation rate | 1.0E-3 |
| recombination | Recombination rate | 3.0 |
| allfreq | Initial allele frequency | 0.2 |
| scaleA | Additive scaling parameter | 1.0 1.0 1.0 |
| scaleD | Dominance scaling parameter | 0.0 0.0 0.0 |
| scaleI | Interaction scaling parameter | 0.0 0.0 0.0 |
| scaleE | Environmental scaling parameter | 0.0 0.0 0.0 |
| skews | Gene network skewness | 1.0 1.0 1.0 |
| effectshape | Shape of the effect-size distribution | 2.0 |
| effectscale | Scale of the effect-size distribution | 1.0 |
| interactionshape | Shape of the interaction weight distribution | 5.0 |
| interactionscale | Scale of the interaction weight distribution | 1.0 |
| dominancevar | Scale of the dominance coefficient distribution | 1.0 |
| tburnin | Duration of the burn-in | 0 |
| tend | Duration of the simulation | 10 |
| tsave | Frequency of data recording | 10 |
| tfreeze | Frequency of individual whole genome sequencing | 100 |
| talkative | Verbose during the simulation | 1 |
| record | Whether to record data | 1 |
| datsave | Whether to save data | 1 |
| choosewhattosave | Whether to choose what variables to save | 0 |
| gensave | Whether to save individual whole genomes | 0 |
| archsave | Whether to save genetic architecture | 0 |
| archload | Whether to load architecture from a file | 0 |
| parsave | Whether to save parameters | 1 |
| archfile | Architecture file | architecture.txt |
| parfile | Parameter log file | paramlog.txt |
| orderfile | File listing variables to save | whattosave.txt |
| seed | Random seed | randomly generated |  
| ntrials | Number of mating trials | 100 |

Useful information about the parameters:

* rdynamics: chemostat (1) or logistic (0) resource dynamics

* capacity and replenish: parameters of the logistic resource dynamics, only used if rdynamics is 0

* inflow and outflow: parameters of the chemostat resource dynamics, only used if rdynamics is 1

* nvertices: must be at least 2 for each trait

* nedges: cannot be higher than n (n - 1) / 2 (complete graph) if n is the number of vertice for a given trait. In practice, the preferential attachment algorithm may fail to produce a network if nedges is very high relative to the number of vertices. Because this is a stochastic process, it is not possible to provide a clear cut-off.

* recombination: this is not a proportion, but the rate of an exponential distribution describing the distance between successive crossover points. It is more-or-less equivalent to the number of crossovers per genome per meiosis event.

* allfreq: the frequency of the "1" allele in the initial population, as opposed to the "0" allele. There are only two possible alleles at any given locus.

* scale parameters: those scale the effects of different components of the genotype-phenotype map. They take one value per trait. Make sure that the scaling coefficients sum up to 1 for each trait, otherwise different traits or simulations may not be comparable.

* skew: a value of 1 produces a power law degree distribution, less than 1 produces a stretched exponential and more than 1 produces a star-shaped network

* effect-size, dominance and interaction weight scale parameters: if the scale or variance of these distributions is 0, all effect-sizes or weights take the value 1.

* tburnin, tend, tsave and tfreeze all are in numbers of generations

* tfreeze determines when to save individual whole genomes. This uses a separate saving frequency from the rest of the data (which use tsave) because it is very space-consuming

* talkative will print to the prompt the generation if set to 1

* if record is 1, the collector module will analyze the data every tsave step. The data will be saved to files only if datsave is 1. All possible variables will be saved to files if hoosewhattosave is 0, otherwise the program will look into the file provided in orderfile for a list of the variables to save.

* gensave needs to be 1 in order to save whole individual genomes (every tfreeze generations)

* setting parsave to 1 will save all the parameters used to run the simulation to the file defined in parfile. This file may be used e.g. to retrieve the exact seed used in a simulation

* archfile is a file containing the genetic architecture. If archsave is 1, the architecture is saved into this file, while if archload is 1, the architecture is loaded from this file.

* The random seed is generated by default using the clock, but you can set your own integer instead

* ntrials is the number of trials used to compute the degree of reproductive isolation

# Architecture file

An architecture file must have a specific organization to be loaded successfully. It is a text file consisting of a series of fields, each directly followed by a line break and a series of numbers to be loaded into that field, separated by blanks. The fields are:

* chromosomes: the location of the end of each chromosome between the start and end of the genome. Make sure that the number of values equals nchrom of the simulation you are going to launch.

* traits: the trait encoded by each locus

* locations: the location of each locus along the genome

* effects: the effect-size of each locus

* dominances: the dominance coefficient of each locus

Make sure that the locus-specific fields have as many values as there are vertices in the simulation you want to run.

Each gene network has its own field named "network" followed by the trait encoded by the network (0, 1 or 2) and the number of edges in that network, all separated by blanks. Those network headers are followed by a line break and several subfields, each followed by their series of blank-separated values.

* weights: the interaction weight of each edge

* edge0: the first interaction partner locus for each edge

* edge1: the second interaction partner locus for each edge

Fields and subfields are separated from each other by line breaks.

You will rarely have to make your own architecture file manually. The most common use case is to run a simulation with the archsave flag, save the genetic architecture in the file, and then re-use this architecture in another simulation with the archload flag and the same file.

# Order file

This optional file contains a list of names of variables to save, if record, datsave and choosewhattosave are all 1. Each variable is saved as a series of 64bit double precision floating point numbers in a binary file with the .dat extension (and so needs to be decoded to be read). This is to optimize space use and read / write speed. There is no saving during the burn-in period. The possible variables are:

| Variable | Meaning | Values per timepoint |
|----|:----:|----:|
| time | Generation | 1 |
| population_size | Number of individuals | 1 |
| ecotype_size | Ecotype population size | 1 per ecotype |
| resources | Equilibrium resource concentrations | 1 per habitat per resource |
| means | Mean trait values | 1 per trait |
| ecotype_means | Ecotype-mean trait values | 1 per trait per ecotype |
| varP | Phenotypic variance | 1 per trait |
| varG | Genetic variance | 1 per trait |
| varA | Additive variance | 1 per trait |
| varD | Dominance variance | 1 per trait |
| varI | Interaction variance | 1 per trait |
| varN | Non-additive variance | 1 per trait |
| varT | Variance in allele frequencies | 1 per trait |
| Pst | Phenotypic differentiation | 1 per trait |
| Gst | Genetic effect differentiation | 1 per trait |
| Qst | Additive effect differentiation | 1 per trait |
| Cst | Non-additive effect differentiation | 1 per trait |
| Fst | Genotypic differentiation | 1 per trait |
| EI | Ecological isolation | 1 |
| SI | Spatial isolation | 1 |
| RI | Reproductive isolation | 1 |
| genome_varP | Locus-specific phenotypic variance | 1 per locus |
| genome_varG | Locus-specific genetic variance | 1 per locus |
| genome_varA | Locus-specific additive variance | 1 per locus |
| genome_varD | Locus-specific dominance variance | 1 per locus |
| genome_varI | Locus-specific interaction variance | 1 per locus |
| genome_varN | Locus-specific non-additive variance | 1 per locus |
| genome_Pst | Locus-specific Pst | 1 per locus |
| genome_Gst | Locus-specific Gst | 1 per locus |
| genome_Qst | Locus-specific Qst | 1 per locus |
| genome_Cst | Locus-specific Cst | 1 per locus |
| genome_Fst | Locus-specific Fst | 1 per locus |
| genome_alpha | Locus-specific average mutational effect | 1 per locus |
| genome_meang | Locus-specific mean genetic value | 1 per locus |
| genome_freq | Locus-specific allele frequency | 1 per locus |
| network_corgen | Edge-specific correlation in genetic values | 1 per edge |
| network_corbreed | Edge-specific correlation in breeding values | 1 per edge |
| network_corfreq | Edge-specific correlation in allele frequencies | 1 per edge |
| network_avgi | First-partner edge-specific epistatic variance in average effect | 1 per edge |
| network_avgj | Second-partner edge-specific epistatic variance in average effect | 1 per edge |
| individual_ecotype | Individual ecotype | 1 per individual |
| individual_habitat | Individual habitat | 1 per individual |
| individual_trait | Individual trait values | 1 per individual per trait |
| individual_midparent | Individual midparent trait values | 1 per individual per trait |

(Edges are ordered by trait.)

Note that if you are choosing what variables to save, it is important to save "time" because the functions provided in the accompanying R package for analyses assume that this file is present. Also note that if you are going to save individual data, you may want to save "population_size" for the same reason.

# Analysis with egssimtools

The program comes with the R package `egssimtools`, which contains useful functions to read and process the data generated during the simulations. It is particularly handy to convert the binary data files saved by the program into workable data sets and vectors. The package comes with a vignette explaining how to use it.

# Use on the cluster

Check out the README in folder `cluster`.

# GUI

GUI Tabs:

![GUI overview](https://github.com/rscherrer/ExplicitGenomeSpeciation/blob/thijs/gui/Screenshot_2019-10-25_combined.png)

Extended coloration of histograms (Update 29-10-2019):
![GUI overview](https://github.com/rscherrer/ExplicitGenomeSpeciation/blob/thijs/gui/Screenshot_2019-10-29_11.37.21.png)

# Downloads

 * [Windows executable](http://richelbilderbeek.nl/EGS_gui.zip)

