# ExplicitGenomeSpeciation

Branch|[![Travis CI logo](ci_setup/pics/TravisCI.png)](https://travis-ci.org)|[![Codecov logo](ci_setup/pics/Codecov.png)](https://www.codecov.io)
---|---|---
master|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=master)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=master)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/master)
develop|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=develop)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=develop)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/develop)
raph|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=raph)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=raph)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/raph)
thijs|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=thijs)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=thijs)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/thijs)
richel|[![Build Status](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation.svg?branch=richel)](https://travis-ci.org/rscherrer/ExplicitGenomeSpeciation)|[![codecov.io](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/coverage.svg?branch=richel)](https://codecov.io/github/rscherrer/ExplicitGenomeSpeciation/branch/richel)


# Run a simulation locally

## Step 1 -- Download the repository

You can download the repository from GitHub, using the green Download tab, or do it from the command line:

```bash
git clone https://github.com/rscherrer/ExplicitGenomeSpeciation
```

## Step 2 -- Update the parameter file

The file `parameters.txt` contains parameters and their values to be passed to the simulation program. The program reads this file. You can change the parameter values to the values you want. Beware that some parameters take multiple values as input. If some parameter names are missing from the parameter file, they will take default values (set from within the program). If parameter names or values are invalid, the program will crash.

## Step 3 -- Launch the simulation

To launch a simulation with the parameter values set in the parameter file, run:

```bash
./EGS parameters.txt
```

It is also possible to run the program with default arguments, in which case the parameter file is not needed:

```bash
./EGS
```

The simulation should print messages as well as generations on the screen as the simulation proceeds. If parameter `record` is set to `1`, the program should write data into multiple output files in the current directory.


# Simulate in batch on the Pergrine cluster

How to run simulations in batches on the cluster.

## Step 1 -- Clone the repository

On your session on the Peregrine cluster:

```bash
git clone https://github.com/rscherrer/ExplicitGenomeSpeciation
```

Then move into the repository:

```bash
cd ExplicitGenomeSpeciation
```

## Step 2 -- Compile the program

To produce the executable use:

```bash
./compile.sh
```

which is tailored to the tools available on the Peregrine cluster. Here we are using Qt5.

## Step 3 -- Write a protocol

Launching multiple simulations to explore the behavior of a model is essentially doing an experiment, and follows a protocol. Write a protocol text file in the current directory, which should look something like:

```
#SBATCH --time=00:05:00
#SBATCH --mem=4Gb
#SBATCH --partition=regular

-mutation 0.1 0.01
-nvertices 50_50_50 100_50_50
```

where each line starting with `#SBATCH` is a SLURM option to be written in the job file that will be submitted to the cluster to run, and each line starting with `-<parameter_name>` is an option to be passed to the program that will set up all combinations of parameters to simulate. The syntax to follow for the SLURM options is the regular syntax used for job files on Peregrine (check the Peregrine wiki for info on that). The syntax for the parameters to test is:

```
-parameter1 v11 v12
-parameter2 v21 v22
```
In this case, all combinations of values for `parameter1` and `parameter2` will be tested. Parameter names should be written exactly how they are expected to be read from the `parameters.txt` file by the simulation program. The program will crash if the parameter names or values are invalid. In the special case of parameters that take multiple values, e.g. the number of vertices underlying each trait in the simulation, the different values for a given parameter set must be separated by underscores, and spaces are reserved to separate different combinations of parameters. For example:

```
-mutation 0.01 0.1
-nvertices 10_10_10 20_20_20
```

will produce four parameter combinations: (1) mutation rate of 0.01 and 10 loci for each trait, (2) mutation rate of 0.01 and 20 loci for each trait, (3) mutation rate of 0.1 and 10 loci for each trait, and (4) mutation rate of 0.1 and 20 loci for each trait.

## Step 4 -- Launch the simulations in batch

A script called `deploy.py` is used to launch simulations for all the combinations defined in the protocol. To make Python 3.6 available as well as other modules needed for the simulations to work on the cluster, first run:

```bash
./load.sh
```

Then, launch the simulations by running:

```bash
python3 deploy.py protocol.txt
```

where `protocol.txt` is the name of the protocol file. The script `deploy.py` will read the protocol file, create a folder in the current directory for each parameter combination to be run, each one with its own updated parameter file, and finally submit each simulation as a job to SLURM from within its own simulation folder. The simulation folders are named after their parameter sets. For example, the folder containing a simulation with mutation set to 0.01 and 10 loci per trait will be named `sim_-mutation_0.01_-nvertices_10_10_10`. The script `deploy.py` creates parameter sets by first creating a template parameter file with default values, then modifying the ones provided in the protocol. You can have a look at what are the default parameter values by looking at the `parameter.txt` present within each simulation folder, or by looking at the Python function `reset_parameters` stored in the file `pytilities.py`. Parameter that do not appear in the protocol stick to their default values and are not involved in the naming of the simulation folders. As stated above, all SLURM options must be specified in the protocol file. Besides these, the job file being submitted only consists of:

```bash
../EGS parameters.txt
```

Here, the executable `EGS` is by default expected to be located in the parent directory (`..`), but the path (relative to within the simulation folder) can be provided as second argument to `deploy.py`. For example:

```bash
python3 deploy.py protocol.txt .
```

will expect the `EGS` executable to be located within the simulation folder itself (`.`).

Running `deploy.py` should be done **with caution**, as it can potentially submit a large number of jobs to SLURM. So, pay attention that the protocol is right before launching many simulations that may significantly decrease your priority score.

## Step 5 -- Check the status of the simulations

You can use the command:

```bash
./check.sh
```

to have a look at the status of the simulations in SLURM and see which ones are done, pending, cancelled etc. This is just a less tedious way to call the bash command `squeue -u $USER`.

The root folder of the repository should now contain multiple folders, one for each simulation and named after the parameters that make them deviate from the default settings. Each simulation folder should contain (1) a parameter file (`parameter.txt`), (2) a job file (`job.sh`), (3) a multitude of data files output by the simulation (if parameter `record` is set to `1`), and (4) a SLURM output file (extension `.out`) summarizing the state of the job and what has been printed to the screen.

## Note

It is possible to launch multiple simulations outside the root folder of the repository. All what is needed in the current working directory are the Python scripts `deploy.py` and `pytilities.py`, and the protocol file. The script `deploy.py` just needs the location of the `EGS` application (which becomes available only after compiling, e.g. with `compile.sh` from the root directory) as second argument to be able to run.

# Quick look at the simulation results

The root directory comes with a few Python scripts that allow visualizing the output of a given simulation. Here we assume that a simulation as run and completed, and we refer to the "simulation folder" as where the data files (in `.dat` format) have been saved. *Note*: the data files are encoded in binary, so they cannot be visualized using a text editor.

The visualization scripts are:

* `read.py`, to print to the screen a series of values of a population-level variable through time
* `plot.py`, to plot lines showing population-level variables through time
* `hist.py`, to plot the distribution of an individual-level variable across the population at the last saved timepoint
* `heatmap.py`, to plot a heatmap showing the distribution of an individual-level variable across the population through time

## Step 1 -- Prepare Python

First you need to move the Python scripts needed for visualization into the simulation folder. Then, some Python libraries might be needed and not installed, so make sure that all dependencies are available by running:

```bash
./pystall.sh
```

This should work from anywhere on the computer.

## Step 2 -- Run the scripts

Either run the Python scripts by calling the Python 3.x interpreter explicitly:

```bash
python3 <name_of_script>.py
```

or 

```
./<name_of_script>.py
```

since all scripts contain a shebang in their first line stating the location of the interpreter as `usr/bin/python3`.

All scripts can take command line arguments, which must be data file names present in the working directory. Except for `read.py`, all scripts can run without argument and will then read default data files. By default, `plot.py` plots the ecological isolation statistic through time, `hist.py` plots a histogram of the ecological trait values in the population, split by ecotype, at the last saving point, and `heatmap.py` plots a heatmap of the distribution of the ecological trait values in the population through time. Data file names are supplied as in the following examples:

```bash
python3 plot.py SI.dat
```

for the spatial isolation statistic through time, or

```bash
python3 heatmap.py population_y.dat
```

for a distribution of mating preferences through time.

`plot.py` can plot more than one variable through time on the same plot. It does so if provided with multiple arguments. For exampe:

```bash
python3 plot.py EI.dat SI.dat RI.dat
```

will plot the ecological, spatial and mating isolation statistics on the same graph.

Note that since `heatmap.py` and `hist.py` plot distributions of individual-level variables, they must be provided with relevant dataset. Those have the prefix `population_` in their file name.

# Downloads

 * [Windows executable](http://richelbilderbeek.nl/EGS_gui.zip)
