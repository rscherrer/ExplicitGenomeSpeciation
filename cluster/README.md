This file explains how to run the program of the Peregrine cluster of the University of Groningen.

# Clone the repository

Enter the cluster. Navigate to where you want the repository to be. Then, type:

```{bash}
git clone https://github.com/rscherrer/ExplicitGenomeSpeciation
```

That will download the repository.

# Build the program

Assuming we are in the repository root (e.g. `cd ExplicitGenomeSpeciation`), build the target using:

```{bash}
./cluster/build_target.sh <path_to_target_folder>
```

I am not sure whether the path to the target folder should be relative to /cluster or to /. I am not sure either whether the shell script should be run from / or from /cluster, and if that changes what the relative path should be. If the target folder does not exist, it will be created. If it exists, it will be overwritten.

The target folder will contain the executable, a protocol.txt file and a run_experiment.sh file.

# Update the protocol

Assuming we are in the target folder, on the cluster, use:

```{bash}
nano protocol.txt
```

To be able to edit the protocol. Here there should a description of what the protocol is and what it should contain. There could be a header in that file, commented using #, also explaining what to write into it.

A protocol file is a file that is read by a shell script in the next step. The information contained in the protocol allows to setup the different jobs that need to be submitted to SLURM. The protocol can contain the following information:

1) SLURM options for running the jobs. Thos include run time, memory, partition etc. Please refer to SLURM or Peregrine documentation for what kind of options can be passed.
The SLURM options should appear in the protocol exactly as they would appear on the job file submitted to SLURM through the sbatch command.
Each SLURM option should be on a separate line in the protocol and start with #SBATCH.

e.g. #SBATCH --time=00:10:00

2) Parameter settings. Those define what values of what parameters should be run.
Parameter names should be followed by their value to be tested, separated by a space.
Each parameter must be on a separate line in the protocol and start with the name of the parameter.
If several values of a parameter are to be tested, add those values at the end of the line of that parameter, separating them from other values using spaces.
Same as explained in the /build/README.md file on how to run a single simulation, the job will crash if some parameter names are invalid or if the values are invalid.
If several parameters are provided to the protocol, meaning that several parameters must be tested,
all combinations of all provided values of the parameters appearing in the protocol will be tested.

For example,

mutation 0.001 0.01
ecosel 0.1 0.2

will launch four simulations, with all combinations of values of mutation rate and ecological selection coefficient provided.

Some parameter are supposed to take multiple values, e.g. nvertices, which is the number of loci underlying each trait in the simulation (3 values).
To enter a value of a multiple-valued parameter in the protocol, separate the different values of the same trial by underscores, as here spaces are reserved for separating parameter combinations.

For example,

nvertices 10_10_10 20_20_20

will launch two simulations, one with 10 loci for each trait, and another one with 20 loci for each trait.

Any line that does not start with either #SBATCH or a valid parameter name will be ignored.

Note that the parameters supplied in the protocol are only those parameters to change from the default parameter setting. You can have a look at the default parameter setting in /default.txt.
All parameters not supplied will therefore stay identical to their default values.
The default parameter values used are in file /cluster/default.txt
You change their values if you want identical values for all combinations you want to test.

# Launch the experiment

```{bash}
./run_experiment.sh
```

This creates as many folders as there are simulations to be run. Every folder contains a job.sh file and a parameters.txt file. All job files are submitted to SLURM and should run on the server. 
As the simulations progress, the simulation folders should fill up with .dat output files. The messages normally output to the screen by the program during the simulation will not appear in the terminal, but instead be printed to .out files, one per simulation folder, by SLURM. You can read those .out files to diagnose the course of a given simulation.

The naming convention for the simulation folders is 

```{bash}
sim-parameter1-value1-parameter2-value2
```

where parameter<i> is the ith parameter provided in the protocol, and
value<i> is the value of the ith parameter for this simulation.
Only parameters that change from the default due to the protocol appear 
in the name of the simulation folder.
