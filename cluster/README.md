This file explains how to run the program of the Peregrine cluster of the University of Groningen.

# Prepare modules

For the simulations to be able to run on Peregrine, make sure some modules are made available, by adding them in the `.bash_profile` file located in the `/home/$USER` directory. Add to this file the following lines:

```{bash}
module add Python/3.6.4-foss-2018a
module add libpng/1.6.37-GCCcore-8.3.0
module add PCRE2/10.33-GCCcore-8.3.0
module add double-conversion/3.1.5-GCCcore-9.3.0
```

# Clone the repository

Navigate to where you want the repository to be. Then, type:

```{bash}
git clone https://github.com/rscherrer/ExplicitGenomeSpeciation
```

That will download the repository.

# Build the program

Assuming we are in the cluster folder (e.g. `cd ExplicitGenomeSpeciation/cluster`), build the target using:

```{bash}
./build_target.sh <path_to_target_folder>
```

If the target folder does not exist, it will be created. Otherwise, it will be overwritten. The target folder will contain the executable `EGS` as well as the necessary files to launch the simulations. All other files than the executable generated during the build (make files, object files) are stored into the `build` folder, within `cluster`.

# Update the protocol

Go to the target folder. A batch of simulations will be launched from here. Here, we consider that a batch of simulations is nothing more than an experiment aiming at exploring various parameter values, and that this experiment requires a protocol. The protocol defines what combinations of parameters are to be tested during the experiment. The target folder contains a file `protocol.txt` containing this information. The file should look something like this:

```
# This is a protocol file

#SBATCH --time=00:30:00
#SBATCH --mem=32Gb
#SBATCH --partition=gelifes

-mutation 0.001 0.01
-ecosel 0.1 0.2

N=1

```

All entries starting with `#SBATCH` are options to be written to the job file that will be passed to SLURM to launch the simulations. For more information about what SLURM options are available on Peregrine, please refer to the Peregrine wiki (https://redmine.hpc.rug.nl/redmine/projects/peregrine/wiki). All entries starting with a dash (e.g. `-mutation`) are interpreted as parameters to be tested. The dash should be immediately followed by a valid parameter name, and then by the values the parameter should take across simulations, separated by spaces. If several parameters with several values are provided, simulations will be launched for all combinations of values provided for these parameters. So, in the above example, simulations will be run for `mutation` = 0.001 -- `ecosel` = 0.1,
`mutation` = 0.001 -- `ecosel` = 0.2, `mutation` = 0.01 -- `ecosel` = 0.1 and `mutation` = 0.01 -- `ecosel` = 0.2. For parameters that take multiple values, such as `nvertices` that takes one value for each trait, please provide values of a single simulation separated by underscores. For example, `-nvertices 10_10_10 20_20_20` will launch simulations with 10 genes per trait as well as simulations with 20 genes per trait. Any line that does not start with `#SBATCH` or `-` will be ignored. Use a text editor (e.g. `nano` on the cluster) to make or edit a protocol file with the information needed to run your experiment. Parameters that are not mentioned in the protocol will take the values defined in the `parameters.txt` file, located in the target folder (a copy of it is in `cluster`). You can edit this file to set values to parameters that should be kept fixed across all simulations, thereby using the protocol only for those parameters that should vary. Any parameter that is neither in the fixed parameter file nor in the protocol will take default values. `N=5` refers to the number of replicate simulations to run per parameter combination. If no line is found starting with `N=` in the protocol, one replicate simulation will be run.

# Launch the simulations

From within the target folder, run:

```{bash}
./run_simulations.sh protocol.txt
```

This creates as many folders as there are simulations to be run, inside the target folder. Every simulation folder contains a `job.sh` file and a `parameters.txt` file. You can copy and paste extra files into the simulation folders using the same command, for example:

```{bash}
./run_simulations.sh protocol.txt whattosave.txt
``` 

will pass file "whattosave.txt" to each of the folders. This is important if you set `choosewhattosave` to 1 when running the simulation, but you can also use it e.g. to pass the same genetic architecture to all the simulation folders.

All job files are submitted to SLURM in one go, and should start running on the server. Use `squeue -u $USER` to check the status and progression of your jobs. The naming convention for the simulation folders is `sim_parameter1_value1_parameter2_value2_r<replicate_number>`, where `parameter1` and `parameter2` are the names of the parameters supplied in the protocol. For the example protocol above, the list of simulation folders should be:

```
sim_mutation_0.001_ecosel_0.1_r1
sim_mutation_0.001_ecosel_0.2_r1
sim_mutation_0.01_ecosel_0.1_r1
sim_mutation_0.01_ecosel_0.2_r1
```

All outputs to screen generated by a given simulation will be stored in a SLURM output file, with extension `.out`, in the simulation folder. You can visualize this file to diagnose the course of a simulation after it is done e.g. whether it went extinct, whether it ran at all or whether it completed succesfully.

# Note

Note that `run_simulations.sh` will run *all* the simulations present in the folder, which means that if you want to run new simulations in the same folder after having already run some, the simulations that already ran will run again and may be overwritten. If you do not want that happening, make sure to remove the old simulations or run your new simulations elsewhere.  

# Relaunch some simulations

Some simulations may have crashed or timed out and may have to be re-run. If the subset of simulations to redo is not small and rather specific, going through submitting a protocol again may be too cumbersome. The script `rerun_simulations.sh` is here for that purpose. Run it using:

```{bash}
./rerun_simulations.sh relaunch.txt
```

From the target folder, where `relaunch.txt` contains a list of names of simulations folders for which the simulation has to be relaunched. Before calling this script, you can update the `job.sh` file located in the target folder, and the new job file will be copied and pasted into each of the simulation folders to be resubmitted.

# Executable

This folder contains a version of the `EGS` executable that was compiled on the Peregrine cluster (add specifications here) using the steps described above. All the results presented in the accompanying study were generated with this version of the program. This is important for reproducibility of the results because stochastic simulations compiled with different compiler or library versions may give different pseudo-random numbers even with the same seed. 
