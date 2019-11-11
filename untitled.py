#!/usr/bin/python3

import pytilities as pyt

# Provide a simulation protocol to submit to peregrine
# The protocol should contain the different parameter values to combine
# As well as slurm flags for the cluster
# From these this script will create folders for all combinations of parameters
# Write a job script for each one of them with the requested options
# And launch the simulations

# Read protocol
args = ["-mutation", "0.01", "0.1", "-ecosel", "0", "1"]

# Make all parameter combinations
dirnames = pyt.combine_parameters(".", args)

# Make a job file
joboptions = ["#SBATCH --partition=gelifes"]
direxec = ".." # assume exectuable is in parent directory
job = "#!/bin/bash\n"
job = job + "#SBATCH --time=00:30:00\n"
job = job + "#SBATCH --mem=32Gb\n"
for opt in joboptions:
	job = job + opt + '\n'
job = job + direxec + "/EGS parameters.txt\n"

# Copy the job file in all directory and submit to SLURM
for dir in dirnames:

	print(dir)

	f = open(dir + "/job.sh", "w+")
	f.write(job)
	f.close()