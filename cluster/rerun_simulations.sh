#!/bin/bash

# This script re-launches a set of simulations
# It needs to be provided with a text file containing the names of the folders to re-run

# Run this script from within the target folder

# Arguments:
# File containing a list of simulation folders to relaunch

if [ "$1" != "" ]; then

	if [ -f "$1" ]; then

		# Read the file containing the list of simulations to relaunch
		input=$1
		while IFS= read -r line
		do

			folder=$line

			# Copy and paste the (updated) job file into each folder
			cp job.sh $folder

			# Submit to SLURM
			cd $folder
			sbatch job.sh
			cd ..	

		done < "$input"		

	else

		echo "Invalid simulation list file"

	fi

else

	echo "Please provide a file with simulations to relaunch"

fi
