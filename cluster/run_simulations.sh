#!/bin/bash

# This script launches multiple simulations on the Peregrine cluster

# Run this script from within the target folder
# Pass it the protocol file as argument

if [ "$1" != "" ]; then

	if [ -f "$1" ]; then

		# Make simulation folders
		python3 deploy.py $1

		for folder in $(ls -d sim*)
		do

			# Pass the job file into the simulation folders
			cp job.sh $folder

			# Launch
			cd $folder
			sbatch job.sh
			cd ..

		done

	else

		echo "Invalid protocol file"

	fi

else

	echo "Please provide protocol file"

fi   
