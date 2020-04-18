#!/bin/bash

# This script launches multiple simulations on the Peregrine cluster

# Run this script from within the target folder
# Pass it the protocol file as argument
# Pass it files to copy into each simulation folder as extra optional arguments

if [ "$1" != "" ]; then

	if [ -f "$1" ]; then

		# Make simulation folders
		python3 deploy.py $1

		# For each simulation folder (they start with "sim")
		for folder in $(ls -d sim*)
		do

			# Pass the job file into the simulation folders
			cp job.sh $folder

			# Pass extra optional files into the simulation folders
			for i in "${@:2}" # loop through input arguments except the first one
			do 
				if[ -f "$i" ]; then
					cp "$i" $folder
				else
					echo "Invalid file to pass to the simulations"
				fi
			done

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
