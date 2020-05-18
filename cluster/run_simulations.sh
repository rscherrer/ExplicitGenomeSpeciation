#!/bin/bash

# This script launches multiple simulations on the Peregrine cluster

# Run this script from within the target folder
# Pass it the protocol file as argument
# Pass it files to copy into each simulation folder as extra optional arguments

# First check that the optional extra files are present
for i in "${@:2}"
do
	if [ ! -f $i ]; then
		echo "Warning: file $i not found"
    fi
done

# If a protocol file is provided
if [ "$1" != "" ]; then

	# Is it found?
	if [ -f "$1" ]; then

		# Make simulation folders
		python3 deploy.py $1

		# For each simulation folder (they start with "sim")
		for folder in $(ls -d sim*)
		do

			# Pass the job file into the simulation folders
			cp job.sh $folder

			# Pass extra optional files into the simulation folders
			for i in "${@:2}"
			do
				if [ -f $i ]; then
					cp $i $folder
	    		fi
			done

			# launches
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
