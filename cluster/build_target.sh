#!/bin/bash

# Make a target folder

if [ "$1" != "" ]; then
	if [ -d "$1" ]; then
		echo "Overwriting target folder"
		rm -rf $1		
	else
		echo "Creating target folder"		
	fi	
	mkdir $1

	# Produce an executable (from here cannot be tested locally)

	module load Qt5
	qmake ../EGS.pro # assuming we are in the cluster folder
	make --silent release

	# Move the executable to the target folder

	mv EGS $1

	# Copy the protocol and the running script to the target folder

	cp protocol.txt $1
	cp run_simulations.sh $1
	
else
	echo "Please provide target folder"
fi

