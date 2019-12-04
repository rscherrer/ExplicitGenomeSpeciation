#!/bin/bash

# Run this script from within the target folder
# Pass it the protocol file as argument

if [ "$1" != "" ]; then

	if [ -d "$1" ]; then

		# Make simulation folders and launch
		python3 deploy.py $1

	else

		echo "Invalid protocol file"

	fi

else

	echo "Please provide protocol file"

fi   
