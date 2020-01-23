#!/usr/bin/python3

# This is a script to change parameter values in the parameter file

# This script only performs string manipulation
# It will check that an argument is provided for any parameter to change
# But will not check its validity, nor will it check whether it is a number or 
# that enough values are supplied
# If the parameter to be matched is not present, will be added at the end
# Again without sanity check
# If the parameter file is invalid, this will be caught by the simulation crashing

import pytilities as pyt

if len(sys.argv) < 3:
	raise Exception("I need at least two arguments: a parameter name and a value")

# Read input arguments
args = sys.argv[1:]

# Set parameters
pyt.set_parameters(".", args)
