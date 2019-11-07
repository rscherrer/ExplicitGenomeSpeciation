#!/usr/bin/python3

# This is a script to change parameter values in the parameter file

# This script only performs string manipulation
# It will check that an argument is provided for any parameter to change
# But will not check its validity, nor will it check whether it is a number or 
# that enough values are supplied
# If the parameter to be matched is not present, will be added at the end
# Again without sanity check
# If the parameter file is invalid, this will be caught by the simulation crashing

import re
import sys
import numpy as np

# Loop through arguments
# If argument starts with a dash, set it as the string to match
# If not followed by a non-dashed argument, error
# Concatenate all subsequent non-dashed arguments as values
# Perform pattern replacement

if len(sys.argv) < 3:
	raise Exception("I need at least two arguments: a parameter name and a value")

# Read input arguments
args = sys.argv[1:]
nargs = len(args)

idpars = [i for i, str in enumerate(args) if re.match(r'\-', str)]

if len(idpars) == 0:
	raise Exception("Please provide parameter names preceded by dashes e.g. -mutation")

# Read parameter file
with open("parameters.txt", "rt") as f:
	data = f.read()	

# For each parameter
for i in range(len(idpars)):

	# Define parameter name and value(s)
	idpar = idpars[i] # parameter index in the list of arguments
	parname = re.sub(r'\-', '', args[idpar]) # parameter name
	idnext = nargs
	if i != len(idpars) - 1: idnext = idpars[i + 1] 
	nvalues = idnext - idpar - 1
	if nvalues == 0:
		raise Exception("No value supplied for parameter " + parname)
	start = idpar + 1
	stop = idpar + nvalues + 1
	values = args[start:stop] # parameter values
	sep = ' '
	values = sep.join(values)

	# Search and replace
	pattern = parname + ".*\n"
	replacement = parname + sep + values + '\n'
	if re.search(pattern, data):		
		data = re.sub(pattern, replacement, data)
	else:
		data = data + replacement

# Overwrite the parameter file
with open("parameters.txt", "wt") as f:
	f.write(data)




