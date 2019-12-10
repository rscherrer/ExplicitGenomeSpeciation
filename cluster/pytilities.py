# This script contains the functions used in deploy.py

import re
import sys
import numpy as np
import os
from shutil import copyfile

# Update parameter values in a parameter file
def set_parameters(dir, args):

	# Loop through arguments
	# If argument starts with a dash, set it as the string to match
	# If not followed by a non-dashed argument, error
	# Concatenate all subsequent non-dashed arguments as values
	# Perform pattern replacement

	# Read input arguments
	nargs = len(args)

	idpars = [i for i, str in enumerate(args) if re.match(r'\-', str)]

	if len(idpars) == 0:
		raise Exception("Please provide parameter names preceded by dashes e.g. -mutation")

	# Read parameter file
	with open(dir + "/parameters.txt", "rt") as f:
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
			data = data + replacement # or add

	# Overwrite the parameter file
	with open(dir + "/parameters.txt", "wt") as f:
		f.write(data)

# Create folders with their corresponding parameter files for multiple combinations of parameters
def combine_parameters(dir, args, nreplicates):

	nargs = len(args)

	idpars = [i for i, str in enumerate(args) if re.match(r'\-', str) and str != "-v"]

	if len(idpars) == 0:
		raise Exception("Please provide parameter names preceded by dashes e.g. -mutation")

	paramentries = []

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
		values = [re.sub("_", " ", val) for val in values]
		entries = ["-" + parname + " " + val for val in values]
		paramentries.append(entries)

	# Produce multiple combinations of parameter entries
	fullentries = paramentries[0]
	i = 1
	while i < len(idpars):
		fullentries = [x + " " + y for x in fullentries for y in paramentries[i]]
		i = i + 1

	dirnames = []

	# For each combination
	for entry in fullentries:

		# Make a folder for each replicate
		dirtemplate = dir + "/sim_" + re.sub(" ", "_", re.sub("-", "", entry))

		for r in range(nreplicates):

			dirname = dirtemplate + "_r" + str(r + 1)
			if not os.path.exists(dirname): os.mkdir(dirname)	

			# Make a parameter file
			copyfile(dir + "/parameters.txt", dirname + "/parameters.txt")
			if os.path.isfile(dir + "/architecture.txt"):
				copyfile(dir + "/architecture.txt", dirname + "/architecture.txt")
			args = entry.split()
			set_parameters(dirname, args)			

			dirnames.append(dirname)

	return dirnames
