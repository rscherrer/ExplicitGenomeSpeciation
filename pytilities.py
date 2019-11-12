import re
import sys
import numpy as np
import os

def reset_parameters(dir):

	# Create file if non extisting
	if not os.path.isfile(dir + "/parameters.txt"):
		f = open(dir + "/parameters.txt", "w+")
		f.close()

	# Read parameter file
	with open(dir + "/parameters.txt", "rt") as f:
		data = f.read()	

	data = "rdynamics 1\n"
	data = data + "trenewal 0.001\n"
	data = data + "capacity 100\n"
	data = data + "replenish 1\n"
	data = data + "hsymmetry 1\n"
	data = data + "ecosel 1.4\n"
	data = data + "dispersal 0.001\n"
	data = data + "birth 2\n"
	data = data + "survival 0.8\n"
	data = data + "sexsel 10\n"
	data = data + "matingcost 0.01\n"
	data = data + "maxfeed 0.0004\n"
	data = data + "demesizes 100 0\n"
	data = data + "nvertices 50 50 50\n"
	data = data + "nedges 0 0 0\n"
	data = data + "nchrom 3\n"
	data = data + "mutation 0.0001\n"
	data = data + "recombination 3\n"
	data = data + "allfreq 0.2\n"
	data = data + "scaleA 1 1 1\n"
	data = data + "scaleD 0 0 0\n"
	data = data + "scaleI 0 0 0\n"
	data = data + "scaleE 0 0 0\n"
	data = data + "skews 1 1 1\n"
	data = data + "effectshape 2\n"
	data = data + "effectscale 1\n"
	data = data + "interactionshape 5\n"
	data = data + "interactionscale 1\n"
	data = data + "dominancevar 1\n"
	data = data + "tburnin 1000\n"
	data = data + "tend 3000\n"
	data = data + "tsave 10\n"
	data = data + "record 1\n"
	data = data + "archsave 0\n"
	data = data + "ntrials 100\n"

	# Overwrite parameter file
	with open(dir + "/parameters.txt", "wt") as f:
		f.write(data)	

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


def combine_parameters(dir, args):

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

		# Make a folder if not already there
		dirname = dir + "/sim_" + re.sub(" ", "_", entry)
		if not os.path.exists(dirname): os.mkdir(dirname)

		# Make a parameter file
		reset_parameters(dirname)
		args = entry.split()
		set_parameters(dirname, args)

		dirnames.append(dirname)

	return dirnames


