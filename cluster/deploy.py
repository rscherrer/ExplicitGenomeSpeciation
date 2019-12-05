#!/usr/bin/python3

import pytilities as pyt
import sys
import re

# This program should read a protocol file (passed as single argument)
# Combine parameter values provided and make a folder for each of them (use combine_parameters)

# There should be at least one argument
if len(sys.argv) == 1:
	raise Exception("Path to protocol file must be provided.")

if len(sys.argv) > 2:
	raise Exception("Too many arguments provided.")

# The first argument is the path to the protocol file
filename = sys.argv[1]

joboptions = ''
paramstring = ''
nreplicates = 1

# Read protocol
with open(filename, "rt") as f:

	# Go through lines
	for cnt, line in enumerate(f):

		# Is there a SBATCH argument?
		if re.match(r'#SBATCH', line):
			joboptions = joboptions + line

		# Is there a parameter argument?
		if re.match(r'\-', line):
			paramstring = paramstring + re.sub('\n', '', line) + ' '

		# Is there a number of replicates?
		if re.match(r'N=', line):
			nreplicates = int(re.sub('N=', '', line))

paramstring = re.sub(' $', '', paramstring)

# Make a template job file
job = job = "#!/bin/bash\n"
job = job + joboptions
job = job + "./EGS parameters.txt\n"

# Make all parameter combinations

params = paramstring.split(' ')

dirnames = pyt.combine_parameters(".", params, nreplicates)

# Save the job file in the target folder
with open("job.sh", "w+") as f:
	f.write(job)