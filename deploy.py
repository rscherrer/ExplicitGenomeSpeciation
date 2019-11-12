#!/usr/bin/python3

import pytilities as pyt
import sys
import re
import subprocess

# This program should read a protocol file (passed as single argument)
# Combine parameter values provided and make a folder for each of them (use combine_parameters)

# There should be at least one argument
if len(sys.argv) == 1:
	raise Exception("Path to protocol file must be provided.")

if len(sys.argv) > 3:
	raise Exception("Too many arguments provided.")

# The first argument is the path to the protocol file
args = sys.argv[1:]
filename = args[0]

# If a second argument is provided, then this is the path to the executable
execpath = ".."
if len(args) > 1:
	execpath = args[1]
execpath = execpath + "/EGS"

joboptions = ''
paramstring = ''
nocluster = 0

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

		# Is there a special string saying we are not on the cluster?
		if re.match("nocluster", line):
			nocluster = 1

paramstring = re.sub(' $', '', paramstring)

# Make a template job file
job = job = "#!/bin/bash\n"
job = job + joboptions
job = job + execpath + ' ' + "parameters.txt\n"

# Make all parameter combinations
params = paramstring.split(' ')
dirnames = pyt.combine_parameters(".", params)

if nocluster:
	print("You are not on the cluster, I will not submit the jobs to SLURM.") 

# Copy the job file in all directories and submit to SLURM
for dir in dirnames:

	with open(dir + "/job.sh", "w+") as f:
		f.write(job)

	if not nocluster: 
		subprocess.run("cd " + dir + "\nsbatch ./job.sh", shell = True)