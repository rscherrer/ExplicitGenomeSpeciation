#!/usr/bin/python3

# Plot a scan of a statistic across the genome at a given generation
# By default plots Fst at the last generation
# Provide an optional data file to plot preceded by -f
# Provide an optional timestep preceded by -t

import numpy as np
from matplotlib import pyplot
from itertools import compress
import os
import re
import sys

# Default arguments
timepoint = -1
filename = "genome_Fst.dat"

# Read in arguments if any
if len(sys.argv) > 1:

	args = sys.argv[1:]

	# Check even number of arguments after the name of the script
	if len(args) % 2 != 0:
		raise Exception("Wrong number of arguments")

	if "-f" in args:
		filename = args[args.index("-f") + 1]
	
	if "-t" in args:
		timepoint = int(args[args.index("-t") + 1])

# Read time
with open("time.dat", "rb") as timefile:
    time = timefile.read()
time = np.frombuffer(time, np.float64)
time = [int(i) for i in time]

if timepoint == -1:
	timepoint = time[-1]

if timepoint not in time:
	raise Exception("Timepoint not found, try another one.")

# Extract timepoint index
t = time.index(timepoint)

# Read the data
with open(filename, "rb") as datafile:
	data = datafile.read()
data = np.frombuffer(data, np.float64)

# I am so stupid
# The number of genes is the number of values divided by the number of time points!

if len(data) % len(time) != 0:
	raise Exception("Wrong number of values in the data file.")

nloci = int(len(data) / len(time))	

# Extract the portion of the data we want to plot
start = t * nloci
end = start + nloci - 1
data = data[start:(end + 1)]

# Plot
pyplot.stem(data)
pyplot.show()
