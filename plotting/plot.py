#!/usr/bin/python3

# Plot simulation output
# Use this script from the command line with the name of
# the data files you want to plot as argument

# Example use: to plot ecological, reproductive and spatial isolation on the same plot, use
# ./plot.py EI.dat RI.dat SI.dat

# Arguments
# -d: directory (should come first)
# -f: data file names

import numpy as np
import matplotlib.pyplot as plt
import sys

# Default arguments
directory = "."
datafilenames = ["EI.dat"]

# Read in arguments
if len(sys.argv) > 1:

	args = sys.argv[1:]

	if "-d" in args:
		directory = args[args.index("-d") + 1]

	if "-f" in args:
		fidx = args.index("-f")
		datafilenames = args[(fidx + 1):len(args)]

# Prepare file names
timefilename = directory + "/time.dat"
datafilenames = [directory + "/" + fname for fname in datafilenames]

# Read time points
with open(timefilename, "rb") as timefile:
    time = timefile.read()
time = np.frombuffer(time, np.float64)

# Loop through program arguments
for fname in datafilenames:

	# Read data from file
	with open(fname, "rb") as datafile:
		data = datafile.read()
	data = np.frombuffer(data, np.float64)

	# Add to plot
	plt.plot(time, data)

# Display
plt.show()

