#!/usr/bin/python3

# Plots a heatmap of a genome scan across loci through time
# Fst through time as a default
# Provide optional custom data file name as a command line argument

# Arguments
# -d: directory
# -f: data file name

import matplotlib.pyplot as plt
import numpy as np
import sys

# Default parameters
directory = "."
datafilename = "genome_Fst.dat"

# Read in arguments if any
if len(sys.argv) > 1:

	args = sys.argv[1:]

	# Check even number of arguments after the name of the script
	if len(args) % 2 != 0:
		raise Exception("Wrong number of arguments")

	if "-d" in args:
		directory = args[args.index("-d") + 1]

	if "-f" in args:
		datafilename = args[args.index("-f") + 1]

# Prepare file names
timefilename = directory + "/time.dat"
datafilename = directory + "/" + datafilename

# Read time
with open(timefilename, "rb") as timefile:
    time = timefile.read()
time = np.frombuffer(time, np.float64)
time = [int(i) for i in time]

# Read the data
with open(datafilename, "rb") as datafile:
	data = datafile.read()
data = np.frombuffer(data, np.float64)

ntimepoints = len(time)
if len(data) % ntimepoints != 0:
	raise Exception("Incorrect number of values in the data")
nloci = int(len(data) / ntimepoints)

# Make a table of zeros
heatmap = np.zeros([nloci, ntimepoints])

# Loop through timepoints
for t in range(ntimepoints):

	start = t * nloci
	end = start + nloci - 1

	scan = data[start:(end + 1)]

	# Replace columns with current genome scan
	heatmap[:,t] = scan

# Plot
plt.imshow(heatmap, cmap='hot', interpolation='nearest')
plt.show()