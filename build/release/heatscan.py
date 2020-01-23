#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import sys

# Plots a heatmap of a genome scan across loci through time
# Fst through time as a default
# Provide optional custom data file name as a command line argument

filename = "genome_Fst.dat"

# Read in arguments if any
if len(sys.argv) > 1:
	filename = sys.argv[1]

# Read time
with open("time.dat", "rb") as binary_file:
    time = binary_file.read()
time = np.frombuffer(time, np.float64)
time = [int(i) for i in time]

# Read the data
with open(filename, "rb") as datafile:
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