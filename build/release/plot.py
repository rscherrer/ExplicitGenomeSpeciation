#!/usr/bin/python3

# Plot simulation output
# Use this script from the command line with the name of
# the data files you want to plot as argument

# Example use: to plot ecological, reproductive and spatial isolation on the same plot, use
# ./plot.py EI.dat RI.dat SI.dat

import numpy as np
import matplotlib.pyplot as plt
import sys

# Read time points
with open("time.dat", "rb") as binary_file:
    time = binary_file.read()
time = np.frombuffer(time, np.float64)

# Read provided arguments
filenames = sys.argv
filenames.pop(0)

# If no argument provided
if len(filenames) == 0:
	filenames = ["EI.dat"]

# Loop through program arguments
for f in filenames:

	# Read data from file
	with open(f, "rb") as binary_file:
		y = binary_file.read()
	y = np.frombuffer(y, np.float64)

	# Add to plot
	plt.plot(time, y)

# Display
plt.show()

