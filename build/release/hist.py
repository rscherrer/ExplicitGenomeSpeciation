#!/usr/bin/python3

# Use to check random stuff

# Plot a histogram of trait values between ecotypes in the last generation
# If called without arguments, will plot the distribution in ecological trait values at the last generation
# Provide an optional data file to plot preceded by -f
# Provide an optional timestep preceded by -t

import numpy as np
from matplotlib import pyplot
from itertools import compress
import sys

# Default arguments
t = -1
datafile = "population_x.dat"

# Read in arguments if any
if len(sys.argv) > 1:

	args = sys.argv[1:]

	# Check even number of arguments after the name of the script
	if len(args) % 2 != 0:
		raise Exception("Wrong number of arguments")

	nargs = int(len(args) / 2)

	for i in range(nargs):

		# If a file name is supplied
		if args[i] == "-f":
			datafile = args[i + 1]

		# Else if a time point is provided
		elif args[i] == "-t":
			t = int(args[i + 1])

# Read time
with open("time.dat", "rb") as binary_file:
    time = binary_file.read()
time = np.frombuffer(time, np.float64)
time = [int(i) for i in time]

# Read counts
countfiles = ["count00.dat", "count01.dat", "count10.dat", "count11.dat"]
counts = np.zeros((len(time), 4))
for f in range(len(countfiles)):
	with open(countfiles[f], "rb") as binary_file:
		n = binary_file.read()
	counts[:, f] = np.frombuffer(n, np.float64)

# Total population size per time step
popsizes = counts.sum(axis=1)
popsizes = [int(i) for i in popsizes]

# Population size in the timestep of interest
if t > 0: 
	t = t - 1 # makes more sense to the user to not count in python style
n = popsizes[t]

# First and last individual of the timepoint we consider
start = sum(popsizes[0:t])
end = start + n - 1

# Read the data to plot
with open(datafile, "rb") as binary_file:
    data = binary_file.read()
data = np.frombuffer(data, np.float64)
data = data[start:(end + 1)]

# Read the ecotype labels of the individuals
with open("population_ecotype.dat", "rb") as binary_file:
    eco = binary_file.read()
eco = np.frombuffer(eco, np.float64)
eco = list(map(bool, eco))
eco = eco[start:(end + 1)]

# Extract vectors of values for both ecotypes
data1 = list(compress(data, eco))
data0 = list(compress(data, [not i for i in eco]))

# Plot the histogram
bins = np.linspace(-1.5, 1.5, 100)
pyplot.hist(data1, bins, alpha=0.5, label='Ecotype 1')
pyplot.hist(data0, bins, alpha=0.5, label='Ecotype 0')
pyplot.legend(loc='upper right')
pyplot.show()

