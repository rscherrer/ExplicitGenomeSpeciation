#!/usr/bin/python3

# Plot a histogram of trait values between ecotypes in the last generation
# If called without arguments, will plot the distribution in ecological trait values at the last generation

# Arguments
# -d: directory
# -f: data file name
# -t: timepoint

import numpy as np
from matplotlib import pyplot
from itertools import compress
import sys

# Default arguments
directory = "."
timepoint = -1
datafilename = "population_x.dat"

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
	
	if "-t" in args:
		timepoint = int(args[args.index("-t") + 1])

# Prepare file names
timefilename = directory + "/time.dat"
datafilename = directory + "/" + datafilename
countfilenames = ["count00.dat", "count01.dat", "count10.dat", "count11.dat"]
countfilenames = [directory + "/" + fname for fname in countfilenames]
ecofilename = directory + "/population_ecotype.dat"

# Read time
with open(timefilename, "rb") as timefile:
    time = timefile.read()
time = np.frombuffer(time, np.float64)
time = [int(i) for i in time]

if timepoint == -1:
	timepoint = time[-1]

if timepoint not in time:
	raise Exception("Timepoint not found, try another one.")

# Extract timepoint index
t = time.index(timepoint)

# Read counts
counts = np.zeros((len(time), 4))
for f in range(len(countfilenames)):
	with open(countfilenames[f], "rb") as countfile:
		n = countfile.read()
	counts[:, f] = np.frombuffer(n, np.float64)

# Total population size per time step
popsizes = counts.sum(axis=1)
popsizes = [int(i) for i in popsizes]

# Population size in the timestep of interest
n = popsizes[t]

# First and last individual of the timepoint we consider
start = sum(popsizes[0:t])
end = start + n - 1

# Read the data to plot
with open(datafilename, "rb") as datafile:
    data = datafile.read()
data = np.frombuffer(data, np.float64)
data = data[start:(end + 1)]

# Read the ecotype labels of the individuals
with open(ecofilename, "rb") as ecofile:
    eco = ecofile.read()
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

