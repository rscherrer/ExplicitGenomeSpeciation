#!/usr/bin/python3

# Plot distribution across population through time
# Use this script from the command line with the data file to plot as argument

# Arguments
# -d: directory
# -f: data file name

import numpy as np
import matplotlib.pyplot as plt
import sys

# Default arguments
directory = "."
datafilename = "population_x.dat"

# Read in arguments
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
countfilenames = ["count00.dat", "count01.dat", "count10.dat", "count11.dat"]
countfilenames = [directory + "/" + fname for fname in countfilenames]

# Read time
with open(timefilename, "rb") as timefile:
    time = timefile.read()
time = np.frombuffer(time, np.float64)
time = [int(i) for i in time]

# Read counts
counts = np.zeros((len(time), 4))
for f in range(len(countfilenames)):
	with open(countfilenames[f], "rb") as countfile:
		n = countfile.read()
	counts[:, f] = np.frombuffer(n, np.float64)

# Total counts
popsizes = counts.sum(axis=1)
popsizes = [int(i) for i in popsizes]

# Repeat time point per indidivual
time = [item for item, rep in zip(time, popsizes) for i in range(rep)]

# Read data
with open(datafilename, "rb") as datafile:
	data = datafile.read()
data = np.frombuffer(data, np.float64)

# Plot
#fig = plt.figure(figsize=(20,0))
#ax1 = fig.add_subplot(111)
plt.hist2d(time, data, bins=30, cmap='Blues')
# plt.hist(time, bins=30)
# plt.savefig('plot.png', dpi=300, bbox_inches='tight')
plt.show()

