#!/usr/bin/python3

# Plot distribution across population through time
# Use this script from the command line with the data file to plot as argument

import numpy as np
import matplotlib.pyplot as plt
import sys

# Read time
with open("time.dat", "rb") as binary_file:
    time = binary_file.read()
time = np.frombuffer(time, np.float64)
time = [int(i) for i in time]

# Read counts
count = np.zeros((len(time), 4))
filenames = ["count00.dat", "count01.dat", "count10.dat", "count11.dat"]
for f in range(len(filenames)):
	with open(filenames[f], "rb") as binary_file:
		n = binary_file.read()
	count[:, f] = np.frombuffer(n, np.float64)

# Total counts
total = count.sum(axis=1)
total = [int(i) for i in total]

# Repeat time point per indidivual
time = [item for item, rep in zip(time, total) for i in range(rep)]

# Read data
filename = "population_x.dat"
if len(sys.argv) == 2: 
	filename = sys.argv[1]
with open(filename, "rb") as binary_file:
	x = binary_file.read()
x = np.frombuffer(x, np.float64)

# Plot
#fig = plt.figure(figsize=(20,0))
#ax1 = fig.add_subplot(111)
plt.hist2d(time, x, bins=30, cmap='Blues')
# plt.hist(time, bins=30)
# plt.savefig('plot.png', dpi=300, bbox_inches='tight')
plt.show()

