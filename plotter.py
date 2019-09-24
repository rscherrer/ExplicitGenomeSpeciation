#!/usr/bin/python3

# Plot simulation output
# Use this script from the command line with the name of
# the data file you want to plot as argument

import matplotlib.pyplot as plt
import numpy as np
import sys


with open("time.dat", "rb") as binary_file:
    # Read the whole file at once
    time = binary_file.read()

time = np.frombuffer(time, np.float64)

with open(sys.argv[1], "rb") as binary_file:
	data = binary_file.read()

data = np.frombuffer(data, np.float64)

fig = plt.figure(figsize=(20,0))
ax1 = fig.add_subplot(111)

plt.plot(time, data)
#plt.savefig('plot.png', dpi=300, bbox_inches='tight')
plt.show()
