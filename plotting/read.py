#!/usr/bin/python3

# Print values to the screen
# Use this script from the command line with the
# name of the binary file you want to read

# Arguments
# -d: directory
# -f: file name

import numpy as np
import sys

# Default arguments
directory = "."
filename = "time.dat"

# Read in arguments
if len(sys.argv) > 1:

	args = sys.argv[1:]

	if "-d" in args:
		directory = args[args.index("-d") + 1]

	if "-f" in args:
		filename = args[args.index("-f") + 1]

# Prepare file names
filename = directory + "/" + filename

with open(filename, "rb") as file:
    data = file.read()

print(np.frombuffer(data, np.float64))
