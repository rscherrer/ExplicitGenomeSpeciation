#!/usr/bin/python3

# This script is to create an element-wise difference of two arrays of floating point numbers
# Contained in two different binary data files

# Arguments:
# -d: directory
# 1: name of the first data file
# 2: name of the second data file
# 3: name of the output result 

import numpy as np
import sys

# Default arguments
directory = "."

# Read in arguments
args = sys.argv[1:]

if len(args) < 3:
	raise Exception("Wrong number of arguments")

inputfilename1 = args[0]
inputfilename2 = args[1]
outputfilename = args[2]

if "-d" in args:
	directory = args[args.index("-d") + 1]

# Prepare filenames
inputfilename1 = directory + "/" + inputfilename1
inputfilename2 = directory + "/" + inputfilename2
outputfilename = directory + "/" + outputfilename

# Read input files
with open(inputfilename1, "rb") as inputfile:
    data1 = inputfile.read()
data1 = np.frombuffer(data1, np.float64)

with open(inputfilename2, "rb") as inputfile:
    data2 = inputfile.read()
data2 = np.frombuffer(data2, np.float64)

# Produce output
output = [data1[i] - data2[i] for i in range(0, len(data1))]

# Save output file
with open(outputfilename, "wb") as outputfile:
	for out in output:
		outputfile.write(out)

