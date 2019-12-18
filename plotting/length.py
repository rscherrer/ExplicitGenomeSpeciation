#!/usr/bin/python3

# Print length of the data
# Use this script from the command line with the
# name of the binary file you want to read

import numpy as np
import sys

with open(sys.argv[1], "rb") as binary_file:
    # Read the whole file at once
    data = binary_file.read()

print(len(np.frombuffer(data, np.float64)))