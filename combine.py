#!/usr/bin/python3

# This script produces multiple combinations of parameters
# It makes one folder per parameter combination
# It takes as arguments a list of parameter names to combine
# And for each parameter a list of values to take
# It combines values of all parameters

import re
import sys
import os
import pytilities as pyt

# Read input arguments
args = sys.argv[1:]

# Make parameter combinations
pyt.combine_parameters(".", args)

