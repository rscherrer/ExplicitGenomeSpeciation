#!/usr/bin/python3

# This script is to reset the parameter file to default parameter values

# Read parameter file
with open("parameters.txt", "rt") as f:
	data = f.read()	

data = "rdynamics 1\n"
data = data + "trenewal 0.001\n"
data = data + "capacity 100\n"
data = data + "replenish 1\n"
data = data + "hsymmetry 1\n"
data = data + "ecosel 1.4\n"
data = data + "dispersal 0.001\n"
data = data + "birth 4\n"
data = data + "survival 0.8\n"
data = data + "sexsel 10\n"
data = data + "matingcost 0.01\n"
data = data + "maxfeed 0.0004\n"
data = data + "demesizes 100 0\n"
data = data + "nvertices 50 50 50\n"
data = data + "nedges 0 0 0\n"
data = data + "nchrom 3\n"
data = data + "mutation 0.00001\n"
data = data + "recombination 3\n"
data = data + "allfreq 0.2\n"
data = data + "scaleA 1 1 1\n"
data = data + "scaleD 0 0 0\n"
data = data + "scaleI 0 0 0\n"
data = data + "scaleE 0 0 0\n"
data = data + "skews 1 1 1\n"
data = data + "effectshape 2\n"
data = data + "effectscale 1\n"
data = data + "interactionshape 5\n"
data = data + "interactionscale 1\n"
data = data + "dominancevar 1\n"
data = data + "tburnin 1000\n"
data = data + "tend 3000\n"
data = data + "tsave 10\n"
data = data + "record 1\n"
data = data + "archsave 0\n"
data = data + "ntrials 100\n"

# Overwrite parameter file
with open("parameters.txt", "wt") as f:
	f.write(data)	