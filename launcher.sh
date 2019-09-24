#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=00:10:00
#SBATCH --mem=10Gb
#SBATCH --partition=short

./EGS parameters.txt
