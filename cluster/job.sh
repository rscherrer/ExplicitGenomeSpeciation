#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=32Gb
#SBATCH --partition=gelifes
./EGS parameters.txt
