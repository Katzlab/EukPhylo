#!/bin/bash
## Last updated Jan 2025 by Auden Cote-L'Heureux

## This script is intended to be used to process genomic CDS with EukPhylo part 1 on an HPC that uses the Slurm workload manager.
## The first part of the script are Slurm-specific parameters that should be adjusted by users to fit their resource allocation
## needs and restrictions, followed by some example commands taken from the GitHub Wiki, more detail for which can be found
## here: https://github.com/Katzlab/EukPhylo/wiki/EukPhylo-Part-1:-GF-assignment


## Example run command

# Start at script 1 and go through script 5 (the final script) using the Universal genetic code
python3 Scripts/wrapper.py -1 1 -2 5 --cds Input -o Output --genetic_code Universal --databases Databases > log.out

