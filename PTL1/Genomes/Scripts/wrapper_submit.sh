#!/bin/bash
## Last updated Jan 2025 by Auden Cote-L'Heureux

## This script is intended to be used to process genomic CDS with EukPhylo part 1 on an HPC that uses the Slurm workload manager.
## The first part of the script are Slurm-specific parameters that should be adjusted by users to fit their resource allocation
## needs and restrictions, followed by some example commands taken from the GitHub Wiki, more detail for which can be found
## here: https://github.com/Katzlab/EukPhylo/wiki/EukPhylo-Part-1:-GF-assignment

## Slurm specific code

#SBATCH --job-name=EukPhylo
#SBATCH --output=EukPhylo.%j.out # Stdout (%j expands to jobId)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=64 # #change to number of srun when running multiple instances
#SBATCH --mem=160G

module purge       #Cleans up any loaded modules
module load slurm
module load tqdm/4.62.3-GCCcore-11.2.0
module load Biopython/1.79-foss-2021b
module load BLAST+/2.12.0-gompi-2021b
module load DIAMOND/2.0.13-GCC-11.2.0
module load VSEARCH/2.22.1-GCC-11.3.0

parent='/Your/Home/Folder'

## Example run command

# Start at script 1 and go through script 5 (the final script) using the Universal genetic code
srun -D ${parent}Scripts python3 ${parent}Scripts/wrapper.py -1 1 -2 5 --cds ${parent}Input -o ${parent}Output --genetic_code Universal --databases ${parent}Databases > log.out

