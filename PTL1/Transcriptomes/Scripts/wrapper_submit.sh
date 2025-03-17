#!/bin/bash
## Last updated Jan 2025 by Auden Cote-L'Heureux

## This script is intended to be used to process genomic CDS with EukPhylo part 1 on an HPC that uses the Slurm workload manager.
## The first part of the script are Slurm-specific parameters that should be adjusted by users to fit their resource allocation
## needs and restrictions, followed by some example commands taken from the GitHub Wiki, more detail for which can be found
## here: https://github.com/Katzlab/EukPhylo/wiki/EukPhylo-Part-1:-GF-assignment



## Example commands

# A simple run that goes from script 1 to script 7 (the last script) using the Universal genetic code

parent='/EukPhylo/PTL1/Transcriptomes/'
out_dir='/Output_data'
in_dir='/Input_data'

## EXAMPLE RUN COMMANDS BELOW

# A simple run of part 2, starting from ReadyToGo files and running through tree building
python3 ${parent}/Scripts/wrapper.py \
    --first_script 1 \
    --last_script 7 \
    --assembled_transcripts ${in_dir} \
    --databases ${parent}Databases \
    --output ${out_dir} > Output.out
