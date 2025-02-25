#!/bin/bash
## Last updated Jan 2025 by Auden Cote-L'Heureux; modified Feb 2025 by Adri K. Grow

## This shell script is used for running EukPhylo part 2, and includes a general setup for use on an HPC that uses 
## the Slurm workload manager. It also includes several example run commands, which correspond to examples explained in more detail in the 
## EukPhylo Wiki (https://github.com/Katzlab/EukPhylo/wiki/EukPhylo-Part-2:-MSAs,-trees,-and-contamination-loop).
## These run commands can also be copied and run in the terminal / command line separately, without a shell script.

parent='/PhyloToL-6/PTL2'
out_dir='/Output_data'
in_dir='/Input_data'

## EXAMPLE RUN COMMANDS BELOW

# A simple run of part 2, starting from ReadyToGo files and running through tree building
python3 ${parent}Scripts/eukphylo.py \
    --start raw \
    --end trees \
    --gf_list ${parent}listofOGs.txt \
    --taxon_list ${parent}taxon_list.txt \
    --data ${in_dir} \
    --output ${out_dir} > Output.out
