#!/bin/bash
## Last updated Jan 2025 by Auden Cote-L'Heureux; modified Feb 2025 by Adri K. Grow

## This shell script is used for running EukPhylo part 2, and includes a general setup for use on an HPC that uses 
## the Slurm workload manager. It also includes several example run commands, which correspond to examples explained in more detail in the 
## EukPhylo Wiki (https://github.com/Katzlab/EukPhylo/wiki/EukPhylo-Part-2:-MSAs,-trees,-and-contamination-loop).
## These run commands can also be copied and run in the terminal / command line separately, without a shell script.
## For the contamination loop, We recommend iterating the sister/subsisters loop multiple times as branches will shift. In contrast, we recommend only running clade grabbing once

## SLURM-SPECIFIC SETUP BELOW

############# UMass HPC (Unity) requirements below ##################### (DELETE section if not applicable):
#SBATCH --job-name=EukPhylo
#SBATCH -n 10 # Number of Cores per Task
#SBATCH --mem=125G # Requested Memory
#SBATCH -p cpu # Partition
#SBATCH -q long # long QOS
#SBATCH -t 334:00:00 # Job time limit
#SBATCH --output=Run_EP.%A_%a.out # Stdout (%j expands to jobId)
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@email.edu
#SBATCH --array=1-600%50
module purge       #Cleans up any loaded modules
module load conda/latest
module load mafft/7.505
module load diamond/2.1.7
conda activate /work/pi_lkatz_smith_edu/Conda_PTL6p2/envs/PTL/

############# Smith HPC (Grid) requirements below ##################### (DELETE section if not applicable):
#SBATCH --job-name=EukPhylo # Job name
#SBATCH --output=Run_EukPhylo.%j.out # Stdout (%j expands to jobId)
#SBATCH --nodes=1
#SBATCH --ntasks=10 ## On the Smith College HPC (Grid), we have to change this to be double the number of task/batches you want to launch
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email@email.edu ##add your email address for job updates
#Load required modules
module purge       #Cleans up any loaded modules
module use /gridapps/modules/all    #make sure module locations is loaded
module load slurm
module load ETE
module load Biopython/1.79-foss-2021b
module load DIAMOND/2.0.13-GCC-11.2.0
module load MAFFT
module load RAxML
module load IQ-TREE/2.1.2-gompi-2021b
module load tqdm/4.64.1-GCCcore-12.2.0
module load Python/3.9.6-GCCcore-11.2.0
module load Guidance_mid #Smith College HPC specific
export PATH=$PATH:/beegfs/fast/katzlab/grid_phylotol_setup/programs/standard-RAxML-master #Smith College HPC specific #export PATH=$PATH:/Path/To/Executable/Files


## Provide your parent path
parent='/Your/Home/Folder/' # The folder where you are running EukPhylo (this should contain the Scripts and input data folders)

## EXAMPLE RUN COMMANDS BELOW

# A simple run of part 2, starting from ReadyToGo files and running through tree building
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/eukphylo.py --start raw --end trees --gf_list ${parent}listofOGs.txt --taxon_list ${parent}taxon_list.txt --data ${parent}Input_folder --output ${parent}Output_folder > Output.out

# Another example starting from ReadyToGo files and running through tree building, with the commonly used similarity filter cutoff, blacklist, and  "sim_taxa_list" arguments (see Wiki)
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/eukphylo.py --start raw --end trees --gf_list ${parent}listofOGs.txt --taxon_list ${parent}taxon_list.txt --data ${parent}Input_folder --output ${parent}Output_folder --similarity_filter --blacklist ${parent}Blacklist.txt --sim_cutoff 0.99 --sim_taxa sim_taxa_list.txt > Output.out

# An example of running just the concatenation step of part 2, starting from trees
srun --exact -n 1 -D ${parent} python3 ${parent}Scripts/eukphylo.py --start trees --concatenate --concat_target_taxa Sr_rh --data ${parent}Output > log.out

# See the Wiki (https://github.com/Katzlab/EukPhylo/wiki/EukPhylo-Part-2:-MSAs,-trees,-and-contamination-loop) for more details!




