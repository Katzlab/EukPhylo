#!/bin/bash
#
#SBATCH --job-name=Gigi_spades
#SBATCH --output=rnaSPAdes_run.%j.out # Stdout (%j expands to jobId)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem=180G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xxx@xxx.edu

module purge       #Cleans up any loaded modules
module load SPAdes

rnaspades.py -m 500 -t 50 -1 ToAssemble/SRR26595464_FPE_q24_minlen75.fastq.gz -2 ToAssemble/SRR26595464_RPE_q24_minlen75.fastq.gz -o Assembled/SRR26595464
rnaspades.py -m 500 -t 50 -1 ToAssemble/SRR26595465_FPE_q24_minlen75.fastq.gz -2 ToAssemble/SRR26595465_RPE_q24_minlen75.fastq.gz -o Assembled/SRR26595465
rnaspades.py -m 500 -t 50 -1 ToAssemble/SRR26595468_FPE_q24_minlen75.fastq.gz -2 ToAssemble/SRR26595468_RPE_q24_minlen75.fastq.gz -o Assembled/SRR26595468

