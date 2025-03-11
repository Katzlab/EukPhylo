# Last updated Apr 2 2024
# Authors: Auden Cote-L'Heureux and Mario Ceron-Romero

# This script runs Guidance in an iterative fashion for more both MSA construction 
# and more rigorous homology assessment than what is offered in EukPhylo part 1.
# Guidance runs until the input number of iterations (--guidance_iters, default = 5) 
# has been reached, or until there are no sequences below the sequence score cutoff.
# All sequences below the score cutoff (--seq_cutoff, default = 0.3) are removed at
# each iteration. By default, EukPhylo does not remove residues that fall below the 
# given residue cutoff (--res_cutoff) and columns that fall below the given column 
# cutoff (--col_cutoff, defaults are 0), though this can be turned on by adjusting 
# these parameters. Outputs at this point are found in the “Guidance_NotGapTrimmed” 
# output folder. We then run MSAs through TrimAl to remove all sites in the alignment 
# that are at least 95% gaps (or --gap_trim_cutoff) generating files in the “Guidance” 
# output folder.

# This step is either intended to be run starting with --start = unaligned (but not raw)
# inputs, meaning one amino acid alignment per OG. It can also be run directly after the
# preguidance step. The run() function is called in two places: in eukphylo.py generally,
# and in contamination.py if the contamination loop is using Guidance as the re-alignment
# method.

#Dependencies
import os, sys, re
import subprocess
import logging
from Bio import SeqIO

# Set up logging configuration
logging.basicConfig(filename='guidance_run.log', level=logging.DEBUG,
                    format='%(asctime)s - %(levelname)s - %(message)s')

# Function to run system commands and capture logs
def run_command(command, log_file):
    try:
        logging.info(f"Running command: {command}")
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        logging.info(f"Command output: {result.stdout}")
        return result
    except subprocess.CalledProcessError as e:
        logging.error(f"Command failed with error: {e.stderr}")
        with open(log_file, 'a') as log:
            log.write(f"Command failed: {command}\nError: {e.stderr}\n")
        return None

# Called in eukphylo.py and contamination.py
def run(params):
    if params.start == 'raw' or params.start == 'unaligned':
        # Checking that pre-Guidance has been run or that unaligned files per OG are provided.
        if params.start == 'raw':
            preguidance_path = params.output + '/Output/Pre-Guidance'
        else:
            preguidance_path = params.data

        if not os.path.isdir(preguidance_path):
            logging.error(f"ERROR: The path {preguidance_path} could not be found when trying to locate pre-Guidance (unaligned) files.")
            exit()

        if len([f for f in os.listdir(preguidance_path) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]) == 0:
            logging.error(f"ERROR: No pre-Guidance (unaligned) files found at {preguidance_path}.")
            exit()

        # Creating intermediate folders
        os.mkdir(params.output + '/Output/Intermediate/Guidance')
        os.mkdir(params.output + '/Output/Intermediate/Guidance/Input')
        os.mkdir(params.output + '/Output/Intermediate/Guidance/Output')

        guidance_input = params.output + '/Output/Intermediate/Guidance/Input/'
        subprocess.run(f'cp -r {preguidance_path}/* {guidance_input}', shell=True)

        guidance_removed_file = open(params.output + '/Output/GuidanceRemovedSeqs.txt', 'w')
        guidance_removed_file.write('Sequence\tScore\n')

        too_many_seqs = False

        # Check if any files have more than 2000 sequences
        for file in [f for f in os.listdir(guidance_input) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]:
            nseqs = len([rec for rec in SeqIO.parse(guidance_input + '/' + file, 'fasta')])
            if nseqs > 2000:
                too_many_seqs = True
                break

        if too_many_seqs and not params.allow_large_files:
            return False

        # Process each unaligned AA fasta file
        for file in [f for f in os.listdir(guidance_input) if f.endswith('.fa') or f.endswith('.faa') or f.endswith('.fasta')]:
            tax_guidance_outdir = params.output + '/Output/Intermediate/Guidance/Output/' + file.split('.')[0].split('_preguidance')[0]
            os.mkdir(tax_guidance_outdir)

            fail = False
            # Iterative process
            for i in range(params.guidance_iters):
                n_recs = len([r for r in SeqIO.parse(guidance_input + '/' + file, 'fasta')])

                if n_recs < 4:
                    logging.warning(f"Gene family {file.split('.')[0].split('_preguidance')[0]} contains fewer than 4 sequences after {i} iterations.")
                    subprocess.run(f'rm -rf {tax_guidance_outdir}', shell=True)
                    if i == 0:
                        fail = True
                    break

                # Determining MAFFT algorithm
                mafft_alg = 'genafpair' if n_recs < 200 else 'auto'

                # Running Guidance
                command = f"perl /guidance/guidance.v2.02/www/Guidance/guidance.pl --seqFile {guidance_input}/{file} --msaProgram MAFFT --seqType aa --outDir {tax_guidance_outdir} --seqCutoff {params.seq_cutoff} --colCutoff {params.col_cutoff} --outOrder as_input --bootstraps 10 --MSA_Param '\\--{mafft_alg} --maxiterate 1000 --thread {params.guidance_threads} --bl 62 --anysymbol' > {params.output}/Output/Intermediate/Guidance/Output/{file[:10]}/log.txt"
                result = run_command(command, f'{params.output}/Output/Intermediate/Guidance/Output/{file[:10]}/log.txt')

                if result is None:
                    fail = True
                    break

                if os.path.isfile(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names'):
                    seqs_below = len([line for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:-1] if float(line.split()[-1]) < params.seq_cutoff])

                    if n_recs - seqs_below < 4:
                        logging.warning(f"Gene family {file.split('.')[0].split('_preguidance')[0]} contains fewer than 4 sequences after {i + 1} iterations.")
                        subprocess.run(f'rm -rf {tax_guidance_outdir}', shell=True)
                        break

                    if seqs_below == 0 or i == params.guidance_iters - 1:
                        logging.info(f"Guidance complete after {i + 1} iterations for gene family {file.split('.')[0].split('_preguidance')[0]}")
                        break

                    for line in [line for line in open(tax_guidance_outdir + '/MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names').readlines()[1:-1] if float(line.split()[-1]) < params.seq_cutoff]:
                        guidance_removed_file.write(line)

                    subprocess.run(f'cp {tax_guidance_outdir}/Seqs.Orig.fas.FIXED.Without_low_SP_Seq.With_Names {guidance_input}/{file}', shell=True)

                    # Handling intermediate files
                    if params.keep_iter:
                        iteration_folder = params.output + '/Output/Intermediate/Guidance/Iterations/' + str(i + 1) + '/' + file.split('.')[0].split('_preguidance')[0]
                        os.makedirs(iteration_folder, exist_ok=True)
                        subprocess.run(f'cp -r {tax_guidance_outdir}/* {iteration_folder}', shell=True)

                        if not params.keep_temp:
                            for gdir_file in os.listdir(iteration_folder):
                                if gdir_file not in ('MSA.MAFFT.Guidance2_res_pair_seq.scr_with_Names', 'MSA.MAFFT.aln.With_Names', 'MSA.MAFFT.Guidance2_res_pair_col.scr', 'log', 'postGuidance_preTrimAl_unaligned.fasta'):
                                    subprocess.run(f'rm -r {iteration_folder}/{gdir_file}', shell=True)
                                else:
                                    if gdir_file == 'MSA.MAFFT.aln.With_Names':
                                        subprocess.run(f'mv {iteration_folder}/{gdir_file} {iteration_folder}/{file.split(".")[0].split("_preguidance")[0]}_{gdir_file}.aln', shell=True)
                                    else:
                                        subprocess.run(f'mv {iteration_folder}/{gdir_file} {iteration_folder}/{file.split(".")[0].split("_preguidance")[0]}_{gdir_file}', shell=True)

                    subprocess.run(f'rm -r {tax_guidance_outdir}/*', shell=True)
                else:
                    fail = True
                    break

            if not fail:
                logging.info(f"Residue and column cutoffs applied for {file.split('.')[0].split('_preguidance')[0]}")
                # Additional steps to apply residue/column cutoffs and finalize the alignment can go here

        guidance_removed_file.close()
        return True
