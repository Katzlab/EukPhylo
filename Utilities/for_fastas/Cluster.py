'''
#Author, date: Godwin Ani and Laura Katz, Feb 9th 2023
#Modified: Adri Grow, April 6th 2025 to allow clustering at 100% (1.0) and output renamed file(s) with id clustered appended to file name
#Dependencies: Python3, CD-Hit
#Intent: For clustering nucleotide or amino acid sequences with the CD-Hit program
#Inputs: A folder of containing AA or DNA fasta files
#Outputs: A folder of clustered files
#Example: python Cluster.py -t dna -id 0.95 -ov 0.67 -i input_folder_dna -o output_folder_dna
'''

import os
import argparse
from tqdm import tqdm
import subprocess

def input_validation(value, error_message):
    try:
        value = float(value)
        if value == 1.0:
            return value
        integer, fractional = str(value).split('.')
        if int(integer) == 0 and len(fractional) == 2:
            return value
    except ValueError:
        pass
    print(error_message)
    exit(1)

def cluster_sequences(program, identity, overlap, input_folder, output_folder):
    for file in tqdm(os.listdir(input_folder)):
        if file.endswith('.fasta'):
            output_name = f"{os.path.splitext(file)[0]}_{int(float(identity) * 100)}clustered.fasta"
            subprocess.run([f'{program}', '-i', f'{input_folder}/{file}', '-o', f'{output_folder}/{output_name}', '-c', f'{identity}', '-d', '0', '-aS', f'{overlap}'])

    for file in os.listdir(output_folder):
        if file.endswith('.clstr'):
            base_name = os.path.splitext(file)[0]  # removes .clstr
            if base_name.endswith('.fasta'):
                base_name = base_name[:-6]  # removes .fasta from end
            new_name = f"{base_name}.txt"
            os.rename(f'{output_folder}/{file}', f'{output_folder}/{new_name}')

def main():
    parser = argparse.ArgumentParser(description='Cluster amino acid or nucleotide sequences using CD-HIT.')
    parser.add_argument('-t', '--type', choices=['aa', 'dna'], required=True, help='Type of sequences (aa for amino acid, dna for nucleotide)')
    parser.add_argument('-id','--identity',  type=str, required=True, help='Sequence identity threshold (e.g. 1.0, 0.99, 0.95)')
    parser.add_argument('-ov', '--overlap', type=str, required=True, help='Sequence alignment overlap value (e.g. 0.67, 0.75)')
    parser.add_argument('-i', '--input_files', type=str, required=True, help='Input folder containing sequences in fasta format')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output folder for clustered sequences ending with -id value')

    args = parser.parse_args()

    if not os.path.isdir(args.input_files):
        print(f'Error: Input folder "{args.input_files}" does not exist.')
        exit(1)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    if args.type == 'aa':
        identity = input_validation(args.identity, 'ERROR! Use format 0.## or 1.0 for amino acid sequence identity threshold.')
        overlap = input_validation(args.overlap, 'ERROR! Use format 0.## for amino acid sequence alignment overlap value.')
        cluster_sequences('cd-hit', identity, overlap, args.input_files, args.output)
    elif args.type == 'dna':
        identity = input_validation(args.identity, 'ERROR! Use format 0.## or 1.0 for nucleotide sequence identity threshold.')
        overlap = input_validation(args.overlap, 'ERROR! Use format 0.## for nucleotide sequence alignment overlap value.')
        cluster_sequences('cd-hit-est', identity, overlap, args.input_files, args.output)
    else:
        print('Invalid sequence type. Choose "aa" for amino acids or "dna" for nucleotides.')
        exit(1)

if __name__ == "__main__":
    main()
