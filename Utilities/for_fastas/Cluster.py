'''
#Author, date: Godwin Ani and Laura Katz, 9th- Feb - 2023.
#Dependencies: Python3, CD-Hit
#Intent: For clustering nucleotide or amino acid sequences with the CD-Hit program.
#Inputs: A folder of containing  Amino acid or DNA fasta files.
#Outputs: A folder of clustered files.
#Example: python Cluster.py --type dna --identity 0.95 --overlap 0.67 --input input_folder_dna --output output_folder_dna
'''



import os
import argparse
from tqdm import tqdm
import subprocess

def input_validation(value, error_message):
    try:
        integer, fractional = value.split('.')
        value = float(value)
        if int(integer) == 0 and len(fractional) == 2:
            return value
    except ValueError:
        pass
    print(error_message)
    exit(1)

def cluster_sequences(program, identity, overlap, input_folder, output_folder):
    for file in tqdm(os.listdir(input_folder)):
        if file.endswith('.fasta'):
            subprocess.run([f'{program}', '-i', f'{input_folder}/{file}', '-o', f'{output_folder}/{file}', '-c', f'{identity}', '-d', '0', '-aS', f'{overlap}'])

    for file in os.listdir(output_folder):
        if file.endswith('.clstr'):
            os.rename(f'{output_folder}/{file}', f'{output_folder}/{file.split("FILE")[0]}Clustered.txt')

def main():
    parser = argparse.ArgumentParser(description='Cluster amino acid or DNA sequences using CD-HIT.')
    parser.add_argument('--type', choices=['aa', 'dna'], required=True, help='Type of sequences (aa for Amino Acids, dna for DNA)')
    parser.add_argument('--identity', type=str, required=True, help='Sequence Identity Threshold (e.g., 0.99, 0.95)')
    parser.add_argument('--overlap', type=str, required=True, help='Sequence Alignment Overlap Value (e.g., 0.67, 0.75)')
    parser.add_argument('--input_files', type=str, required=True, help='Input folder containing sequences in fasta format')
    parser.add_argument('--output', type=str, required=True, help='Output folder for clustered sequences')

    args = parser.parse_args()

    if not os.path.isdir(args.input_files):
        print(f'Error: Input folder "{args.input_files}" does not exist.')
        exit(1)

    if not os.path.isdir(args.output):
        os.mkdir(args.output)

    if args.type == 'aa':
        identity = input_validation(args.identity, 'ERROR! Use format 0.## for Amino acids sequence identity threshold.')
        overlap = input_validation(args.overlap, 'ERROR! Use format 0.## for Amino acids sequence alignment overlap value.')
        cluster_sequences('cd-hit', identity, overlap, args.input_files, args.output)
    elif args.type == 'dna':
        identity = input_validation(args.identity, 'ERROR! Use format 0.## for DNA sequence identity threshold.')
        overlap = input_validation(args.overlap, 'ERROR! Use format 0.## for DNA sequence alignment overlap value.')
        cluster_sequences('cd-hit-est', identity, overlap, args.input_files, args.output)
    else:
        print('Invalid sequence type. Choose "aa" for Amino Acids or "dna" for DNA.')
        exit(1)

if __name__ == "__main__":
    main()
