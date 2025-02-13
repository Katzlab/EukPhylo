'''
Author & Date: Adri K. Grow + ChatGPT, Nov 11th 2024
    - Updated 02/13/25 to accept either transcriptome and genome assembled data in command line
Motivation: assess and rename assembled transcript or genome files for PTL6p1 (or EPv1p1)
Intention: warn if any 'transcripts.fasta' or 'contigs.fasta' files are missing or empty for an LKH, otherwise rename and copy them with their assigned 10-digit code by LKH
Input:
    - a base directory containing subdirectories for each LKH assembled file, named 'WTA_LKH<xxxx>' or 'WGA_LKH<xxxx>', each containing a 'transcripts.fasta' or 'contigs.fasta' file
    - a mapping .txt file with LKH#s tab-separated with corresponding 10-digit codes
Output: 
    - a folder named 'renamed_transcripts|contigs' with assembled files now named by 10-digit codes; e.g. "Sr_rh_Ro04_assembledTranscripts.fasta"
Dependencies: python3
Usage: 
    - for transcriptomes: python3 ProcessAndRenameAssembledData.py <assembled transcriptomes directory> <mapping_file.txt> transcriptomes
    - for genomes: python3 ProcessAndRenameAssembledData.py <assembled genomes directory> <mapping_file.txt> genomes
'''

import os
import shutil
import sys

def read_lkh_mapping(mapping_file):
    """Reads the LKH number to 10-digit code mapping from a file."""
    mapping = {}
    with open(mapping_file, 'r') as file:
        for line in file:
            lkh_number, code = line.strip().split('\t')
            mapping[lkh_number] = code
    return mapping

def process_directory(base_dir, mapping, output_dir, data_type):
    """Iterates over all subdirectories in base_dir, processes relevant fasta files."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create output directory if it doesn't exist
    
    # Set file naming patterns based on data type
    folder_prefix = "WTA_LKH" if data_type == "transcriptomes" else "WGA_LKH"
    fasta_filename = "transcripts.fasta" if data_type == "transcriptomes" else "contigs.fasta"
    output_suffix = "_assembledTranscripts.fasta" if data_type == "transcriptomes" else "_assembledContigs.fasta"
    
    for folder_name in os.listdir(base_dir):
        folder_path = os.path.join(base_dir, folder_name)

        if os.path.isdir(folder_path) and folder_name.startswith(folder_prefix):
            lkh_number = folder_name.split('_')[1]  # Extract LKH number from folder name
            fasta_file = os.path.join(folder_path, fasta_filename)
            
            if not os.path.isfile(fasta_file):
                print(f"    WARNING: file '{fasta_filename}' is missing in folder {folder_name}.")
                continue
            
            if os.path.getsize(fasta_file) == 0:
                print(f"    WARNING: file '{fasta_filename}' is empty in folder {folder_name}.")
                continue
            
            if lkh_number in mapping:
                new_name = f"{mapping[lkh_number]}{output_suffix}"
                output_path = os.path.join(output_dir, new_name)
                shutil.copy(fasta_file, output_path)
            else:
                print(f"Notification: No 10-digit code found for LKH number {lkh_number} in folder {folder_name}.")

def main():
    if len(sys.argv) != 4 or sys.argv[3] not in ["transcriptomes", "genomes"]:
        print("Usage: python script.py <base_dir> <mapping_file> <transcriptomes|genomes>")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    mapping_file = sys.argv[2]
    data_type = sys.argv[3]
    
    if not os.path.isdir(base_dir):
        print(f"Error: The directory '{base_dir}' does not exist.")
        sys.exit(1)

    if not os.path.isfile(mapping_file):
        print(f"Error: The file '{mapping_file}' does not exist.")
        sys.exit(1)
    
    output_dir = os.path.join(os.getcwd(), "renamed_transcripts" if data_type == "transcriptomes" else "renamed_contigs")
    
    mapping = read_lkh_mapping(mapping_file)
    process_directory(base_dir, mapping, output_dir, data_type)

if __name__ == "__main__":
    main()
