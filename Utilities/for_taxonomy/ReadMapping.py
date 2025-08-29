'''
#Author, date: ?
#Uploaded: updated by Adri Grow, Aug 2025
#Intent: map a group of trimmed reads to a reference 
#Dependencies: Python, HISAT2, samtools, (optional: sambamba)
#EDIT LINES: 19 & 36
#Inputs: Folder named 'TrimmedReads' containing the forward and reverse trimmed reads that start with the same unique identifier for each sample/cell
#Outputs: Folders with the names of the unique identifier (e.g. LKHs) containing the bam files
#Usage: python3 ReadMapping.py
#IMPORTANT: Lines 34-42 manipulate the output files in several different ways including converting .sam to .bam, sorting, optional deduplicating, optional quality filtering, and retaining only mapped reads. It is the responsibility of the user to determine exactly which commands are needed for their dataset.
'''

import os
from Bio import SeqIO

#This first command builds your reference with HISAT
#If you've already done this, DON'T run this command! Instead, comment it out (use a # in front of it)
#It will output several files. Don't worry about them, HISAT will know what to do
os.system("hisat2-build Foram_reference.fasta Foram_Index") #Replace "Foram_reference.fasta" with your reference fasta name, and optionally change "Foram_Index" to your preferred index name

folder = os.listdir("TrimmedReads") #Replace "TrimmedReads" with the name of the folder containing your trimmed reads, if different than TrimmedReads

folder.sort() #This sorts the trimmed reads folder so that all the files are passed in order

for x in folder:
	#This is specific for file names starting with 'LKH' unqiue identifiers formatted similar to 'LKH###_FPE.fastq.gz'
	if "LKH" in x and "FPE" in x: #Assigning a variable to forward reads. Make sure you have both forward and reverse reads for each cell!
		FPE = x
		sample_id = FPE.split("_FPE")[0]
	if "LKH" in x and "RPE" in x: #Assigning a variable to reverse reads
		RPE = x

		if FPE.split("_FPE")[0] == RPE.split("_RPE")[0]: #Match sample IDs dynamically 
			#The next few lines are several HISAT commands that will create new files
			#If necessary, EDIT the name of the index and the name of the trimmed reads folder in the very next line only
			os.system("hisat2 -x Foram_Index -1 TrimmedReads/" +FPE+ " -2 TrimmedReads/" +RPE+ " -S sample.sam") #running HISAT2
			os.system("samtools view -bS sample.sam > sample.bam") #converts .sam file to .bam file
			os.remove("sample.sam") #remove the .sam file (already converted to .bam, sam files are large and unnecessary to keep)
			#os.system("samtools fixmate -O bam sample.bam  fixmate_sample.bam") #use this command if you will be using the sambamba markdup command to remove duplicate reads (Katzlab default for transcriptomics and amplicon is to not remove duplicates)
			os.system("samtools sort -O bam -o sorted_sample.bam sample.bam") #sorts the .bam file alignments by leftmost coordinates
			#os.system("sambamba markdup -r sorted_sample.bam sorted_sample.dedup.bam") #removes duplicate reads - may not be appropriate for your study or protocols, user will need to determine if this is best practice for their study
			#os.system("samtools view -h -b -q 40 sorted_sample.dedup.bam > sorted_sample.q40.bam") #only keeps reads with mapping quality ≥ 40, input is the dedup file but can easily be modified to use the sorted .bam file
			#os.system("samtools view -h -b -q 20 sorted_sample.dedup.bam > sorted_sample.q20.bam") #only keeps reads with mapping quality ≥ 20, input is the dedup file but can easily be modified to use the sorted .bam file
			os.system("samtools view -h -F 4 -b sorted_sample.bam > sorted_mapped_sample.bam") #only keeps mapped reads, using the sorted .bam file as input - this is the Katzlab transcriptomic and amplicon final output that should be used for continued analyses

			if not os.path.isdir(sample_id):
				os.mkdir(sample_id) #making folders with the names of the LKHs or unique identifiers
				
			for file in os.listdir('.'): #These lines move the bam files created into the new LKH/unique identifier folders
				if(file.endswith('.sam') or file.endswith('.bam')):
					os.rename(file, f"{sample_id}/{file}")
				
print("~~~~~~~~~~~:>~") #When the snake appears in terminal, the script has finished running for all samples/cells!
		
	
	
