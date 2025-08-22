'''
#Author, date: ?
#Uploaded: updated by Adri Grow, Aug 2025
#Intent: map a group of trimmed reads to a reference 
#Dependencies: Python3, hisat2, samtools, sambamba
#EDIT LINES: 18 & 33
#Inputs: Folder named 'TrimmedReads' containing all the trimmed reads
#Outputs: Folders with the names of the LKHs containing the sam/bam files
#Example: python ReadMapping.py
'''

import os
from Bio import SeqIO

#this first command builds your reference with HISAT
#If you've already done this, DON'T run this command! Instead, comment it out (use a # in front of it)
#It will output several files. Don't worry about them, HISAT will know what to do
os.system("hisat2-build Foram_reference.fasta Foram_Index") #change to your reference.fasta and (optionally) rename the index

folder = os.listdir("TrimmedReads") #Insert the name of the folder which has your trimmed reads inside the quotes
folder.sort() #This sorts the folder so that all the LKHs are in order

for x in folder:
	if "LKH" in x and "FPE" in x: #assigning a variable to forward reads. Make sure you have both forward and reverse reads for each cell!
		FPE = x
		sample_id = FPE.split("_FPE")[0]
	if "LKH" in x and "RPE" in x: #assigning a variable to reverse reads
		RPE = x

		if FPE.split("_FPE")[0] == RPE.split("_RPE")[0]: # Match sample IDs dynamically 
			#The next few lines are several HISAT commands that will create new files
			#EDIT the name of the index and the name of the trimmed reads folder in the first command below (line 33) if necessary
			os.system("hisat2 -x Foram_Index -1 TrimmedReads/" +FPE+ " -2 TrimmedReads/" +RPE+ " -S sample.sam") 
			os.system("samtools view -bS sample.sam > sample.bam")
			os.remove("sample.sam")
			os.system("samtools fixmate -O bam sample.bam  fixmate_sample.bam")
			os.system("samtools sort -O bam -o sorted_sample.bam fixmate_sample.bam")
			os.system("sambamba markdup -r sorted_sample.bam sorted_sample.dedup.bam")
			os.system("samtools view -h -b -q 40 sorted_sample.dedup.bam > sorted_sample.q40.bam")
			os.system("samtools view -h -b -q 20 sorted_sample.dedup.bam > sorted_sample.q20.bam")
			os.system("samtools view -h -F 4 -b sorted_sample.dedup.bam > defaultparameters_sample.bam")

			if not os.path.isdir(sample_id):
				os.mkdir(sample_id) #making folders with the names of the LKHs
				
			for file in os.listdir('.'): #These lines move the sam/bam files that Hisat creates into the new LKH folders.
				if(file.endswith('.sam') or file.endswith('.bam')):
					os.rename(file, f"{sample_id}/{file}")
				
print("~~~~~~~~~~~:>~") #When the snake appears in terminal, your script has finished running!
		
	
	
