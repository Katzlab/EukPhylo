#!/usr/bin/env python3

#Author, date: Giulia Magri Ribeiro and Adri K. Grow updated from Xyrus Maurer-Alcala and Ying Yan; June 13 2025
#Motivation: Trim adaptors from reads and quality trimming before Assembly
#Intent: clean up reads
#Dependencies: biopython and bbmap folder
#Inputs: parameters.txt, fastq.gz forward and reverse reads
#Outputs:trimmed reads in ToAssemble folder
#Example: python3 Trim_Reads.py parameter.txt
#Katzlab parameters are 24 for quality trimming and 75 for minimum length as of June 2025


from Bio import SeqIO
import sys,os
import time

#------------------------------ Checks the Input Arguments ------------------------------#

if len(sys.argv) == 1:
	print ('\n\nThis script will remove Adapters, do quality trimming and length trimming on given score and assembly from your raw reads')
	print ('\n\nChecking the overall quality and reads size on FastQC is recommended\n\n')
	print ('Example Usage:\n\n\t' + 'katzlab$ python3 Trim_Reads.py parameter.txt\n\n')
	print ('\t\tQuestions/Comments? Email Giulia (author) at gribeiro@smith.edu\n\n')
	sys.exit()


elif len(sys.argv) != 2:
	print ('\n\nDouble check that you have added all the necessary command-line inputs! (see usage below for an example)\n\n')
	print  ('Example Usage:\n\n\t' + 'katzlab$ python3 Trim_Reads.py parameter.txt\n\n')	
	print ('Please also check that you have a parameter.txt (tab separated values) file which should contain your current filename, new filename, score of quality trimming, and minimum length (see an example below)\n\n')
	print  ('parameter.txt example:\n\n\t' + 'XKATZ_20161110_K00134_IL100076423_S41_L005\tLKH001_Spirostomum\t24\t100\n\tXKATZ_20161110_K00134_IL100076416_S17_L005\tLKH002_Loxodes\t28\t100\n')
	sys.exit()

else:
	parameter_file = sys.argv[1]
	mailaddress = 'your_email@xxx.edu'  # default email
	if not os.path.isdir('ToAssemble/'):
		os.system('mkdir ToAssemble')
	
### takes your downloaded data and renames the file so that it has taxonomic information in the filename
def rename(code):
	for filename in os.listdir(os.curdir):
		if filename.endswith('.fastq.gz'):
		### check name code here for forward reads
			if '_FwdPE' in filename or '_R1' in filename:
				cur_name = filename.split('_FwdPE')[0] if '_FwdPE' in filename else filename.split('_R1')[0]
				if cur_name in code:
					new_name = code[cur_name]
					print(cur_name, new_name)
					os.system(f'mv {filename} {new_name}_FwdPE.fastq.gz')
					os.system(f'mkdir -p {new_name}')
		### check name code here for reverse reads
			elif '_RevPE' in filename or '_R2' in filename:
				cur_name = filename.split('_RevPE')[0] if '_RevPE' in filename else filename.split('_R2')[0]
				if cur_name in code:
					new_name = code[cur_name]
					print(cur_name, new_name)
					os.system(f'mv {filename} {new_name}_RevPE.fastq.gz')


### Uses the adapters.fa file in the bbtools resources folder (and BBDuK) to remove adapter sequences -- update if necessary
### Uses BBDuK to quality trim reads so the average is q24 and the min length is 100 -- adjust if needed ... flags will be added eventually			
def QualityTrim(qtrim, minlen):
	for filename in os.listdir(os.curdir):
		if 'FwdPE' in filename:
			new_name = filename.split('_FwdPE')[0]
			qscore = qtrim[new_name]
			lscore = minlen[new_name]
			qtrimcmd = '_q'+qscore+'_minlen'+lscore
			log_file = filename.split('_Fwd')[0] + '/' + filename.split('_Fwd')[0] + qtrimcmd + '_bbduk.log'
			os.system('./bbmap/bbduk.sh -Xmx20g in1=./' + filename + ' in2=./' + filename.replace('Fwd','Rev') + ' out1=ToAssemble/'+filename.replace('FwdPE','FPE'+qtrimcmd) + ' out2=ToAssemble/' + filename.split('Fwd')[0]+'RPE'+qtrimcmd+'.fastq.gz qtrim=rl trimq='+qscore+' minlen='+lscore+' mink=11 k=23 hdist=1 ktrim=r ref=bbmap/resources/adapters.fa stats=' + filename.split('_Fwd')[0] +'/'+ filename.split('_Fwd')[0] + qtrimcmd + '_Stats.txt overwrite=true'+ ' > ' + log_file + ' 2>&1')	
			

### Calls on rnaSPAdes to do the transcriptome assembly on the quality trimmed files.				
#def rnaSPAdesAssembly():
#	for filename in os.listdir(os.curdir+'/ToAssemble'):						
#		if 'LKH' in filename:
#		if 'FPE_q' in filename:
#			os.system('python rnaSPAdes-0.1.1/bin/rnaspades.py -m 26 -k 21,33,55,77 --min-complete-transcript 300 -1 ToAssemble/' + filename + ' -2 ToAssemble/' + filename.replace('FPE','RPE')+' -o ' + filename.split('_FPE')[0] + '/; echo "Finished assembling ' + filename.split('_FPE')[0] + '" | mail -s "Finished Transcriptome Assembly ' + (time.strftime("%d/%m/%y")) + '" ' + mailaddress) > out.txt
	

def main():
	code = {}
	qtrim = {}
	minlen = {}
	for line in open(parameter_file,'r'):
		code[line.split('\t')[0]] = line.split('\t')[1].split('\n')[0]
		qtrim[line.split('\t')[1]] = line.split('\t')[2].split('\n')[0]
		minlen[line.split('\t')[1]] = line.split('\t')[3].split('\n')[0]
	rename(code)
	QualityTrim(qtrim, minlen)
#	rnaSPAdesAssembly()
main()
