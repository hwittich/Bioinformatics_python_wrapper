#Python wrapper for Henry Wittich's COMPBio mini project
'''This python script does performs a series of analyses on HCMV transcriptome data taken from 2 different donors at 2 different time points: 2 days and 6 days post-infection'''

import argparse
import sys
import os
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO

def check_arg(args=None):
	'''This function parses through the arguments included when running
	this script from command line.'''
	parser = argparse.ArgumentParser(description='Script to perform analysis on HCMV transcriptomes.')
	parser.add_argument('--test','-t',action='store_true',help='Run flag if you want to run the program on test data. Default is false and the program will download the full dataset.')
	return parser.parse_args(args)

#Retrieve command line arguments
args=check_arg(sys.argv[1:])
test=args.test

#Setting up Entrez
Entrez.email = "hwittich@luc.edu"
#Retrieving directories
root = os.getcwd() #root directory
data = root+"/data" #data stored here
#Making list of sample IDs for later use
sample_IDs = ['SRR5660030','SRR5660033','SRR5660044','SRR5660045']

#Retrieving the data
if (test == False): #Only if it isn't a test run #Need to download data from online
   	#Creating the output directory
	os.system("mkdir "+root+"/miniProject_Henry_Wittich")
	outdir = root+"/miniProject_Henry_Wittich"
	#Opening the output file
	outfile = open(outdir+"/miniProject.log","w")

	#Move into data directory for downloads
	os.chdir("data")
	sample_paths = {} #empty dictionary to store the paths to each sample file
	for sample in sample_IDs: #loop through list of IDs to retrieve data
		os.system("wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-11/"+sample+"/"+sample+".1")
		#Making paired-end fastq files
		fastq_dump_command = "fastq-dump --split-files "+sample+".1 -O "+data #store paired-end fastq files in data directory
		os.system(fastq_dump_command)
		sample_paths[sample] = [data+"/"+sample+".1_1.fastq",data+"/"+sample+".1_2.fastq"]

	os.chdir(root) #move back to root

else: #If it is a test run, data is already stored in data directory
	#Creating the output directory
	os.system("mkdir "+root+"/test_outputs")
	outdir = root+"/test_outputs"
	#Opening the output file
	outfile = open(outdir+"/miniProject.log","w")

	#Defining the paths to the samples
	sample_paths = {}
	for sample in sample_IDs:
        	sample_paths[sample] = [data+"/"+sample+".1_1_test.fastq",data+"/"+sample+".1_2_test.fastq"]

#Retrieving the reference transcriptome for kallisto
os.chdir("data")
CDS_count = 0 #store the number of coding sequences
handle = Entrez.efetch(db="nucleotide",id="EF999921",rettype="gb",retmode="text")
record=SeqIO.read(handle,"genbank") #read in genbank entry
handle.close()
transcriptome_path = data+"/HCMV_transcriptome.fa" #path to the transcriptome file
transcriptome = open(transcriptome_path,"w") #open file to write all CDS
for feature in record.features: #iterate through list of features
        if feature.type == "CDS":
                CDS_count+=1
                transcriptome.write(">"+str(record.id)+"\n") #write the sequence ID
                transcriptome.write(str(feature.extract(record.seq))+"\n") #write the coding sequence
transcriptome.close()
#Print number of CDS to outfile file
outfile.write("The HCMV genome (EF999921) has "+str(CDS_count)+" CDS.\n")

#Retrieving the reference genome for bowtie2
handle=Entrez.efetch(db="nucleotide",id="EF999921",rettype="fasta",retmode="text")
record=SeqIO.read(handle,"fasta")
handle.close()

genome_path = data+"/HCMV_genome.fa" #path to the genome file
genome = open(genome_path,"w")
genome.write(">"+str(record.id)+"\n")
genome.write(str(record.seq)+"\n")
genome.close()

#Unzip the Betaherpesvirinae sequences file
os.system("gunzip sequence.fasta.gz")
subfamily_path = data+"/sequence.fasta"

os.chdir(root) #move out of data directory to perform rest of analysis

#Building kallisto index
os.system("mkdir "+outdir+"/kallisto") #making kallisto output directory
kallisto_index_command = "kallisto index -i "+outdir+"/kallisto/HCMV_transcriptome.idx --make-unique "+transcriptome_path
os.system(kallisto_index_command)

#Loop through list of samples and run quantification step
#Using deafult k-mer size, 31, because reads are minimum length 50
sleuth_infile = open(outdir+"/kallisto/sample_table.tsv","w") #file to write table for sleuth
sleuth_infile.write("sample\tdays post-infection\tpath\n") #headers
for sample in sample_IDs:
	kallisto_quantification_command="kallisto quant -i "+outdir+"/kallisto/HCMV_transcriptome.idx -o "+outdir+"/kallisto/"+sample+" -b 30 -t 2 "+sample_paths[sample][0]+" "+sample_paths[sample][1]
	os.system(kallisto_quantification_command)
	if sample == sample_IDs[0] or sample == sample_IDs[2]: #1st and 3rd samples are 2 days post-infection
		sleuth_infile.write(sample+"\t2\t"+outdir+"/kallisto/"+sample+"\n")
	else: #Other samples are 6 days post-infection
		sleuth_infile.write(sample+"\t6\t"+outdir+"/kallisto/"+sample+"\n")
sleuth_infile.close()

outfile.close() #close output file so that R can write into file
#Analyzing the kallisto quantification output with sleuth in R
os.system("Rscript sleuth_analysis.R "+str(test)) #Adding test argument so R knows where to find the output file
outfile = open(outdir+"/miniProject.log","a") #Reopen output file for Python to write to
#append mode so I don't overwrite the previous lines

#Creating a bowtie2 index from the HCMV genome
os.system("mkdir "+outdir+"/bowtie") #making bowtie output directory
bowtie_build_command = "bowtie2-build "+genome_path+" "+outdir+"/bowtie/HCMV"
os.system(bowtie_build_command)

#Aligning transcriptomes to HCMV reference genome with bowtie2 to filter out host transcript reads
for sample in sample_IDs: #loop through samples
        #Calculate number of read pairs before filtering
	before = len(open(sample_paths[sample][0],"r").readlines())/4 #each read in fastq format takes up 4 lines

	#Perform bowtie2 alignment, storing the mapped reads in a new fastq file
	bowtie_align_command = "bowtie2 --quiet -x "+outdir+"/bowtie/HCMV -1 "+sample_paths[sample][0]+" -2 "+sample_paths[sample][1]+" -S "+outdir+"/bowtie/HCMVmap.sam --al-conc "+outdir+"/bowtie/"+sample+".1_mapped_%.fq"
	os.system(bowtie_align_command)

	#Add the mapped reads to sample paths dictionary
	sample_paths[sample].append(outdir+"/bowtie/"+sample+".1_mapped_1.fq") #at index 2
	sample_paths[sample].append(outdir+"/bowtie/"+sample+".1_mapped_2.fq") #at index 3

	#Calculate number of read pairs after filtering
	after = len(open(sample_paths[sample][2],"r").readlines())/4

        #Writing to output file
	if sample == sample_IDs[0] or sample == sample_IDs[1]: #1st and 2st samples are from donor 1
		if sample == sample_IDs[0]: #first sample is 2dpi
			outfile.write("Donor 1 (2dpi) had "+str(before)+" read pairs before Bowtie2 filtering and "+str(after)+" read pairs after.\n")
		else: #2nd sample is 6dpi
			outfile.write("Donor 2 (6dpi) had "+str(before)+" read pairs before Bowtie2 filtering and "+str(after)+" read pairs after.\n")
	else: #3rd and 4th samples are from donor 3
		if sample == sample_IDs[2]: #3rd sample is 2dpi
			outfile.write("Donor 3 (2dpi) had "+str(before)+" read pairs before Bowtie2 filtering and "+str(after)+" read pairs after.\n")
		else: #4th sample is 6dpi
			outfile.write("Donor 4 (6dpi) had "+str(before)+" read pairs before Bowtie2 filtering and "+str(after)+" read pairs after.\n")

#Assembling all four transcriptomes with SPAdes
os.system("mkdir "+outdir+"/SPAdes") #making SPAdes output directory
SPAdes_command = "spades -k 55,77,99,127 -t 2 --only-assembler --pe1-1 "+sample_paths[sample_IDs[0]][2]+" --pe1-2 "+sample_paths[sample_IDs[0]][3]+" --pe2-1 "+sample_paths[sample_IDs[1]][2]+" --pe2-2 "+sample_paths[sample_IDs[1]][3]+" --pe3-1 "+sample_paths[sample_IDs[2]][2]+" --pe3-2 "+sample_paths[sample_IDs[2]][3]+" --pe4-1 "+sample_paths[sample_IDs[3]][2]+" --pe4-2 "+sample_paths[sample_IDs[3]][3]+" -o "+outdir+"/SPAdes/HCMV_assembly"
os.system(SPAdes_command)
outfile.write(SPAdes_command+"\n")

#Evaluating the assembly
contigs_file = open(outdir+"/SPAdes/HCMV_assembly/contigs.fasta") #making a handle for SeqIO to parse
longest_contig = Seq("") #finding the longest_contig for BLAST
assembly_length = 0
long_contigs = 0 #counting contigs >1000bp
#Parsing the contigs file
records = SeqIO.parse(contigs_file,"fasta")
#Looping through the records
for record in records:
	assembly_length = assembly_length+len(record.seq) #add sequence length to total
	if len(record.seq) > 1000:
		long_contigs+=1
	if len(record.seq) > len(longest_contig):
		longest_contig=record.seq
outfile.write("There are "+str(long_contigs)+" contigs > 1000 bp in the assembly.\n")
outfile.write("There are "+str(assembly_length)+" bp in the assembly.\n")
contigs_file.close()

#Generate the query file for BLAST with the longest contig
os.system("mkdir "+outdir+"/BLAST") #making BLAST output directory
query_file = open(outdir+"/BLAST/query.fasta","w")
query_file.write(">HCMV\n")
query_file.write(str(longest_contig)+"\n")
query_file.close()

#Creating the BLAST db
make_blast_command = "makeblastdb -in "+subfamily_path+" -out "+outdir+"/BLAST/Betaherpesvirinae -title Betaherpesvirinae -dbtype nucl"
os.system(make_blast_command)

#Running BLAST
blast_command = 'blastn -query '+outdir+'/BLAST/query.fasta -db '+outdir+'/BLAST/Betaherpesvirinae -out '+outdir+'/BLAST/myresults.csv -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"'
os.system(blast_command)

#Parsing the BLAST output
outfile.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")
blast_output = open(outdir+"/BLAST/myresults.csv","r")
lines = blast_output.readlines()
for i in range(0,10): #Pull the top 10 hits from the top of the file
	line = lines[i][:-1].split(",")
	for i in range(0,10): #We want to write the 10 attributes to the output file
		outfile.write(str(line[i])+"\t")
	outfile.write("\n")

outfile.close()
