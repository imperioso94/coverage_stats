#!/usr/bin/env python
# this program gets statistics about coverage from tab delimited files created using "bedtools genomecov"
# bedtools genomecov -bga -ibam example_1.bam > example_1.coverage.bed\n\n\t and than:\n\t\tgrep -w \"0$\" example_1.coverage.bed > 0_cov_file.bed\n\n"
import os
import sys, getopt
import csv
import Bio.SeqIO as IO
import math



stdout_fileno = sys.stdout
file=""
genome=""


if len(sys.argv)==5 :
    
    try:
        opts, args = getopt.getopt(sys.argv[1:],"hi:g:",["input_bam=","genome="])
        
    except getopt.GetoptError:
        stdout_fileno.write("\nusage: ./surface_cov_stats.py -i <bam_file> -g <genome.fasta> \n")
        
        sys.exit(2)
else:
    stdout_fileno.write("\nusage: ./surface_cov_stats.py -i <bam_file> -g <genome.fasta> \n")
    
    sys.exit(2)


for opt, arg in opts:     
    if opt in ("-i", "--input"):
        
        file = arg

    elif opt in ("-g", "--genome"):
        
        genome = arg
        
        
        
        
    elif opt in ("-h"):
        stdout_fileno.write("\nusage: ./surface_cov_stats.py -i <bam_file> -g <genome.fasta> \n\n\n\t <bam_file is the file to be analised\n\t\ <genome.fasta> is the reference genome on which reads where alligned")
        sys.exit(2)
        
    else:
        stdout_fileno.write("\nusage: ./surface_cov_stats.py -i <bam_file> -g <genome.fasta> \n")
        
        sys.exit(2)

# create coverage_file and 0_cov_file

os.system("mkdir -p ./coverage")
os.system("mkdir -p ./0_cov")

string_1= "bedtools genomecov -bga -ibam "+file+" > ./coverage/"+ file + ".coverage.bed"
os.system(string_1)

string_2= "grep -w \"0$\" ./coverage/"+file+".coverage.bed > ./0_cov/"+file+".0_cov.bed" 
os.system(string_2)

 

def calculate_surface_cov(no_cov,genome):

    record_dict = IO.to_dict(IO.parse(genome, "fasta"))
    for key in record_dict.items():
        print("\n\t",key[0]," lenght: \n\t",len(key[1].seq),"\n")
        genome_lenght = len(key[1].seq)
    
    total_count=0
     
    with open(no_cov) as no_cov_lines:
        no_cov_reader = csv.reader(no_cov_lines, delimiter='\t')
        for line in no_cov_reader:
            total_count += int(line[2]) - int(line[1])
            
    
    surface_coverage_perc = round(((genome_lenght - total_count) / (genome_lenght))*100,2)

    print("\tsurface coverage percentage: ",surface_coverage_perc,"%\n")


def calculate_medium_stdv_per_base_coverage(coverage,genome):

    record_dict = IO.to_dict(IO.parse(genome, "fasta"))
    genome_lenght = 0
    
    for key in record_dict.items():
        
        genome_lenght = len(key[1].seq)
        
    total_count=0

    # mean calculation
    
    with open(coverage) as coverage_lines:
        file_reader = csv.reader(coverage_lines, delimiter='\t')
        for line in file_reader:
            total_count += (int(line[2]) - int(line[1]))*int(line[3])
            
    mean = round(total_count / genome_lenght, 2)

    count=0
    stdv=0

    # stdv calculation
    
    with open(coverage) as coverage_lines:
        file_reader = csv.reader(coverage_lines, delimiter='\t')
        for line in file_reader:
            count = (int(line[2]) - int(line[1]))
            
            for i in range(count):
                
                stdv += (int(line[3])-mean)**2
                
        stdv = round(math.sqrt(stdv / genome_lenght), 2)


    print("\tmean per base coverage: ", mean)
    print("\tstdv per base coverage: ", stdv,"\n")
        
    
    
file_0_coverage="./0_cov/"+file+".0_cov.bed"
file_coverage="./coverage/"+file+".coverage.bed"

calculate_surface_cov(file_0_coverage,genome)
calculate_medium_stdv_per_base_coverage(file_coverage,genome)
