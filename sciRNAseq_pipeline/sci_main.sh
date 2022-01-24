#!/bin/bash

# define the fastq folder including all the fastq files
fastq_folder="/scratch/bclab/szhao/220124_sciRNA_tech/fastq"

# define the PCR group sample id for each fastq file
# Not sure what it does yet
sample_ID="./sample_ID.txt"

# output folder
all_output_folder="${fastq_folder}/output"

# define the core number for parallele processing 
core=15
samtools_core=4

# define the number of UMI cutoff
# right now I am just doing it for the tech dev part
cutoff=10

# define the location of index file for read alignment for STAR
index="/scratch/ref/star/star_genome_d1_vd1_gtfv22"

# define the gtf file for gene counting
gtf_file="/net/shendure/vol1/home/cao1025/reference/gtf_reference/mm10/gencode.vM12.chr_patch_hapl_scaff.annotation.gtf.gz"

#define the mismatch rate for removing duplicates:
mismatch=1

#define the script locations
script_folder="./scripts"

#load modules
module load samtools/1.9
module load star/2.7.2b
module load miniconda3
module load cutadapt/1.9.1

#Activate ocnda environment
eval "$(conda shell.bash hook)"
conda activate mascot




