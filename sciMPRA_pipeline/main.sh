#!/bin/bash
###################### Define general parameters #####################

# define the fastq folder for the raw sequencing data for sciMPRA. 
fastq_folder="~/barcode_fastq/"

# define the PCR group sample id for each fastq file
# A list of simplified sample names for fastq files
sample_ID="./sample_ID.txt"

# define the output folder
all_output_folder="${fastq_folder}/../output"

# define the script locations
script_folder="../sciRNAseq_pipeline/scripts"


# define the location of the RT barcodes
### Note: RT_384_bc.txt is a list of RT barcode sequences
RT_barcode=$script_folder//RT_384_bc.pickle2
# define the location of the ligation barcodes (they are in the script folder)
### Note: lig_384_bc.txt is a list of ligation barcode sequences
ligation_barcode=$script_folder/lig_384_bc.pickle2
# define the location of the combined RT and ligation barcodes
### Note: combined_384_bc.txt contains 384x384 barcode sequences
barcodes=$script_folder//combined_384_bc.txt

# define the core number for parallele processing 
core=4
#samtools_core=4
# define the location of the R script for multi-core processing
#R_script=$script_folder/sci3_bash_input_ID_output_core.R
script_path=$script_folder


# define the number of UMI cutoff for splitting single cell; 
# cells with UMIs less than this number will be discarded (this is a filtering parameter)
#cutoff=200

# define the mismatch rate for removing duplicates:
# If two reads only have one base pair mismatch, they are treated as duplicates
#mismatch=1

# define the location of index file for read alignment for STAR
#index="/scratch/ref/star/star_genome_d1_vd1_gtfv22"

# define the gtf file for gene counting
#gtf_file="/scratch/bclab/yawei/ref/gtf_reference/mm10/gencode.vM12.chr_patch_hapl_scaff.annotation.gtf.gz"


# load required modules
#module load samtools/1.12
### Note star version needs to be compatible with genome index
#module load star/2.5.2b
# module load miniconda3
### Note: I use my own conda
#module load cutadapt/1.9.1
#module load fastqc
# manually install trimgalore in /home/yawei.wu/tools/TrimGalore-0.6.6/trim-galore
# Note: trimgalore is a wrapper of fastqc and cutadapt
#module load r
### Note: BiocParallel is required in R

# Activate conda environment
eval "$(conda shell.bash hook)"
# need to figure out which packages to install for this env
conda activate sciMPRA

now=$(date)
echo "Current time : $now"
######################### UMI attach ###################################
# the script take an input folder, a sample ID list, an output folder, the RT barcode list, the ligation barcode list and core number. Then it extract the RT and ligation barcode from read1, correct them to the nearest RT and ligation barcode (with edit distance <= 1), and attach the RT and ligation barcode and UMI sequence to the read name of read2. Reads with unmatched RT or ligation barcodes are discarded.
# Note: the current barcode mapping script does not do fuzzy match, you cannot really change the mismatch paramter.
# Note: this function outputs a modified read2 fastq file with barcode information in the first lines, separated by ","

input_folder=$fastq_folder
output_folder=$all_output_folder/UMI_attach
script=$script_path/UMI_barcode_attach_gzipped_with_dic.py
# simplify file names
#echo "Changing the name of the fastq files..."
for sample in $(cat $sample_ID); do echo changing name $sample; mv $input_folder/*$sample*1.fastq.gz $input_folder/$sample.R1.fastq.gz; mv $input_folder/*$sample*2.fastq.gz $input_folder/$sample.R2.fastq.gz; done

echo "Attaching barcode and UMI...."
mkdir -p $output_folder
python3 $script $input_folder $sample_ID $output_folder $ligation_barcode $RT_barcode $core
echo "Barcode transformed and UMI attached."


######################### Find quads ###################################
# here we find quads for each sample with the slower python script 
# TODO: implement the faster GO script for finding quads.
echo "Find quads ..."
UMI_attached_R2=$all_output_folder/UMI_attach
# Use the original R1 fastq file
R1_folder = $input_folder
output_folder = $all_output_folder/quads
python3 ./script/find_quads.py $UMI_attached_R2 $R1_folder $output_folder $sample_ID $core
