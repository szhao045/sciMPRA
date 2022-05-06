#!/bin/bash
###################### Define general parameters #####################

# define the fastq folder including all the fastq files
fastq_folder="/scratch/bclab/yawei/sciMPRA/sciMPRA_pipeline/test_fastq"

# define the PCR group sample id for each fastq file
# A list of simplified sample names for fastq files
sample_ID="./sample_ID.txt"

# define the output folder
all_output_folder="${fastq_folder}/../output"

# define the script locations
script_folder="./scripts"


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
samtools_core=4
# define the location of the R script for multi-core processing
R_script=$script_folder/sci3_bash_input_ID_output_core.R
script_path=$script_folder


# define the number of UMI cutoff for splitting single cell; 
# cells with UMIs less than this number will be discarded (this is a filtering parameter)
cutoff=200

# define the mismatch rate for removing duplicates:
# If two reads only have one base pair mismatch, they are treated as duplicates
mismatch=1

# define the location of index file for read alignment for STAR
index="/scratch/ref/star/star_genome_d1_vd1_gtfv22"

# define the gtf file for gene counting
gtf_file="/scratch/bclab/yawei/ref/gtf_reference/mm10/gencode.vM12.chr_patch_hapl_scaff.annotation.gtf.gz"


# load required modules
spack load samtools/1.12
spack load star/2.5.2b
### Note star version needs to be compatible with genome index
spack load fastqc
# manually install trimgalore in /home/yawei.wu/tools/TrimGalore-0.6.6/trim-galore
# Note: trimgalore is a wrapper of fastqc and cutadapt
spack load r
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
for sample in $(cat $sample_ID); do echo changing name $sample; mv $input_folder/*$sample*R1*.fastq.gz $input_folder/$sample.R1.fastq.gz; mv $input_folder/*$sample*R2*.fastq.gz $input_folder/$sample.R2.fastq.gz; done

echo "Attaching barcode and UMI...."
mkdir -p $output_folder
python3 $script $input_folder $sample_ID $output_folder $ligation_barcode $RT_barcode $core
echo "Barcode transformed and UMI attached."

########################## Trimming read2 #########################
# the script take UMI attached R2 and output trimmed fastq files
echo
echo "Start trimming the read2 file..."
echo $(date)

trimmed_fastq=$all_output_folder/trimmed_fastq
UMI_attached_R2=$all_output_folder/UMI_attach
bash_script=$script_path/sci3_trim.sh
### Rscript needs to have bplapply
Rscript $R_script $bash_script $UMI_attached_R2 $sample_ID $trimmed_fastq $core

########################## STAR mapping  #######################
#align the reads with STAR, 
#define the output folder for mapping
echo
input_folder=$trimmed_fastq
STAR_output_folder=$all_output_folder/STAR_alignment
filtered_sam_folder=$all_output_folder/filtered_sam
rmdup_sam_folder=$all_output_folder/rmdup_sam

#align read2 to the index file using STAR
### Need to figure out what are the output files
echo "Start alignment using STAR..."
echo input folder: $input_folder
echo sample ID file: $sample_ID
echo index file: $index
echo output_folder: $STAR_output_folder
#make the output folder
mkdir -p $STAR_output_folder
#load the genome
STAR --genomeDir $index --genomeLoad LoadAndExit
#start the alignment
for sample in $(cat $sample_ID); do echo Aligning $sample;STAR --runThreadN $core --outSAMstrandField intronMotif --genomeDir $index --readFilesCommand zcat --readFilesIn $input_folder/$sample*gz --outFileNamePrefix $STAR_output_folder/$sample --genomeLoad LoadAndKeep; done
#remove the index from the memory
STAR --genomeDir $index --genomeLoad Remove
echo "All alignment done."

### Note: only 10% of the reads uniquely mapped after this step
########################## Filtering&Sorting SAM files #######################
#Filter sam file based on q > 30, and sort the files 

echo
echo "Start filtering and sorting the sam files..."
echo input folder: $STAR_output_folder
echo output folder: $filtered_sam_folder
bash_script=$script_path/sci3_filter.sh
### Again they call bash script from R
Rscript $R_script $bash_script $STAR_output_folder $sample_ID $filtered_sam_folder $samtools_core

####################### Removing duplicates ###################################
# remove duplicates based on UMI sequence and tagmentation site

echo
echo "Start removing duplicates..."
echo input folder: $filtered_sam_folder
echo output folder: $rmdup_sam_folder
mkdir -p $rmdup_sam_folder

# remove exactly matched duplicates
bash_script=$script_path/sci3_rmdup_nomismatch.sh 
## bash_script=$script_path/sci3_rmdup.sh
Rscript $R_script $bash_script $filtered_sam_folder $sample_ID $rmdup_sam_folder $core $mismatch

# repeat the rmdup process to remove duplicates based on edit distance of UMI sequence
echo
echo "Start removing duplicates..."
echo input folder: $all_output_folder/rmdup_sam
echo output folder: $all_output_folder/rmdup_sam_2
mkdir -p $all_output_folder/rmdup_sam_2

bash_script=$script_path/sci3_rmdup.sh
filtered_sam_folder=$all_output_folder/rmdup_sam
rmdup_sam_folder=$all_output_folder/rmdup_sam_2
Rscript $R_script $bash_script $filtered_sam_folder $sample_ID $rmdup_sam_folder $core $mismatch

######################### demultiplexing sam files #########################
# split the sam file based on the barcode, and mv the result to the report folder
sam_folder=$all_output_folder/rmdup_sam_2
output_folder=$all_output_folder/sam_splitted

echo
echo "Start splitting the sam file..."
echo samfile folder: $sam_folder
echo sample list: $sample_ID
echo ouput folder: $output_folder
echo barcode file: $barcodes
echo cutoff value: $cutoff


bash_script=$script_path/sci3_split.sh
Rscript $R_script $bash_script $sam_folder $sample_ID $output_folder $core $barcodes $cutoff


cat $output_folder/*sample_list.txt>$output_folder/All_samples.txt
cp $output_folder/All_samples.txt $output_folder/../barcode_samples.txt
# output the report the report/barcode_read_distribution folder
mkdir -p $output_folder/../report/barcode_read_distribution
mv $output_folder/*.txt $output_folder/../report/barcode_read_distribution/
mv $output_folder/*.png $output_folder/../report/barcode_read_distribution/
echo
echo "All sam file splitted."

########################### counting genes #############################
# count reads mapping to genes
output_folder=$all_output_folder/report/human_mouse_gene_count/
input_folder=$all_output_folder/sam_splitted
script=$script_path/sciRNAseq_count.py
sample_ID=$all_output_folder/barcode_samples.txt
echo "Start the gene count...."
python3 $script $gtf_file $input_folder $sample_ID $core

echo "Make the output folder and transfer the files..."
mkdir -p $output_folder
cat $input_folder/*.count > $output_folder/count.MM
rm $input_folder/*.count
cat $input_folder/*.report > $output_folder/report.MM
rm $input_folder/*.report
mv $input_folder/*_annotate.txt $output_folder/
echo "All output files are transferred~"

now=$(date)
echo "Current time : $now"

