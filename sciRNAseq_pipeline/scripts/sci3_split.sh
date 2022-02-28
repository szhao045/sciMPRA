
input_folder=$1
sample=$2
output_folder=$3
barcode_file=$4
cutoff=$5

mismatch=1

python_script="/scratch/bclab/yawei/sciMPRA/sciMPRA_pipeline/scripts/sam_split.py"

python3 $python_script $input_folder/$sample.sam $barcode_file $output_folder $cutoff
echo splitting sample done: $sample