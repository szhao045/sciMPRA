input_folder=$1
sample=$2
output_folder=$3
mismatch=$4

python_script="/scratch/bclab/yawei/sciMPRA/sciMPRA_pipeline/scripts/rm_dup_barcode_UMI_no_mismatch.py"

echo Filtering sample: $sample

python3 $python_script $input_folder/$sample.sam $output_folder/$sample.sam $mismatch

echo Filtering $sample finished