
input_folder=$1
sample=$2
output_folder=$3

# The original version load python 2 here

echo Trimming sample: $sample
/home/yawei.wu/tools/TrimGalore-0.6.6/trim_galore $input_folder/$sample*.gz -a AAAAAAAA --three_prime_clip_R1 1 -o $output_folder

echo Trimming $sample done.