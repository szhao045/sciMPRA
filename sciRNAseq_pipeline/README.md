# sci-RNA-seq pipeline
This pipeline is adapted from Junyue Cao's sci-RNA-seq3 pipeline. 
## Dependencies
1. miniconda3 (locally installed)
2. fastqc (module load)
3. R (module load)
    * BiocParallelS
    * tidyverse
    * data.table
    * Matrix
4. samtools/1.12 (module load)
5. star/2.5.2b (module load)
6. cutadapt/1.9.1 (module load)
7. trim_galore (locally installed)
    * Note: copy-paste the path of your local trim_galore into scirpts/sci3_trim.sh

## Preparations
1. Download sciRNAseq_pipeline folder. Always keep sciRNAseq_pipeline/ as the current working directory.
2. create "fastq" folder under sciRNAseq_pipeline and put all fastq.gz files there (note: suffix of fastq files must be .fastq.gz)
3. put lig and RT barcode files under scripts folder
    * sci_main.sh only takes pickle files as barcode input, if you only have txt files, run generate_pickle.py to get pickle file.
4. put sample ID file under sciRMAseq_pipeline/


## Create conda env
```
create env -n sciMPRA
```



## Generate pickle file for RT and ligation barcodes (optional)

Use generate_pickle.py to generate pickle files from txt files. 
You can skip this step if you already have lig_384_bc.pickle2 and RT_384_bc.pickle2 files under sciRNAseq_pipeline/scripts folder

```
conda activate sciMPRA
python3 scripts/generate_pickle.py [ligation barcode txt file directory] [RT barcode txt file directory]
```

## Run sci_main.sh (barcode dissection + STAR mapping + deduplication + gene counting)

### Modify sci_main.sh 
1. Define all the paramaters in sci_main.sh. 
2. Check sci3_trim.sh file, make sure you have right path to trim_galore tool. 

### Run sci_main.sh under sciRNAseq_pipeline/
1. Write a sbatch file, specify the number of cores you want (should be comparable to core number you write in sci_main.sh). The command to run sci_main.sh is:
```
bash sci_main.sh
```
2. submit sbatch job

The pipeline will generate barcode-annotated fastq files, aligned and filtered SAM files and gene count reports. 

## Generate gene count sparse matrix
Run gene_count_processing_sciRNAseq.R to generate Rdata that can be directly loaded in downstream scRNA-seq analysis. 

Script to run this script:

```
module load R
Rscript gene_count_processing_sciRNAseq.R
```
