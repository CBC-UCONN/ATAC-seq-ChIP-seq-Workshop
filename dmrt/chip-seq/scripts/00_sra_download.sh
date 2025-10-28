#!/bin/bash
#SBATCH --job-name=00_download
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-21

# Please do not run this script in the workshop.
# Instead just run 01_symlink.sh to symlink to the pre-downloaded data

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load sratoolkit/3.0.5

# Store some paths as variables
meta_data=../meta/chip-sra-meta.csv
outdir=../data/raw-fastq
tmp=/scratch/$USER/sra-cache

mkdir -p $tmp

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Download the sample fastq files using the SRA toolkit
fasterq-dump $sample --temp $tmp --split-files --outdir $outdir 

# Compress the downloaded fastq files to save space
gzip $outdir/${sample}_*.fastq

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"