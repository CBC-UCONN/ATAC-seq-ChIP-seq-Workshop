#!/bin/bash
#SBATCH --job-name=08_multiqc_mapping
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=mcbstudent
#SBATCH --qos=mcbstudent
#SBATCH --output=logs/%x_%j.out

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load MultiQC/1.9

# Store some paths as variables
meta_data=../meta/chip-sra-meta.csv
indir="../results/07_samstats ../results/07_markdups"
outdir=../results/08_multiqc_mapping

# Create output directory if it doesn't exist
mkdir -p $outdir

# Run MultiQC
multiqc -o $outdir $indir

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"