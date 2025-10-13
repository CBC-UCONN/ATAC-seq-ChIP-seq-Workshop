#!/bin/bash
#SBATCH --job-name=03_multiqc_raw
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%j.out


echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load requrired modules
module load MultiQC/1.9

# Store some paths as variables
indir=../results/02_fastqc
outdir=../results/03_multiqc_raw

# Create output directory if it doesn't exist
mkdir -p $outdir

# Run MultiQC
multiqc -o $outdir $indir

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"