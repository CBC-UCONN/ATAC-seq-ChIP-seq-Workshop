#!/bin/bash
#SBATCH --job-name=11_atac_qc
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-10

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

source ~/.bashrc
mamba activate atac

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
unfiltered=../results/06_map
filtered=../results/08_filter
outdir=../results/11_atac_qc

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Run QC R script
Rscript 11_atac_qc.R \
  $unfiltered/$sample.sorted.bam \
  $filtered/$sample.filtered.sorted.bam \
  $outdir

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"


