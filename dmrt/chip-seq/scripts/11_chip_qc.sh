#!/bin/bash
#SBATCH --job-name=11_chip_qc
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-12

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

source ~/.bashrc
mamba activate r_atac

# Store some paths as variables
meta_data=../meta/chip-sra-meta.csv
filtered=../results/08_filter
peaks=../results/10_macs

# Create output directory if it doesn't exist 

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

outdir=../results/11_chip_qc/$sample
mkdir -p $outdir

# Run QC R script
Rscript 11_chip_qc.R \
  $filtered/$sample.filtered.sorted.bam \
  $peaks/${sample}_peaks.narrowPeak \
  $outdir


echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"
