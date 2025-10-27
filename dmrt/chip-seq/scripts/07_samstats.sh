#!/bin/bash
#SBATCH --job-name=07_samstats
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-12

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load samtools/1.20

# Store some paths as variables
meta_data=../meta/chip-sra-meta.csv
indir=../results/06_map
outdir=../results/07_samstats

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Run samstats on the mapped BAM file 
samtools stats $indir/$sample.sorted.bam > $outdir/$sample.samstats.txt

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"