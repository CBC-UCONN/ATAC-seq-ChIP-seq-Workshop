#!/bin/bash
#SBATCH --job-name=07_markdups
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-21

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load requrired modules
module load GATK/4.3.0.0

# Store some paths as variables
meta_data=../meta/chip-sra-meta.csv
indir=../results/06_map
outdir=../results/07_markdups

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Get library layout (SINGLE or PAIRED) for the sample
layout=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="LibraryLayout")col=i}NR==row+1&&col{print $col}' $meta_data)

# Run MarkDuplicates on the mapped BAM file
gatk MarkDuplicates \
  --INPUT $indir/${sample}.sorted.bam \
  --OUTPUT $outdir/${sample}.marked.bam \
  --METRICS_FILE $outdir/${sample}.metrics.txt

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"


