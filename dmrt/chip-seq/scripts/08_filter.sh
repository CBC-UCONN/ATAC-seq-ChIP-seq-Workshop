#!/bin/bash
#SBATCH --job-name=08_filter
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=mcbstudent
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-12

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load deeptools/3.5.0
module load samtools/1.20

export TMPDIR=/scratch/$USER/deeptools_tmp
mkdir -p $TMPDIR

# Store some paths as variables
meta_data=../meta/chip-sra-meta.csv
indir=../results/06_map
outdir=../results/08_filter

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Filter BAM file using deeptools alignmentSieve
alignmentSieve \
  -b $indir/$sample.sorted.bam \
  -o $outdir/$sample.filtered.bam \
  --filterMetrics $outdir/$sample.filtered.stats \
  --minMappingQuality 30 \
  --ignoreDuplicates

# Create index for filtered BAM
samtools sort -o $outdir/$sample.filtered.sorted.bam $outdir/$sample.filtered.bam
samtools index $outdir/$sample.filtered.sorted.bam

rm $outdir/$sample.filtered.bam

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"