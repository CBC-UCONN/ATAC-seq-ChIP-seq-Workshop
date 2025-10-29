#!/bin/bash
#SBATCH --job-name=14_counts
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-10

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

module load subread/2.0.3

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
bamdir=../results/08_filter
peakdir=../results/12_macs
outdir=../results/14_counts

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Convert peak file to SAF format
peakfile=$peakdir/${sample}_peaks.broadPeak
saf=$outdir/${sample}.saf
awk 'OFS="\t" {print $1"."$2+1"."$3, $1, $2+1, $3, "."}' $peakfile > $saf 

# Count reads in peaks
featureCounts \
  -p \
  -a $saf \
  -F SAF \
  -o $outdir/${sample}.counts \
  $bamdir/$sample.filtered.sorted.bam

rm $saf

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"

