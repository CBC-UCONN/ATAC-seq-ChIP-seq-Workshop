#!/bin/bash
#SBATCH --job-name=12_motif
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-12

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load homer/4.11

# Store some paths as variables
meta_data=../meta/chip-sra-meta.csv
indir=../results/10_macs
outdir=../results/12_motif
ref=../../../resources/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

findMotifsGenome.pl \
  $indir/${sample}_peaks.narrowPeak \
  $ref \
  $outdir/${sample}_motifs \
  -size 50 \
  -S 6


echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"

