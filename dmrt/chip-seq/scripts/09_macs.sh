#!/bin/bash
#SBATCH --job-name=10_macs
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-21

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load macs3/3.0.0a6 

# Store some paths as variables
meta_data=../meta/chip-sra-meta.csv
indir=../results/06_map
outdir=../results/10_macs

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Get library layout (SINGLE or PAIRED) for the sample
layout=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="LibraryLayout")col=i}NR==row+1&&col{print $col}' $meta_data)

# Run fastp on the sample reads
if [ "$layout" == "PAIRED" ]; then

  macs3 callpeak \
    --treatment $indir/$sample.sorted.bam \
    --format BAMPE \
    --gsize 2654621783 \
    --name $outdir/$sample \
    --nolambda \
    --bdg \
    --trackline \
    --qvalue 0.05

elif [ "$layout" == "SINGLE" ]; then

  macs3 callpeak \
    --treatment $indir/$sample.sorted.bam \
    --format BAM \
    --gsize 2654621783 \
    --name $outdir/$sample \
    --nolambda \
    --bdg \
    --trackline \
    --qvalue 0.05

else 
  echo "Error: Unknown library layout '$layout' for sample '$sample'"
  exit 1
fi

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"

