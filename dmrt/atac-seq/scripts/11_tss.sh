#!/bin/bash
#SBATCH --job-name=11_tss
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-16

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

module load R/4.2.2

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
indir=../results/08_filter
outdir=../results/11_frag_dist

# Create output directory if it doesn't exist 
mkdir -p $outdir

sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)


Rscript 11.tss.R \
  ${indir}/${sample}_filtered.bam \
  ${outdir}



echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"


