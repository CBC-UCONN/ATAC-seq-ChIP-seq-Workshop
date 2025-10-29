#!/bin/bash
#SBATCH --job-name=10_frag_dist
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-10

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

module load picard/3.1.1

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
indir=../results/08_filter
outdir=../results/10_frag_dist

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

java -Xmx31G -jar $PICARD CollectInsertSizeMetrics \
-I $indir/$sample.filtered.sorted.bam \
-O $outdir/$sample.fraglen.stats \
-H $outdir/$sample.fraglen.pdf \
-M 0.5


echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"


