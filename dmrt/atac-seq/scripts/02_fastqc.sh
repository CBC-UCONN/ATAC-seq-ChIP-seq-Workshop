#!/bin/bash
#SBATCH --job-name=02_fastqc
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-16

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load fastqc/0.12.1

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
indir=../data/raw-fastq
outdir=../results/02_fastqc

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Get library layout (SINGLE or PAIRED) for the sample
layout=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="LibraryLayout")col=i}NR==row+1&&col{print $col}' $meta_data)

echo $layout

# Run FastQC on the sample reads
if [ "$layout" == "PAIRED" ]; then

  fastqc -o $outdir $indir/${sample}_1.fastq.gz $indir/${sample}_2.fastq.gz

elif [ "$layout" == "SINGLE" ]; then

  fastqc -o $outdir $indir/${sample}_1.fastq.gz

else 
  echo "Error: Unknown library layout '$layout' for sample '$sample'"
  exit 1
fi


echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"