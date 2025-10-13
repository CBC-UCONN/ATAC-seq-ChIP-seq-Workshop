#!/bin/bash
#SBATCH --job-name=04_trim
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-16

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load requrired modules
module load fastp/0.23.2 

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
indir=../data/raw-fastq
outdir=../results/04_trim

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

  fastp \
    -i $indir/${sample}_1.fastq.gz \
    -I $indir/${sample}_2.fastq.gz \
    -o $outdir/${sample}_1.trimmed.fastq.gz \
    -O $outdir/${sample}_2.trimmed.fastq.gz \
    -h $outdir/${sample}_fastp.html \
    -j $outdir/${sample}_fastp.json

elif [ "$layout" == "SINGLE" ]; then

  fastp \
    -i $indir/${sample}_1.fastq.gz \
    -o $outdir/${sample}_1.trimmed.fastq.gz \
    -h $outdir/${sample}_fastp.html \
    -j $outdir/${sample}_fastp.json

else 
  echo "Error: Unknown library layout '$layout' for sample '$sample'"
  exit 1
fi

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"