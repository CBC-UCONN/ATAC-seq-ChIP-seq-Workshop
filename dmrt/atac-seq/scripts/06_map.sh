#!/bin/bash
#SBATCH --job-name=06_map
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-16

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load requrired modules
module load bwa-mem2/2.2.1
module load samtools/1.20

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
indir=../results/04_trim
outdir=../results/06_map
reference=../../../resources/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa.gz

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

  bwa-mem2 mem \
    -t $SLURM_CPUS_PER_TASK \
    $reference \
    $indir/${sample}_1.trimmed.fastq.gz \
    $indir/${sample}_2.trimmed.fastq.gz | samtools sort -o $outdir/${sample}.sorted.bam

elif [ "$layout" == "SINGLE" ]; then

  bwa-mem2 mem \
    -t $SLURM_CPUS_PER_TASK \
    $reference \
    $indir/${sample}_1.trimmed.fastq.gz | samtools sort -o $outdir/${sample}.sorted.bam

else 
  echo "Error: Unknown library layout '$layout' for sample '$sample'"
  exit 1
fi

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"