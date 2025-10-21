#!/bin/bash
#SBATCH --job-name=08_filter
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-16

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
# module load GATK/4.3.0.0 
# module load samtools/1.20
module load deeptools/3.5.0

export TMPDIR=/scratch/$USER/deeptools_tmp
mkdir -p $TMPDIR

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
indir=../results/06_map
outdir=../results/08_filter

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# # Get library layout (SINGLE or PAIRED) for the sample
# layout=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
#     'NR==1{for(i=1;i<=NF;i++)if($i=="LibraryLayout")col=i}NR==row+1&&col{print $col}' $meta_data)

alignmentSieve \
  -b $indir/$sample.sorted.bam \
  -o $outdir/$sample.filtered.bam \
  --filterMetrics $outdir/$sample.filtered.stats \
  --minMappingQuality 30 \
  --ignoreDuplicates \
  --ATACshift 



# gatk MarkDuplicates \
#   --INPUT $indir/$sample.sorted.bam \
#   --OUTPUT $outdir/$sample.dedup.tmp.bam \
#   --METRICS_FILE $outdir/$sample.dedup_stats.txt \
#   --REMOVE_DUPLICATES true


# if [ "$layout" == "PAIRED" ]; then

#   samtools view -h $outdir/$sample.dedup.tmp.bam | grep -v chrM \
#   | samtools view -q 30 -b -h -f 2 \
#   | samtools sort -o $outdir/$sample.filtered.bam

# elif [ "$layout" == "SINGLE" ]; then

#   samtools view -h $outdir/$sample.dedup.tmp.bam | grep -v chrM \
#   | samtools view -q 30 -b -h \
#   | samtools sort -o $outdir/$sample.filtered.bam

# else 
#   echo "Error: Unknown library layout '$layout' for sample '$sample'"
#   exit 1
# fi


# rm $outdir/$sample.dedup.tmp.bam



echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"