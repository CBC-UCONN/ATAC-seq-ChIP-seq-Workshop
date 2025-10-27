#!/bin/bash
#SBATCH --job-name=13_frip
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-10

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

module load bedtools/2.29.0
module load samtools/1.9

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
bamdir=../results/08_filter
peakdir=../results/12_macs
outdir=../results/13_frip

# Create output directory if it doesn't exist 
mkdir -p $outdir

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Count reads in peaks
reads_in_peaks=$(bedtools intersect \
  -a $peakdir/${sample}_peaks.broadPeak \
  -b $bamdir/$sample.filtered.bam \
  -bed \
  -c \
  | awk '{sum+=$NF} END {print sum}')

# Count total mapped reads
total_reads=$(samtools view -c $bamdir/$sample.filtered.bam)

# Calculate FRiP
frip=$(awk -v rip=$reads_in_peaks -v total=$total_reads \
       'BEGIN {printf "%.4f", rip/total}')

outfile=$outdir/${sample}_frip.txt
echo "Reads in peaks: $reads_in_peaks" > $outfile
echo "Total reads: $total_reads" >> $outfile
echo "FRiP score: $frip" >> $outfile
echo ""

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"


