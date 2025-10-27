#!/bin/bash
#SBATCH --job-name=17_footprints
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --array=1-10

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

source ~/.bashrc
mamba activate tobias

# Store some paths as variables
meta_data=../meta/atac-sra-meta.csv
bamdir=../results/08_filter
peakdir=../results/12_macs

# Get Nth SRA run accession based on SLURM_ARRAY_TASK_ID from metadata csv file
sample=$(awk -F, -v row=${SLURM_ARRAY_TASK_ID} \
    'NR==1{for(i=1;i<=NF;i++)if($i=="Run")col=i}NR==row+1&&col{print $col}' $meta_data)

# Create output directory if it doesn't exist 
outdir=../results/17_footrprints/$sample
mkdir -p $outdir

genome=../../../resources/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa
peaks=$outdir/${sample}_peaks.filt.broadPeak
footprints=$outdir/$sample.footprints.bw
motifs_db=../../../resources/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt
tgt_motifs=$outdir/target_motifs.txt

# Get motifs of interest
grep -E -A 4 \
  --no-group-separator \
  "Foxl2|ESR2|Dmrt1|SOX9" \
  $motifs_db \
  > $tgt_motifs

# Filter peaks to remove mitochondrial peaks
grep -v "^MT" \
  $peakdir/${sample}_peaks.broadPeak \
  > $peaks 

TOBIAS ATACorrect \
  --bam $bamdir/$sample.filtered.sorted.bam \
  --genome $genome \
  --peaks $peaks \
  --outdir $outdir \
  --read_shift 0 0

TOBIAS ScoreBigwig \
  --signal $outdir/$sample.filtered.sorted_corrected.bw \
  --regions $peaks \
  --output $outdir/$sample.footprints.bw

TOBIAS BINDetect \
  --signals $outdir/$sample.footprints.bw \
  --motifs $tgt_motifs \
  --genome $genome \
  --peaks $peaks \
  --outdir $outdir

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"




