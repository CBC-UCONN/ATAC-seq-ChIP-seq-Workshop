#!/bin/bash
#SBATCH --job-name=02_nf-core_hic
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%j.out

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"


module load nextflow/25.04.6

nextflow run nf-core/hic \
   -profile xanadu \
   -work-dir /scratch/$USER/nf-core_hic_work \
   -resume \
   --input 02_samplesheet.csv \
   --outdir ../results/02_nf-hic_out \
   --fasta /core/cbc/tutorials/workshopdirs/Chip-ATAC/resources/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa \
   --restriction_site ^GATC,G^ANTC \
   --ligation_site GATCGATC,GANTGATC,GANTANTC,GATCANTC


echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"