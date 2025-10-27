#!/bin/bash
#SBATCH --job-name=02_index_mouse_bwa
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=%x_%j.out


echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load seqtk/

seqkit fx2tab -l Mus_musculus.GRCm39.dna_rm.primary_assembly.fa.gz | awk '{print "chr" $1 "\t" $NF}' | \
    sort -k2 -nr > Mus_musculus.GRCm39.dna_rm.primary_assembly_sizes.txt

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"