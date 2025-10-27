#!/bin/bash
#SBATCH --job-name=01_download
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=%x_%j.out


echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"


# Download Mouse Genome Reference (GRCm39) primary assembly with masked repeats
wget https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_rm.primary_assembly.fa.gz

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"