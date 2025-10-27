#!/bin/bash
#SBATCH --job-name=02_index_mouse_bwa
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=%x_%j.out


echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load bwa/0.7.17


pref=Mus_musculus.GRCm39.dna_rm.primary_assembly
mkdir ${pref}_bwa 
bwa index -p ${pref}_bwa/${pref}_bwa $pref.fa


echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"