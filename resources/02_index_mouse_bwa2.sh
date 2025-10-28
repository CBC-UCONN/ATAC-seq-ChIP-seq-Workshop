#!/bin/bash
#SBATCH --job-name=02_index_mouse_bwa2
#SBATCH --cpus-per-task=1
#SBATCH --mem=80G
#SBATCH --partition=mcbstudent
#SBATCH --qos=general
#SBATCH --output=%x_%j.out


echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

# Load required modules
module load bwa-mem2/2.2.1


pref=Mus_musculus.GRCm39.dna_rm.primary_assembly
mkdir ${pref}_bwa2 
bwa-mem2 index -p ${pref}_bwa2/${pref}_bwa2 $pref.fa


echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"