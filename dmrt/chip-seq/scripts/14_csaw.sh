#!/bin/bash
#SBATCH --job-name=14_csaw
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%j.out

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

source ~/.bashrc
# TODO: Change the name to match slides
mamba activate r_atac

Rscript 14_csaw.R

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"