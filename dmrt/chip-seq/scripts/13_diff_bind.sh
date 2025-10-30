#!/bin/bash
#SBATCH --job-name=13_diff_bind
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --partition=general
#SBATCH --qos=general
#SBATCH --output=logs/%x_%j.out

echo "Job running on: $(hostname)"
start=$(date +%s)
echo "Start time: $(date)"

source ~/.bashrc
mamba activate r_atac

Rscript 13_diff_bind.R

echo "End time: $(date)"
echo "Elapsed time: $(date -ud "@$(($(date +%s)-start))" +'%H hr %M min %S sec')"