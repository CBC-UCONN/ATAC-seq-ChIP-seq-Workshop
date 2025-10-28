#!/bin/bash

name="${1%.sh}"
results_dir="/core/cbc/tutorials/workshopdirs/Chip-ATAC/dmrt/chip-seq/results"
logs_dir="/core/cbc/tutorials/workshopdirs/Chip-ATAC/dmrt/chip-seq/scripts/logs"

# Check if argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No argument provided"
    echo "Usage: $0 <script>"
    exit 1
fi

mkdir -p ../results
mkdir -p logs

if [ -d "../results/$name"]; then
  if [ -L "../results/$name" ]; then
    echo "Removing existing symlink ../results/$name"
    rm "../results/$name"
  else
    echo "../results/$name already exists and is not a symlink"
    echo "If you wish to remove it, run: rm -rf ../results/$name"
    exit 1
  fi
fi

if [ ! -d "../results/$name" ]; then
  ln -s "$results_dir/$name" "../results/$name"
  if [ $? -eq 0 ]; then
    echo "Symlink created: ../results/$name"
  fi

  source_logs="$logs_dir/$name"*
  log_count=0
  for logfile in $source_logs; do
    basename_log=$(basename "$logfile")
    ln -s "$logfile" logs/"$basename_log"
  done
  if [ $log_count -eq 0 ]; then
    echo "No log files found matching"
  else
    echo "Symlinks to log files created: logs/$name*.out"
  fi
else
  echo "../results/$name already exists"
  echo "If you wish to remove it, run: rm -rf ../results/$name"
  exit 1
fi







