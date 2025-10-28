#!/bin/bash

name="${filename%.sh}"
results_dir="/core/cbc/tutorials/workshopdirs/Chip-ATAC/dmrt/atac-seq/results/"
logs_dir="/core/cbc/tutorials/workshopdirs/Chip-ATAC/dmrt/atac-seq/scripts/logs/"

# Check if argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No argument provided"
    echo "Usage: $0 <script>"
    exit 1
fi

mkdir -p ../results
mkdi -p logs

if [ ! -d "../results/$name" ]; then
  ln -s "$results_dir/$name" "../results/$name"

  source_logs="$logs_dir/$name"*
  log_count=0
  for logfile in $source_logs; do
    basename_log=$(basename "$logfile")
    ln -s "$logfile" logs/"$basename_log"
  done
  if [ $log_count -eq 0 ]; then
      echo "Warning: No log files found matching $source_logs"
  fi
else
  echo "../results/$name already exists"
  echo "If you wish to remove it, run: rm -r ../results/$name"
  exit 1
fi







