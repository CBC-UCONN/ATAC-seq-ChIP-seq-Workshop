#!/bin/bash

# Check if argument is provided
if [ $# -eq 0 ]; then
    echo "Error: No argument provided"
    echo "Usage: $0 <directory_name>"
    exit 1
fi

# Create the results symlink
ln -s /core/cbc/tutorials/workshopdirs/Chip-ATAC/dmrt/atac-seq/results/$1 ../results/$1

# Check if results symlink was created successfully
if [ $? -eq 0 ]; then
    echo "Symlink created successfully: ../results/$1 -> /core/cbc/tutorials/workshopdirs/Chip-ATAC/dmrt/atac-seq/results/$1"
else
    echo "Error: Failed to create results symlink"
    exit 1
fi

# Create symlinks for all matching log files
source_logs="/core/cbc/tutorials/workshopdirs/Chip-ATAC/dmrt/atac-seq/scripts/logs/$1"*
log_count=0

for logfile in $source_logs; do
    if [ -e "$logfile" ]; then
        basename_log=$(basename "$logfile")
        ln -s "$logfile" logs/"$basename_log"
        if [ $? -eq 0 ]; then
            echo "Log symlink created: logs/$basename_log -> $logfile"
            ((log_count++))
        else
            echo "Warning: Failed to create symlink for $logfile"
        fi
    fi
done

if [ $log_count -eq 0 ]; then
    echo "Warning: No log files found matching $source_logs"
fi

echo "Done! Created 1 results symlink and $log_count log symlinks."