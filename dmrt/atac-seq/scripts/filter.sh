# Filter non mtDNA reads from BAM file

samtools idxstats input.bam | cut -f 1 | grep -v MT | xargs samtools view -b input.bam > output.bam