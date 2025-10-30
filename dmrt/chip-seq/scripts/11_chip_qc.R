library(ChIPQC)

args <- commandArgs(trailingOnly=TRUE)
filtered_bam <- args[1]
peaks <- args[2]
outdir <- args[3]

# Extract sample name from BAM file path
name <- gsub(".sorted.bam", "", basename(args[1]))

sample <- ChIPQCsample(filtered_bam)

ChIPQCreport(sample, peaks, reportName=paste0(name, "_ChIP_QC_report"), reportFolder=outdir)