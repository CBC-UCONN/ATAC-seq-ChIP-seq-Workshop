
## Load libraries
library(ChIPQC)

# Setup
out_dir <- "../results/11_chip_qc/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_df <- read.csv("../meta/chip-sra-meta.csv")
meta_df <- meta_df[1:4,]
meta_df
bam_files <- file.path("../results/08_filter", paste0(meta_df$Run, ".filtered.sorted.bam"))
peak_files <- file.path("../results/10_macs", paste0(meta_df$Run, "_peaks.narrowPeak"))

stopifnot(all(file.exists(bam_files)))
stopifnot(all(file.exists(peak_files)))

sample_sheet <- data.frame(
  SampleID = meta_df$Run,
  Condition = meta_df$source_name,
  bamReads = bam_files,
  Peaks = peak_files,
  PeakCaller = "macs",
  PeakFormat = "narrow")

## Create ChIPQC object
chipObj <- ChIPQC(sample_sheet, chromosomes=NULL) 

## Create ChIPQC report
ChIPQCreport(chipObj, reportName="ChIP_QC_report", reportFolder=outdir)

# View summary statistics
chipqc_obj

# Access specific QC metrics
QCmetrics(chipqc_obj)

# Plot quality metrics
plotCC(chipqc_obj)  # Cross-coverage plot
plotRegi(chipqc_obj)  # Reads in genomic regions
plotFrip(chipqc_obj)  # Fraction of reads in peaks

# Export QC metrics to CSV
qc_metrics <- QCmetrics(chipqc_obj)
write.csv(qc_metrics, file.path(out_dir, "qc_metrics.csv"), row.names = FALSE)
