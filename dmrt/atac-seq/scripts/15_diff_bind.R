
library(DiffBind)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
library(ChIPseeker)
library(GenomeInfoDb)
library(rtracklayer)

# # Setup
out_dir <- "../results/15_diff_bind/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
meta_df <- read.csv("..//meta/atac-sra-meta.csv")

# Create vectors of file paths for BAM and peak files
bam_files <- file.path("../results/08_filter", paste0(meta_df$Run, ".filtered.sorted.bam"))
peak_files <- file.path("../results/12_macs", paste0(meta_df$Run, "_peaks.broadPeak"))

# Verify that all files exist, stop execution if any are missing
stopifnot(all(file.exists(bam_files)))
stopifnot(all(file.exists(peak_files)))

# Create sample sheet data frame
sample_sheet <- data.frame(
  SampleID = meta_df$Run,
  Condition = meta_df$source_name,
  bamReads = bam_files,
  Peaks = peak_files,
  PeakCaller = "macs",
  PeakFormat = "bed")

# Create DBA object
db <- dba(sampleSheet = sample_sheet)
cat("\nDBA object:\n")
db

# Count reads in peaks
# Very time intensive step!
cat("\nCounting reads in peaks...\n")
db <- dba.count(db, 
  summits = 250,
  minOverlap = 2,
  bUseSummarizeOverlaps = TRUE)


# Can use the following to save and reload DBA object if needed to avoid re-running time intensive 
# counting function 
# Save DBA object
cat("\nSaving DBA object...\n")
dba.save(db, "count", dir = out_dir)

# Reload DBA object (if needed)
db <- dba.load("count", dir = out_dir)


cat("\nCounts:\n")
dba.show(db)

# Define contrasts
db <- dba.contrast(db, categories = DBA_CONDITION, minMembers = 2)

# View contrasts
cat("\nShow contrasts:\n")
dba.show(db, bContrasts = TRUE)

# Perform differential analysis
cat("\nPerforming differential analysis...\n")
db <- dba.analyze(db)

# View results summary
dba.show(db, bContrasts = TRUE)

# Get results with FDR < 0.05
diff_peaks <- dba.report(db, th = 0.05)

# Number of differential peaks
cat("\nNumber of differential peaks (FDR < 0.05): \n")
length(diff_peaks)
head(diff_peaks)

# Generate diagnostic plots
pdf(file.path(out_dir, "diagnostic_plots.pdf"))
dba.plotHeatmap(db)
dba.plotPCA(db, DBA_CONDITION, label=DBA_ID)
dba.plotHeatmap(db, correlations = FALSE)
dba.plotVolcano(db, contrast = 1)
dba.plotVolcano(db, contrast = 2)
dba.plotVolcano(db, contrast = 3)
dba.plotMA(db)
dev.off()

# Get results for each contrast
res1 <- dba.report(db, contrast = 1, th=1)
res2 <- dba.report(db, contrast = 2, th=1)
res3 <- dba.report(db, contrast = 3, th=1)

cat("\nAnnotating regions...\n")
# Ensure that chromosome naming style matches taxonomic database
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
seqlevelsStyle(res1) <- seqlevelsStyle(txdb)[1]
seqlevelsStyle(res2) <- seqlevelsStyle(txdb)[1]
seqlevelsStyle(res3) <- seqlevelsStyle(txdb)[1]

# Annotate
anno_res1 <- annotatePeak(res1, TxDb = txdb, annoDb = "org.Mm.eg.db", tssRegion = c(-3000, 3000))
anno_res2 <- annotatePeak(res2, TxDb = txdb, annoDb = "org.Mm.eg.db", tssRegion = c(-3000, 3000))
anno_res3 <- annotatePeak(res3, TxDb = txdb, annoDb = "org.Mm.eg.db", tssRegion = c(-3000, 3000))

# Output annotated results to CSV
cat("\nWriting annotated results to CSV...\n")
write.csv(as.data.frame(anno_res1), file.path(out_dir, "contrast1_all.csv"), row.names=FALSE)
write.csv(as.data.frame(anno_res2), file.path(out_dir, "contrast2_all.csv"), row.names=FALSE)
write.csv(as.data.frame(anno_res3), file.path(out_dir, "contrast3_all.csv"), row.names=FALSE)

sig_res1 <- res1[res1$FDR <= 0.05, ]
sig_res2 <- res2[res2$FDR <= 0.05, ]
sig_res3 <- res3[res3$FDR <= 0.05, ]

# Output significant peaks to a BED file
export.bed(sig_res1, file.path(out_dir, "contrast1_significant_peaks.bed"))
export.bed(sig_res2, file.path(out_dir, "contrast2_significant_peaks.bed"))
export.bed(sig_res3, file.path(out_dir, "contrast3_significant_peaks.bed"))


