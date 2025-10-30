# Load required libraries
library(csaw)
library(edgeR)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(org.Mm.eg.db)
library(rtracklayer)

# Setup
out_dir <- "../results/14_csaw/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

meta_df <- read.csv("../meta/chip-sra-meta.csv")
meta_df <- meta_df[meta_df$antibody == "Abcam ab4729"  & !is.na(meta_df$antibody), ]

bam_files <- file.path("../results/08_filter", paste0(meta_df$Run, ".filtered.sorted.bam"))
stopifnot(all(file.exists(bam_files)))

# Setup design matrix
groups <- factor(meta_df$genotype)
design <- model.matrix(~0 + groups)
colnames(design) <- c("eGfp", "CreERT2") # Rename columns to remove spaces

# Setup contrasts
contrast1 <- makeContrasts(eGfp_vs_CreERT2 = eGfp - CreERT2, levels=design)

# Compute average fragment length
cat("\nEstimating fragment length...\n")
x <- correlateReads(bam_files)
frag.len <- maximizeCcf(x)
frag.len

# Count reads in windows across the genome
cat("\nCounting reads in windows...\n")
win.data <- windowCounts(bam_files, width=10, ext=frag.len) 
win.data$totals
win.data

# Filter windows based on background read counts 
cat("\nFiltering windows...\n")
bins <- windowCounts(bam_files, bin=TRUE, width=10000)
filter.stat <- filterWindowsGlobal(win.data, bins)
min.fc <- 3 # Minimum fold-change above background
keep <- filter.stat$filter > log2(min.fc)
summary(keep)

# Plot to show impact of filtering
pdf(file.path(out_dir, "filtering_histogram.pdf"))
hist(filter.stat$filter, main="", breaks=50,
    xlab="Background abundance (log2-CPM)")
abline(v=log2(min.fc), col="red")
dev.off()

# Do the actual filtering
filtered.data <- win.data[keep,]

# Compute normalization factors
cat("\nComputing normalization factors...\n")
filtered.data <- normFactors(bins, se.out = filtered.data)
filtered.data$norm.factors

# Fit model
cat("\nPerforming differential binding analysis...\n")
y <- asDGEList(filtered.data)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# Test for differential accessibility
results1 <- glmQLFTest(fit, contrast=contrast1)


##############
# Would want to repeat this for each contrast
# Examine results
cat("\nPreview of results:\n")
head(results1$table)

merged1 <- mergeResults(filtered.data, results1$table, tol=100, merge.args=list(max.width=5000))
merged1 # Has regions, combined p-values, and the best window among merged windows

is.sig <- merged1$combined$FDR <= 0.05
cat("\nDifferentially open regions (GC vs SC):\n")
table(merged1$combined$direction[is.sig])


# Annotate merged regions
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
seqlevelsStyle(merged1$regions) <- seqlevelsStyle(txdb)[1]
anno <- detailRanges(merged1$regions, txdb=txdb,
  orgdb=org.Mm.eg.db, promoter=c(3000, 3000), dist=5000)
merged1$regions$overlap <- anno$overlap
merged1$regions$left <- anno$left
merged1$regions$right <- anno$right
write.csv(as.data.frame(merged1), file.path(out_dir, "contrast1_results.csv"), row.names = FALSE)

# Output significant results to bed file
test <- merged1$regions[is.sig]
test$score <- -10*log10(merged1$combined$FDR[is.sig])
names(test) <- paste0("region", 1:sum(is.sig))
export(test, file.path(out_dir, "clusters.bed"))
