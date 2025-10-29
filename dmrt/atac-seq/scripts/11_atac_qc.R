library(ATACseqQC)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
unfiltered_bam <- args[1]
filtered_bam <- args[2]
outdir <- args[3]

# Extract sample name from BAM file path
name <- gsub(".sorted.bam", "", basename(args[1]))

# Compute library complexity
pdf(paste0(outdir, "/", name, "_libComplexity.pdf"))
# Use unfiltered BAM for library complexity estimation
estimateLibComplexity(readsDupFreq(unfiltered_bam))
dev.off()


# Compute fragment size distribution
pdf(paste0(outdir, "/", name, "_fragSizeDist.pdf"))
fragSize <- fragSizeDist(filtered_bam, name)
dev.off()

