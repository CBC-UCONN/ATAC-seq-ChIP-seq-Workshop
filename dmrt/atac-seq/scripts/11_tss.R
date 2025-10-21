# library(BiocManager)
# BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments",
#            "BSgenome.mm39.UCSC.mm39", "TxDb.Mmusculus.UCSC.mm39.knownGene",
#            "phastCons100way.UCSC.mm39"))


library(ATACseqQC)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)


# input the bamFile from the ATACseqQC package 
bamfile <- system.file("extdata", args[1], 
                        package="ATACseqQC", mustWork=TRUE)
bamfile.labels <- gsub(".bam", "", basename(bamfile))
alignment <- readGAlignments(bamfile)


# Estimate library complexity
estimateLibComplexity(readsDupFreq(bamfile))


# Compute fragment size distribution
pdf(paste0(bamfile.labels, "_fragSizeDist.pdf"))
fragSize <- fragSizeDist(bamfile, bamfile.labels)
dev.off()


# Transcription Start Site enrichment analysis
txs <- transcripts(TxDb.Mmusculus.UCSC.mm39.knownGene)
tsse <- TssEscore(alignment, step=10)
print(tsse)
plot(100*(-9:10-.5), tsse$values, type="l", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")





# library(ATACseqQC)
# library(TxDb.Mmusculus.UCSC.mm39.refGene)  # For GRCm39/mm39

# # Load your BAM file
# bamFile <- "your_mouse_sample.bam"

# # Get TSS positions from mm39
# txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
# tss <- promoters(txdb, upstream=0, downstream=1)

# # Calculate TSS enrichment
# tssEnrich <- TSSEscore(bamFile, tss)

# library(ATACseqQC)
# library(TxDb.Mmusculus.UCSC.mm39.refGene)  # For GRCm39/mm39

# # Load your BAM file
# bamFile <- "your_mouse_sample.bam"

# # Get TSS positions from mm39
# txdb <- TxDb.Mmusculus.UCSC.mm39.refGene
# tss <- promoters(txdb, upstream=0, downstream=1)

# # Calculate TSS enrichment
# tssEnrich <- TSSEscore(bamFile, tss)