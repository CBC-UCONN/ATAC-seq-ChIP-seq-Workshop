library(ATACseqQC)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)

# Get command line arguments
args <- commandArgs(trailingOnly=TRUE)
unfiltered_bam <- args[1]
filtered_bam <- args[2]
outdir <- args[3]

# Extract sample name from BAM file path
name <- gsub(".sorted.bam", "", basename(args[1]))

# # Compute library complexity
# tryCatch({
# pdf(paste0(outdir, "/", name, "_libComplexity.pdf"))
# estimateLibComplexity(readsDupFreq(unfiltered_bam))
# dev.off()
# }, error=function(e) {
#   message("Error in library complexity calculation: ", e)
# })

# # Compute fragment size distribution
# tryCatch({
# pdf(paste0(outdir, "/", name, "_fragSizeDist.pdf"))
# fragSize <- fragSizeDist(filtered_bam, name)
# dev.off()
# }, error=function(e) {
#   message("Error in fragment size distribution calculation: ", e)
# })

# Transcription Start Site enrichment analysis
tryCatch({
aln <- readBamFile(filtered_bam, asMates=TRUE, bigFile=FALSE)
txs <- transcripts(TxDb.Mmusculus.UCSC.mm39.knownGene)
tsse <- TSSEscore(aln, txs, step=10)
print(tsse)
pdf(paste0(outdir, "/", name, "_TSSenrichment.pdf"))
plot(100*(-9:10-.5), tsse$values, type="l", 
     xlab="distance to TSS",
     ylab="aggregate TSS score")
dev.off()
}, error=function(e) {
  message("Error in TSS enrichment analysis: ", e)
})