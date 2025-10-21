

outdir=../data/raw-fastq

# Make output directory path
mkdir -p $outdir

for f in /core/cbc/tutorials/workshopdirs/Chip-ATAC/dmrt/hi-c/data/raw-fastq/*fastq.gz; do
    ln -s $f $outdir/$(basename $f)
done
