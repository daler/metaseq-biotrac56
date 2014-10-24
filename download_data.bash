#!/bin/bash

cat << EOF
  This script downloads and prepares the following data:

    - BAMs for RNA-seq replicates in two different cell types (K562 and H1-hESC)

    - narrowPeak format BED-like files for some erythroid-specific
    transcription factors in K562 cells

    - bigWig files for the same TFs as well as for some histone modifications

    - annotations from the most recent GENCODE release for hg19.


  To speed up all analysis, only data from chr11 is downloaded for each data
  type.


  [NOTE]: When downloading BAM files, the warning message:

    [knet_seek] SEEK_END is not supported for HTTP. Offset is unchanged.

  is expected, harmless, and can be ignored.

EOF

set -e

# Pipe-delimited list of URL|destination.
#
# Add more data files here
cat > filelist << EOF
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqH1hescR1x75dAlignsRep1V2.bam|data/H1-hESC_1.chr11.bam
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqH1hescR1x75dAlignsRep2V2.bam|data/H1-hESC_2.chr11.bam
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bam|data/K562_1.chr11.bam
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCaltechRnaSeq/wgEncodeCaltechRnaSeqK562R1x75dAlignsRep2V2.bam|data/K562_2.chr11.bam
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsK562Gata1bIggmusPk.narrowPeak.gz|data/K562_GATA1.narrowPeak
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsK562Tal1sc12984IggmusPk.narrowPeak.gz|data/K562_TAL1.narrowPeak
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsK562Gata1bIggmusSig.bigWig|data/K562_GATA1.chr11.bigWig
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhTfbs/wgEncodeSydhTfbsK562Tal1sc12984IggmusSig.bigWig|data/K562_TAL1.chr11.bigWig
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHistone/wgEncodeBroadHistoneK562H3k4me3StdSig.bigWig|data/K562_H3K4me3.chr11.bigWig
http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeSydhHistone/wgEncodeSydhHistoneK562H3k27me3bUcdSig.bigWig|data/K562_H3K27me3.chr11.bigWig
ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz|data/gencode.v19.annotation.chr11.gtf
EOF


mkdir -p data

log () {
    echo "[ $(date) ] $1"
}

# Download just chr11; we have to index later.
for i in $(grep ".bam" filelist); do
    url=$(echo $i | cut -f1 -d "|")
    dest=$(echo $i | cut -f2 -d "|")
    if [ -e $dest ]; then
        log "$dest exists, skipping"
        continue
    fi
    log "$url > $dest"
    samtools view -b $url chr11 > $dest

    # samtools downloads a local copy of the index.  We don't need it anymore,
    # so clean that up here.
    rm $(basename $url).bai

    log "indexing $dest"
    samtools index $dest
done


# K562 peaks (TFs; expect them to be at K562-specific genes).  Only save those
# on chr11.
for i in $(grep ".narrowPeak" filelist); do
    url=$(echo $i | cut -f1 -d "|")
    dest=$(echo $i | cut -f2 -d "|")
    if [ -e $dest ]; then
        log "$dest exists, skipping"
        continue
    fi
    log "$url > $dest"
    curl $url | zcat | grep "^chr11" > $dest
done


# K562 bigWigs (for profiling). Use bigWigToBedGraph to only extract chr11, and
# then re-create subsetted bigWigs
bigWigSubset () {
    if [ ! -e hg19.genome ]; then
        fetchChromSizes hg19 > hg19.genome
    fi
    url=$1
    tmp=$(basename $url).tmp
    bigWigToBedGraph $url $tmp -chrom=chr11
    bedGraphToBigWig $tmp hg19.genome $2
    rm $tmp
}

for i in $(grep ".bigWig" filelist); do
    url=$(echo $i | cut -f1 -d "|")
    dest=$(echo $i | cut -f2 -d "|")
    if [ -e $dest ]; then
        log "$dest exists, skipping"
        continue
    fi
    log "$url > $dest"
    bigWigSubset $url $dest
done

# Download the chr11 features for gencode v19 annotations.
for i in $(grep "gencode.v19" filelist); do
    url=$(echo $i | cut -f1 -d "|")
    dest=$(echo $i | cut -f2 -d "|")
    if [ -e $dest ]; then
        log "$dest exists, skipping"
        continue
    fi
    log "$url > $dest"
    curl $url | zcat | grep "^chr11" > $dest
done

# Create a gffutils db
if [ ! -e data/gencode.v19.annotation.chr11.db ]; then
    log "building gffutils database from annotations"
    cat > gffutils-clean-gencode.py << EOF
import gffutils
gffutils.create_db(
    "data/gencode.v19.annotation.chr11.gtf",
    "data/gencode.v19.annotation.chr11.db",
    keep_order=False,
    sort_attribute_values=False,
    id_spec={'gene': 'gene_id', 'transcript': 'transcript_id'},
    verbose=True, merge_strategy='merge', infer_gene_extent=False,)
EOF
    python gffutils-clean-gencode.py
fi

log "Complete. Please see the 'data' directory for files"
