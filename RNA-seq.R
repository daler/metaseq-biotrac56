library(DESeq2)

library(Rsubread)
df <- featureCounts(
    c(
        "data/K562_1.chr11.bam",
        "data/K562_2.chr11.bam",
        "data/H1-hESC_1.chr11.bam",
        "data/H1-hESC_2.chr11.bam"),
    annot.ext="data/gencode.v19.annotation.chr11.gtf",
    isGTFAnnotationFile=TRUE)

write.table(df$counts, file='data/featureCounts.counts.txt', sep='\t')
write.table(df$stat, file='data/featureCounts.stats.txt', sep='\t')
count.data <- read.table('data/featureCounts.counts.txt', sep='\t')

colData <- data.frame(
    celltype=c('k562', 'k562', 'h1-hesc', 'h1-hesc')
)

dds <- DESeqDataSetFromMatrix(
    countData=count.data,
    colData=colData,
    design=~celltype)
colnames(dds) <- colnames(count.data)


dds$celltype <- relevel(dds$celltype, 'h1-hesc')
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
print(summary(res))

write.table(res, file='data/DESeq-results.txt', sep='\t')
