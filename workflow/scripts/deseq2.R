log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

contrast <- c("condition", snakemake@params[["contrast"]])
norm_counts = counts(dds, normalized=T)
res <- results(dds, contrast=contrast, parallel=parallel)

coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample_name", check.names=FALSE)
coldata <- coldata[order(row.names(coldata)), , drop=F]
condition <- snakemake@params[["model"]]
condition <- sub("~", "", condition)
for (i in unique(coldata[[condition]]))
{
    mean <- data.frame(mean=rowMeans(norm_counts[,c(rownames(coldata[coldata[[condition]]==i,]))]))
    mean <- cbind(rownames(mean), mean)
    rownames(mean) <- NULL
    colnames(mean) <- c("row", paste("baseMean_", i, sep=""))
    res<-merge(res, mean)
}
# shrink fold changes for lowly expressed genes
# use ashr so we can use `contrast` as conversion to coef is not trivial
# see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#extended-section-on-shrinkage-estimators
res <- lfcShrink(dds, contrast=contrast, res=res, type="ashr")

# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(data.frame("gene"=rownames(res),res), file=snakemake@output[["table"]], row.names=FALSE, sep='\t')
