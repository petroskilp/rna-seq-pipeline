log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("KernSmooth")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

contrast <- c("condition", "Huntington", "TD")
res <- results(dds, contrast = contrast)

alpha <- 0.01 # Threshold on the p-value

# par(mfrow=c(1,2))

# Compute significance, with a maximum of 320 for the p-values set to 0 due to limitation of computation precision
res$sig <- -log10(res$padj)
sum(is.infinite(res$sig))

res[is.infinite(res$sig),"sig"] <- 350
# View(resultDESeq2[is.na(resultDESeq2$pvalue),])

# Select genes with a defined p-value (DESeq2 assigns NA to some genes)
genes.to.plot <- !is.na(res$pvalue)
# sum(genes.to.plot)
range(res[genes.to.plot, "log2FoldChange"])

## Volcano plot of adjusted p-values
cols <- densCols(res$log2FoldChange, res$sig)
cols[res$pvalue ==0] <- "purple"
res$pch <- 19
res$pch[res$pvalue ==0] <- 6

svg(snakemake@output[[1]])

plot(res$log2FoldChange, 
     res$sig, 
     col=cols, panel.first=grid(),
     main="Volcano plot", 
     xlab="Effect size: log2(fold-change)",
     ylab="-log10(adjusted p-value)",
     pch=res$pch, cex=0.4)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

## Plot the names of a reasonable number of genes, by selecting those begin not only significant but also having a strong effect size
gn.selected <- abs(res$log2FoldChange) > 2 & res$padj < alpha 
if(length(res$log2FoldChange[gn.selected])>0){
     text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.6)
}


dev.off()
