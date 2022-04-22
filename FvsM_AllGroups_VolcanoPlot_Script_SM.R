BiocManager::install("calibrate")
library("calibrate")
FvsM_Imm_NonInflamed_UniqueALL <- read.csv("femalevmale_Imm_noninflamed_uniqueALL.csv", header=TRUE)
with(FvsM_Imm_NonInflamed_UniqueALL, plot(log2FoldChange , -log10(pvalue), pch=20, main="Female vs. Male - Immune Group - Non-Inflamed", xlim=c(-5,5), ylim=c(0,15)))
with(subset(FvsM_Imm_NonInflamed_UniqueALL, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(FvsM_Imm_NonInflamed_UniqueALL, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log10(padj), labs=X, cex=.9))

FvsM_Imm_Inflamed_UniqueALL <- read.csv("femalevmale_Imm_inflamed_uniqueALL.csv", header=TRUE)
with(FvsM_Imm_Inflamed_UniqueALL, plot(log2FoldChange , -log10(pvalue), pch=20, main="Female vs. Male - Immune Group - Inflamed", xlim=c(-5,5), ylim=c(0,15)))
with(subset(FvsM_Imm_Inflamed_UniqueALL, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(FvsM_Imm_Inflamed_UniqueALL, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log10(padj), labs=X, cex=.9))



FvsM_Fib_NonInflamed_UniqueALL <- read.csv("femalevmale_Fib_noninflamed_uniqueALL.csv", header=TRUE)
with(FvsM_Fib_NonInflamed_UniqueALL, plot(log2FoldChange , -log10(pvalue), pch=20, main="Female vs. Male - Fibroblast Group - Non-Inflamed", xlim=c(-5,5), ylim=c(0,10)))
with(subset(FvsM_Fib_NonInflamed_UniqueALL, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(FvsM_Fib_NonInflamed_UniqueALL, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log10(padj), labs=X, cex=.9))

FvsM_Fib_Inflamed_UniqueALL <- read.csv("femalevmale_Fib_inflamed_uniqueALL.csv", header=TRUE)
with(FvsM_Fib_Inflamed_UniqueALL, plot(log2FoldChange , -log10(pvalue), pch=20, main="Female vs. Male - Fibroblast Group - Inflamed", xlim=c(-5,5), ylim=c(0,10)))
with(subset(FvsM_Fib_Inflamed_UniqueALL, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(FvsM_Fib_Inflamed_UniqueALL, padj<.05 & abs(log2FoldChange)>1.5), textxy(log2FoldChange, -log10(padj), labs=X, cex=.9))



FvsM_Epi_NonInflamed_UniqueALL <- read.csv("femalevmale_Epi_noninflamed_uniqueALL.csv", header=TRUE)
with(FvsM_Epi_NonInflamed_UniqueALL, plot(log2FoldChange , -log10(pvalue), pch=20, main="Female vs. Male - Epithelial Group - Non-Inflamed", xlim=c(-5,5), ylim=c(0,10)))
with(subset(FvsM_Epi_NonInflamed_UniqueALL, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(FvsM_Epi_NonInflamed_UniqueALL, padj<.05 & abs(log2FoldChange)>2.5), textxy(log2FoldChange, -log10(padj), labs=X, cex=.9))

FvsM_Epi_Inflamed_UniqueALL <- read.csv("femalevmale_Epi_inflamed_uniqueALL.csv", header=TRUE)
with(FvsM_Epi_Inflamed_UniqueALL, plot(log2FoldChange , -log10(pvalue), pch=20, main="Female vs. Male - Epithelial Group - Inflamed", xlim=c(-5,5), ylim=c(0,10)))
with(subset(FvsM_Epi_Inflamed_UniqueALL, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(FvsM_Epi_Inflamed_UniqueALL, padj<.05 & abs(log2FoldChange)>2.5), textxy(log2FoldChange, -log10(padj), labs=X, cex=.9))
