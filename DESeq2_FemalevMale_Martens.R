library(DESeq2)
library(pheatmap)
library(dplyr)
library(Seurat)
library(scuttle)

setwd("C:/Users/marte/Desktop/School/University of Michigan/Winter 2022/BIOINF545/Project")

###################################################################################################
## Immune Original Seurat Object
Imm <- readRDS("train.Imm.seur.rds")
Imm <- metadata(Imm)

## Subset into comparison groups
colitis <- subset(Imm, subset = Disease %in% c("Colitis"))
colitis_inflamed <- subset(colitis, subset = Health %in% c("Inflamed"))
colitis_noninflamed <- subset(colitis, subset = Health %in% c("Non-inflamed"))
HC <- subset(Imm, subset = Disease %in% c("HC"))
## Aggregate each group
col_inf_summed <- aggregate(colitis_inflamed)
col_noninf_summed <- aggregate(colitis_noninflamed)
HC_summed <- aggregate(HC)
## DEseq2 by gender
res_col_inf <- resDEseq2(col_inf_summed)
res_col_noninf <- resDEseq2(col_noninf_summed)
res_HC <- resDEseq2(HC_summed)
Imm_inflamed <- res_col_inf[!(rownames(res_col_inf) %in% c(rownames(res_col_noninf), rownames(res_HC))),]


summary(res_col_inf)
summary(res_col_noninf)
summary(res_HC)

write.csv(res_col_inf, "Gender Differences/Immune/colitis_inf_femalevmale_IMM.csv")
write.csv(res_col_noninf, "Gender Differences/Immune/colitis_noninf_femalevmale_IMM.csv")
write.csv(res_HC, "Gender Differences/Immune/HC_femalevmale_IMM.csv")
write.csv(Imm_inflamed, "Gender Differences/Immune/HC_femalevmale_Imm_inflamed_unique.csv")

###################################################################################################
## Epithelial Original Seurat Object
Epi <- readRDS("train.Epi.seur.rds")
Epi <- metadata(Epi)

## Subset into comparison groups
colitis <- subset(Epi, subset = Disease %in% c("Colitis"))
colitis_inflamed <- subset(colitis, subset = Health %in% c("Inflamed"))
colitis_noninflamed <- subset(colitis, subset = Health %in% c("Non-inflamed"))
HC <- subset(Epi, subset = Disease %in% c("HC"))
## Aggregate each group
col_inf_summed <- aggregate(colitis_inflamed)
col_noninf_summed <- aggregate(colitis_noninflamed)
HC_summed <- aggregate(HC)
## DEseq2 by gender
res_col_inf <- resDEseq2(col_inf_summed)
res_col_noninf <- resDEseq2(col_noninf_summed)
res_HC <- resDEseq2(HC_summed)
Epi_inflamed <- res_col_inf[!(rownames(res_col_inf) %in% c(rownames(res_HC))),]


range <- max(abs(data.matrix(Epi_inflamed$log2FoldChange)))
pheatmap(data.matrix(Epi_inflamed$log2FoldChange), breaks = seq(-range, range, length.out = 100), cluster_rows = T, cluster_cols = F)

summary(res_col_inf)
summary(res_col_noninf)
summary(res_HC)

write.csv(res_col_inf, "Gender Differences/Epithelial/colitis_inf_femalevmale_Epi.csv")
write.csv(res_col_noninf, "Gender Differences/Epithelial/colitis_noninf_femalevmale_Epi.csv")
write.csv(res_HC, "Gender Differences/Epithelial/HC_femalevmale_Epi.csv")
write.csv(Epi_inflamed, "Gender Differences/Epithelial/HC_femalevmale_Epi_inflamed_unique.csv")

###################################################################################################
## Fibroblasts Original Seurat Object
Fib <- readRDS("train.Fib.seur.rds")
Fib <- metadata(Fib)

## Subset into comparison groups
colitis <- subset(Fib, subset = Disease %in% c("Colitis"))
colitis_inflamed <- subset(colitis, subset = Health %in% c("Inflamed"))
colitis_noninflamed <- subset(colitis, subset = Health %in% c("Non-inflamed"))
HC <- subset(Fib, subset = Disease %in% c("HC"))
## Aggregate each group
col_inf_summed <- aggregate(colitis_inflamed)
col_noninf_summed <- aggregate(colitis_noninflamed)
HC_summed <- aggregate(HC)
## DEseq2 by gender
res_col_inf <- resDEseq2(col_inf_summed)
res_col_noninf <- resDEseq2(col_noninf_summed)
res_HC <- resDEseq2(HC_summed)

summary(res_col_inf)
summary(res_col_noninf)
summary(res_HC)

Fib_inflamed <- res_col_inf[!(rownames(res_col_inf) %in% c(rownames(res_col_noninf), rownames(res_HC))),]

write.csv(res_col_inf, "Gender Differences/Fibroblasts/colitis_inf_femalevmale_Fib.csv")
write.csv(res_col_noninf, "Gender Differences/Fibroblasts/colitis_noninf_femalevmale_Fib.csv")
write.csv(res_HC, "Gender Differences/Fibroblasts/HC_femalevmale_Fib.csv")
write.csv(Fib_inflamed, "Gender Differences/Fibroblasts/HC_femalevmale_Fib_inflamed_unique.csv")


#####################################################################################################
#####################################################################################################
## Functions
metadata <- function(x){
  meta <- read.table("all.meta2.txt", header = T, sep = '\t')
  gender <- read.csv("gen2.csv")
  merged <- left_join(meta, gender, by = c("Subject" = "Subject.ID"))
  row.names(merged) <- merged$NAME
  x <- AddMetaData(
    object = x,
    metadata = merged
  )
  x
}

aggregate <- function(x){
  agg <- c("Subject", "Cluster")
  x.sce <- as.SingleCellExperiment(x)
  summed <- aggregateAcrossCells(x.sce, id=colData(x.sce)[,agg])
  colnames(summed) <- apply(colData(summed)[,agg], 1, function(x) paste(x, collapse="_"))
  summed
}

resDEseq2 <- function(x){
  levels(x$Gender) <- c("Female", "Male")
  dds <- DESeqDataSet(x, design = ~Gender)
  dds <- dds[rowSums(counts(dds))>=dim(dds)[2],]
  dds <- DESeq(dds)
  female_v_male <- results(dds,contrast=c("Gender", "Female", "Male"),pAdjustMethod="fdr")
  sig <- female_v_male[which(female_v_male$padj < 0.05),]
  sig
}
