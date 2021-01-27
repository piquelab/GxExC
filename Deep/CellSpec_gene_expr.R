# Anthony Findley
# 9/25/2020

# Purpose: Plot cpm expression values for genes knkown to be expressed specificaly in LCLs, IPSCs, or CMs to show that we have different cell types. (Supplemental figure 18)

R

library(DESeq2)
library(tidyverse)

load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/LCL_combined/P1_R2/out_data_DLCL_P1R2_04.19.19/data_objects/DESeq2_DLCL_P1R2_04.19.19.RData")

# Try the variance stablizing transformation (vst) first
LCL_vsd <- vst(ddsFull, blind=FALSE)

# To get the matrix of normalized values
LCL_normCounts <- assay(LCL_vsd)

# Calculate the variance for each transcript
RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}

LCL_MeanVar <- data.frame(Mean = rowMeans(LCL_normCounts), Variance = RowVar(LCL_normCounts))

# Add ensg
anno <- read.table("/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.tsv.gz", header=T, sep="\t", stringsAsFactors=F)
LCL_MeanVar <- merge(LCL_MeanVar, anno[,c("Transcript.stable.ID", "Gene.stable.ID", "Gene.name", "Gene.type")], by.x=0, by.y="Transcript.stable.ID")

LCL_MeanVar$CellType <- "LCL"

#################
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/IPSC_combined/P1_R2/out_data_DIPSC_P1R2_04.19.19/data_objects/DESeq2_DIPSC_P1R2_04.19.19.RData")
IPSC_vsd <- vst(ddsFull, blind=FALSE)
IPSC_normCounts <- assay(IPSC_vsd)
IPSC_MeanVar <- data.frame(Mean = rowMeans(IPSC_normCounts), Variance = RowVar(IPSC_normCounts))
IPSC_MeanVar <- merge(IPSC_MeanVar, anno[,c("Transcript.stable.ID", "Gene.stable.ID", "Gene.name", "Gene.type")], by.x=0, by.y="Transcript.stable.ID")
IPSC_MeanVar$CellType <- "IPSC"

load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/CM_combined/R2_R3/out_data_DCM_R2R3_04.19.19/data_objects/DESeq2_DCM_R2R3_04.19.19.RData")
CM_vsd <- vst(ddsFull, blind=FALSE)
CM_normCounts <- assay(CM_vsd)
CM_MeanVar <- data.frame(Mean = rowMeans(CM_normCounts), Variance = RowVar(CM_normCounts))
CM_MeanVar <- merge(CM_MeanVar, anno[,c("Transcript.stable.ID", "Gene.stable.ID", "Gene.name", "Gene.type")], by.x=0, by.y="Transcript.stable.ID")
CM_MeanVar$CellType <- "CM"

# Put these all together again
all_MeanVar <- rbind(LCL_MeanVar, IPSC_MeanVar, CM_MeanVar)
all_MeanVar$CellType <- factor(all_MeanVar$CellType, levels=c("LCL", "IPSC", "CM"))

# Pick out representative genes
all_MeanVar %>% filter(Gene.name %in% c("IGKC", "HLA-B", "CD74", "SOX2", "NANOG", "POU5F1", "TNNT2", "MYH7", "MYH6")) %>% group_by(CellType, Gene.name) %>% summarise(mean = mean(Mean)) -> repGenes_mean

###
# Actually, I think I want the expression values for each sample, so use the "normCounts" files, select just the genes we want, and average the transcripts for each gene in a library.

LCL_normCounts_gene <- merge(LCL_normCounts, anno[,c("Transcript.stable.ID", "Gene.stable.ID", "Gene.name")], by.x=0, by.y="Transcript.stable.ID")
LCL_normCounts_gene %>% filter(Gene.name %in% c("IGKC", "HLA-B", "CD74", "SOX2", "NANOG", "POU5F1", "TNNT2", "MYH7", "MYH6")) %>% pivot_longer(`DLCL1R1-HT10`:`DLCL2R2-HT96`, names_to = "library", values_to = "counts") %>% group_by(library, Gene.name) %>% summarise(gene_count = mean(counts)) -> LCL_repGenes_byLib
LCL_repGenes_byLib$CellType <- "LCL"

IPSC_normCounts_gene <- merge(IPSC_normCounts, anno[,c("Transcript.stable.ID", "Gene.stable.ID", "Gene.name")], by.x=0, by.y="Transcript.stable.ID")
IPSC_normCounts_gene %>% filter(Gene.name %in% c("IGKC", "HLA-B", "CD74", "SOX2", "NANOG", "POU5F1", "TNNT2", "MYH7", "MYH6")) %>% pivot_longer(`DIPSC1R1-HT10`:`DIPSC2R2-HT96`, names_to = "library", values_to = "counts") %>% group_by(library, Gene.name) %>% summarise(gene_count = mean(counts)) -> IPSC_repGenes_byLib
IPSC_repGenes_byLib$CellType <- "IPSC"

CM_normCounts_gene <- merge(CM_normCounts, anno[,c("Transcript.stable.ID", "Gene.stable.ID", "Gene.name")], by.x=0, by.y="Transcript.stable.ID")
CM_normCounts_gene %>% filter(Gene.name %in% c("IGKC", "HLA-B", "CD74", "SOX2", "NANOG", "POU5F1", "TNNT2", "MYH7", "MYH6")) %>% pivot_longer(`DCM1R1-HT10`:`DCM2R2-HT96`, names_to = "library", values_to = "counts") %>% group_by(library, Gene.name) %>% summarise(gene_count = mean(counts)) -> CM_repGenes_byLib
CM_repGenes_byLib$CellType <- "CM"

# Combine cells and plot
all_repGenes_byLib <- rbind(LCL_repGenes_byLib, IPSC_repGenes_byLib, CM_repGenes_byLib)

all_repGenes_byLib$Gene.name <- factor(all_repGenes_byLib$Gene.name, levels = c("IGKC", "HLA-B", "CD74", "SOX2", "NANOG", "POU5F1", "TNNT2", "MYH7", "MYH6"))
all_repGenes_byLib$CellType <- factor(all_repGenes_byLib$CellType, levels = c("LCL", "IPSC", "CM"))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/repGenes_all.box.pdf")
ggplot(all_repGenes_byLib, aes(x=Gene.name, y=gene_count, fill=CellType)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + theme_bw() + ylab("VST average counts") + scale_fill_manual(values=c("#C77CFF", "#619CFF", "#00BA38"))
dev.off()

# Try breaking up the genes into 3 plots based on the cell types they represent
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/repGenes_all_LCLgenes.box.pdf")
ggplot(all_repGenes_byLib %>% filter(Gene.name %in% c("IGKC", "HLA-B", "CD74")), aes(x=Gene.name, y=gene_count, fill=CellType)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + theme_bw() + ylab("VST average counts") + scale_fill_manual(values=c("#C77CFF", "#619CFF", "#00BA38"))
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/repGenes_all_IPSCgenes.box.pdf")
ggplot(all_repGenes_byLib %>% filter(Gene.name %in% c("SOX2", "NANOG", "POU5F1")), aes(x=Gene.name, y=gene_count, fill=CellType)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + theme_bw() + ylab("VST average counts") + scale_fill_manual(values=c("#C77CFF", "#619CFF", "#00BA38"))
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/repGenes_all_CMgenes.box.pdf")
ggplot(all_repGenes_byLib %>% filter(Gene.name %in% c("TNNT2", "MYH7", "MYH6")), aes(x=Gene.name, y=gene_count, fill=CellType)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + theme_bw() + ylab("VST average counts") + scale_fill_manual(values=c("#C77CFF", "#619CFF", "#00BA38"))
dev.off()
