# Anthony Findley
# 1/3/20

# Purpose: Run GO, KEGG, and Reactome on DEGs from CMs. This compares each treatment to the background separately.
# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/CM_combined/R2_R3

# this script runs {GO, KEGG, Reactome} enrichment analyses on DEGs from DESeq2
library(tidyverse)
library(knitr)
#library(DESeq2)
# the following 2 packages work only in R/3.5.2 in new environment
library(annotables)
library(clusterProfiler)
#library(qqman)
#biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)
library(ReactomePA)
library(pathview)

sessionInfo()

# treatment="T13C1"

# First get list of treatments to run
treatments <- list.files("out_data_DCM_R2R3_04.19.19/stats/")
treatments <- sapply(strsplit(treatments, "_"), "[", 6)
treatments <- sapply(strsplit(treatments, "\\.."), "[", 1)
treatments <- treatments[grep("T", treatments)] # Gets rid of control vs control

# open background genes:
# background <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/CM_combined/R2_R3/out_data_DCM_R2R3_04.19.19/stats/DCM_R2R3_04.19.19_DEG_stats_CO1.txt", header=TRUE, stringsAsFactors=F)
# background <- unique(background[,"g.id"])
# background.df <- bitr(background, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
# background <- as.character(background.df[,2])

for (i in 1:length(treatments)){
        treatment <- treatments[i]

        # read in the genes:
        # gene_list is a vector of DE genes
        gene_list <- read.table(paste0("/nfs/rprdata/Anthony/GxE_iPSC/Deep/CM_combined/R2_R3/out_data_DCM_R2R3_04.19.19/stats/DCM_R2R3_04.19.19_DEG_stats_", treatment, ".txt"), stringsAsFactors=F, header=TRUE)

        # First get the background (all genes tested with p-value)
        background <- unique(gene_list[,"g.id"])
        background.df <- bitr(background, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
        background <- as.character(background.df[,2])

        # Now get upregulated genes
        up_DEG <- gene_list[which(gene_list$padj < 0.1 & gene_list$logFC > 0.25),]
        up_DEG <- unique(up_DEG[,"g.id"])
        up_DEGs.df <- bitr(up_DEG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
        up_DEGs_ENTREZ <- as.character(up_DEGs.df[,2])

        # GO over-representation test (takes a couple minutes):
        GO <- enrichGO(up_DEGs_ENTREZ, OrgDb="org.Hs.eg.db", keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.1, pAdjustMethod = "BH", universe=background,qvalueCutoff = 0.1, minGSSize = 10, maxGSSize = 500,readable = FALSE, pool = FALSE)

        # KEGG over-representation test:
        KEGG <- enrichKEGG(gene=up_DEGs_ENTREZ, organism='hsa', keyType = "kegg",pvalueCutoff = 0.1, pAdjustMethod = "BH", universe=background,qvalueCutoff = 0.1, minGSSize = 10)

        # REACTOME analysis:
        REACT <- enrichPathway(gene=up_DEGs_ENTREZ, organism="human",pvalueCutoff = 0.1, pAdjustMethod = "BH", universe=background,qvalueCutoff=0.1, minGSSize=10)

        #save the results:
        system("mkdir -p GO_Results")
        write.table(GO@result, file=paste0("./GO_Results/GO_", treatment, "_upDEGs.txt"), sep="\t", quote=FALSE, col.names=T)
        write.table(KEGG@result, file=paste0("./GO_Results/KEGG_", treatment, "_upDEGs.txt"), sep="\t", quote=FALSE, col.names=T)
        write.table(REACT@result, file=paste0("./GO_Results/REACTOME_", treatment, "_upDEGs.txt"), sep="\t", quote=FALSE, col.names=T)

        # Just the significant terms in one DF:
        newGO <- GO@result[GO@result$Count>1 & GO@result$qvalue<0.1,]
        newKEGG <- KEGG@result[KEGG@result$Count>1 & KEGG@result$qvalue<0.1,]
        newREACT <- REACT@result[REACT@result$Count>1 & REACT@result$qvalue<0.1,]

        merged <- rbind(newGO, newKEGG, newREACT)
        system("mkdir -p GO_Results/significant")
        write.table(merged, file=paste0("./GO_Results/significant/", treatment, "_enrichment_summary_upDEGs.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=T)


        ### Now do downregulated
        down_DEG <- gene_list[which(gene_list$padj < 0.1 & gene_list$logFC < -0.25),]
        down_DEG <- unique(down_DEG[,"g.id"])
        down_DEGs.df <- bitr(down_DEG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
        down_DEGs_ENTREZ <- as.character(down_DEGs.df[,2])
        GO <- enrichGO(down_DEGs_ENTREZ, OrgDb="org.Hs.eg.db", keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.1, pAdjustMethod = "BH", universe=background,qvalueCutoff = 0.1, minGSSize = 10, maxGSSize = 500,readable = FALSE, pool = FALSE)
        KEGG <- enrichKEGG(gene=down_DEGs_ENTREZ, organism='hsa', keyType = "kegg",pvalueCutoff = 0.1, pAdjustMethod = "BH", universe=background,qvalueCutoff = 0.1, minGSSize = 10)
        REACT <- enrichPathway(gene=down_DEGs_ENTREZ, organism="human",pvalueCutoff = 0.1, pAdjustMethod = "BH", universe=background,qvalueCutoff=0.1, minGSSize=10)
        system("mkdir -p GO_Results")
        write.table(GO@result, file=paste0("./GO_Results/GO_", treatment, "_downDEGs.txt"), sep="\t", quote=FALSE, col.names=T)
        write.table(KEGG@result, file=paste0("./GO_Results/KEGG_", treatment, "_downDEGs.txt"), sep="\t", quote=FALSE, col.names=T)
        write.table(REACT@result, file=paste0("./GO_Results/REACTOME_", treatment, "_downDEGs.txt"), sep="\t", quote=FALSE, col.names=T)
        newGO <- GO@result[GO@result$Count>1 & GO@result$qvalue<0.1,]
        newKEGG <- KEGG@result[KEGG@result$Count>1 & KEGG@result$qvalue<0.1,]
        newREACT <- REACT@result[REACT@result$Count>1 & REACT@result$qvalue<0.1,]
        merged <- rbind(newGO, newKEGG, newREACT)
        system("mkdir -p GO_Results/significant")
        write.table(merged, file=paste0("./GO_Results/significant/", treatment, "_enrichment_summary_downDEGs.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=T)


        ### Now do either up or downregulated
        all_DEG <- gene_list[which(gene_list$padj < 0.1 & abs(gene_list$logFC) > 0.25),] # all DE genes
        all_DEG <- unique(all_DEG[,"g.id"])
        all_DEGs.df <- bitr(all_DEG, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
        all_DEGs_ENTREZ <- as.character(all_DEGs.df[,2])
        GO <- enrichGO(all_DEGs_ENTREZ, OrgDb="org.Hs.eg.db", keyType = "ENTREZID", ont = "BP", pvalueCutoff = 0.1, pAdjustMethod = "BH", universe=background,qvalueCutoff = 0.1, minGSSize = 10, maxGSSize = 500,readable = FALSE, pool = FALSE)
        KEGG <- enrichKEGG(gene=all_DEGs_ENTREZ, organism='hsa', keyType = "kegg",pvalueCutoff = 0.1, pAdjustMethod = "BH", universe=background,qvalueCutoff = 0.1, minGSSize = 10)
        REACT <- enrichPathway(gene=all_DEGs_ENTREZ, organism="human",pvalueCutoff = 0.1, pAdjustMethod = "BH", universe=background,qvalueCutoff=0.1, minGSSize=10)
        system("mkdir -p GO_Results")
        write.table(GO@result, file=paste0("./GO_Results/GO_", treatment, "_allDEGs.txt"), sep="\t", quote=FALSE, col.names=T)
        write.table(KEGG@result, file=paste0("./GO_Results/KEGG_", treatment, "_allDEGs.txt"), sep="\t", quote=FALSE, col.names=T)
        write.table(REACT@result, file=paste0("./GO_Results/REACTOME_", treatment, "_allDEGs.txt"), sep="\t", quote=FALSE, col.names=T)
        newGO <- GO@result[GO@result$Count>1 & GO@result$qvalue<0.1,]
        newKEGG <- KEGG@result[KEGG@result$Count>1 & KEGG@result$qvalue<0.1,]
        newREACT <- REACT@result[REACT@result$Count>1 & REACT@result$qvalue<0.1,]
        merged <- rbind(newGO, newKEGG, newREACT)
        system("mkdir -p GO_Results/significant")
        write.table(merged, file=paste0("./GO_Results/significant/", treatment, "_enrichment_summary_allDEGs.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=T)

}

