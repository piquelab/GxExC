# Anthony Findley

# This is Alan's script to do the variance partitioning of splicing on all cells together, and then on each cell type separately. I've modified it to ensure that our colors and labels are identical.
# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/All/Variance_Partition

library(tidyverse)
library(variancePartition)
library(limma)
library(edgeR)
library(data.table)
library(doParallel)

cl <- makeCluster(8)
registerDoParallel(cl)

# I can skip the first part and load in the object he saved in line 25
#CM <- fread("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/CM_juncs.txt.gz") %>% select(-V1)
#IPSC <- fread("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/IPSC_juncs.txt.gz") %>% select(-V1)
#LCL <- fread("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/LCL_juncs.txt.gz") %>% select(-V1)

#geneCounts <- full_join(CM, IPSC, by = "Intron") %>% full_join(LCL, by = "Intron") %>%
#  column_to_rownames("Intron") %>% as.matrix()

#geneCounts[is.na(geneCounts)] <- 0
#save(geneCounts, "geneCounts_allcells_full.RData")
load("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/geneCounts_allcells_full.RData")

# identify genes that pass expression cutoff.
isexpr <- rowSums(cpm(geneCounts)>1) >= 0.5 * ncol(geneCounts)

# create data structure with only expressed genes.
gExpr <- DGEList(counts=geneCounts[isexpr,])

# Perform TMM normalization
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of uncertainty.
# Recommend including variables with a small number of categories that explain a substantial amount of variation.
meta <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/VarPart_meta.txt", header = T)

# I need to make meta have names which match gExpr
meta$Filename <- paste0(meta$File, ".junc")
rownames(meta) <- meta$Filename

meta <- meta[order(match(meta$Filename, colnames(gExpr$counts))),] # Puts meta in same order as gExpr
design <- model.matrix( ~ CellType, meta)

# Estimate precision weights for each gene and sample.
# This models uncertainty in expression measurements.
vobjGenes <- voom(gExpr, design)

# Define formula.
form <- ~ (1|Individual) + (1|CellType) + (1|Plate.ID) + (1|Treatment.Name)

# variancePartition seamlessly deals with the result of voom() by default, it seamlessly models the precision weights.
# This can be turned off with useWeights=FALSE.
varPart <- fitExtractVarPartModel(vobjGenes, form, meta)

# It took forever to get the varPart object, so save it
write.table(varPart, file="splicing_allCells.varPart_corrected.tab", row.names=T, col.names=T, quote=F, sep="\t")

summary(varPart)
   Cell.Type        Individual        Treatment            Plate
 Min.   :0.0000   Min.   :0.00000   Min.   :0.000000   Min.   :0.00000
 1st Qu.:0.4678   1st Qu.:0.01389   1st Qu.:0.001823   1st Qu.:0.01531
 Median :0.7032   Median :0.03235   Median :0.007075   Median :0.03707
 Mean   :0.6393   Mean   :0.05927   Mean   :0.014836   Mean   :0.05741
 3rd Qu.:0.8525   3rd Qu.:0.07131   3rd Qu.:0.018012   3rd Qu.:0.07713
 Max.   :0.9999   Max.   :0.97977   Max.   :0.628445   Max.   :0.68228
   Residuals
 Min.   :0.000023
 1st Qu.:0.086833
 Median :0.181788
 Mean   :0.229226
 3rd Qu.:0.331942
 Max.   :0.996374

 # For how many transcripts is individual component the greatest:
sum(dat$Individual > dat$Cell.Type & dat$Individual > dat$Plate & dat$Individual > dat$Treatment & dat$Individual > dat$Residuals) # 1274

sum(dat$Treatment > dat$Cell.Type & dat$Treatment > dat$Plate & dat$Treatment > dat$Individual & dat$Treatment > dat$Residuals) # 50

# Make the plot look the way I want it
colnames(varPart) <- c("Cell Type", "Individual", "Plate", "Treatment", "Residuals")
varPart <- varPart[,c("Cell Type", "Individual", "Treatment", "Plate", "Residuals")]

# --------------------------------------

# Color scheme is:
	# Cell = #e41a1c
	# Individual = #377eb8
	# Treatment = #4daf4a
	# Individual x Treatment (not in this plot) = #984ea3
	# Plate = #ff7f00
	# Residuals = grey85

# violin plot of contribution of each variable to total variance.
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/splicing_allCells_corrected.VarPlot.pdf")
plotVarPart(varPart, reorder=FALSE, col=c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "grey85")) +
  theme(axis.text.x=element_text(size=18), axis.text.y = element_text(size = 16)) #+
 # stat_summary(geom = "text", vjust = 1.1, fontface = "bold", size = 7)
dev.off()


##################################################
##################################################

# Make the gene expression variance partitioning plots look the same
load("VarPart_R2_R3_corrected.Rd")

colnames(vp) <- c("Cell Type", "Individual", "Plate", "Treatment", "Residuals")
vp <- vp[,c("Cell Type", "Individual", "Treatment", "Plate", "Residuals")]

# violin plot of contribution of each variable to total variance.
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/final_GE_allCells_corrected.VarPlot.pdf")
plotVarPart(vp, reorder=FALSE, col=c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "grey85")) +
  theme(axis.text.x=element_text(size=18), axis.text.y = element_text(size = 16)) #+
 # stat_summary(geom = "text", vjust = 1.1, fontface = "bold", size = 7)
dev.off()

###
# Now do gene expression for the cells analyzed separately

load("VarPart_R2_R3_CellTypes.testedForSplice_corrected.Rd")

colnames(vp_LCL) <- c("Individual", "Plate", "Treatment", "Ind:Treat", "Residuals")
vp_LCL <- vp_LCL[ ,c("Individual", "Treatment", "Ind:Treat", "Plate", "Residuals")]

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/final_GE_LCL_corrected.VarPlot.pdf")
plotVarPart(vp_LCL, reorder=FALSE, col=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "grey85")) +
  theme(axis.text.x=element_text(size=18), axis.text.y = element_text(size = 16)) #+
 # stat_summary(geom = "text", vjust = 1.1, fontface = "bold", size = 7)
dev.off()

# CMs
colnames(vp_CM) <- c("Individual", "Ind:Treat", "Plate", "Treatment", "Residuals")
vp_CM <- vp_CM[ ,c("Individual", "Treatment", "Ind:Treat", "Plate", "Residuals")]

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/final_GE_CM_corrected.VarPlot.pdf")
plotVarPart(vp_CM, reorder=FALSE, col=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "grey85")) +
  theme(axis.text.x=element_text(size=18), axis.text.y = element_text(size = 16)) #+
 # stat_summary(geom = "text", vjust = 1.1, fontface = "bold", size = 7)
dev.off()

# IPSCs
colnames(vp_IPSC) <- c("Individual", "Ind:Treat", "Treatment", "Plate", "Residuals")
vp_IPSC <- vp_IPSC[ ,c("Individual", "Treatment", "Ind:Treat", "Plate", "Residuals")]

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/final_GE_IPSC_corrected.VarPlot.pdf")
plotVarPart(vp_IPSC, reorder=FALSE, col=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "grey85")) +
  theme(axis.text.x=element_text(size=18), axis.text.y = element_text(size = 16)) #+
 # stat_summary(geom = "text", vjust = 1.1, fontface = "bold", size = 7)
dev.off()

### Francesca wants to know for how many genes the Treat or Ind:Treat variables explain the most variance, after removing residuals
vp_LCL %>% select(-Residuals) -> vp_LCL_noResid
vp_LCL_noResid_percent <- vp_LCL_noResid / rowSums(vp_LCL_noResid)
colnames(vp_LCL_noResid_percent) <- c("Plate", "Treat", "Ind", "Ind_Treat")

sum(vp_LCL_noResid_percent$Treat + vp_LCL_noResid_percent$Ind_Treat < 0.5, na.rm=T) # 54569
sum(vp_LCL_noResid_percent$Treat + vp_LCL_noResid_percent$Ind_Treat > 0.5, na.rm=T) # 32780

# CMs
vp_CM %>% select(-Residuals) -> vp_CM_noResid
vp_CM_noResid_percent <- vp_CM_noResid / rowSums(vp_CM_noResid)
colnames(vp_CM_noResid_percent) <- c("Plate", "Treat", "Ind", "Ind_Treat")

sum(vp_CM_noResid_percent$Treat + vp_CM_noResid_percent$Ind_Treat < 0.5, na.rm=T) # 72120
sum(vp_CM_noResid_percent$Treat + vp_CM_noResid_percent$Ind_Treat > 0.5, na.rm=T) # 18607

# IPSCs
vp_IPSC %>% select(-Residuals) -> vp_IPSC_noResid
vp_IPSC_noResid_percent <- vp_IPSC_noResid / rowSums(vp_IPSC_noResid)
colnames(vp_IPSC_noResid_percent) <- c("Plate", "Treat", "Ind", "Ind_Treat")

sum(vp_IPSC_noResid_percent$Treat + vp_IPSC_noResid_percent$Ind_Treat < 0.5, na.rm=T) # 45701
sum(vp_IPSC_noResid_percent$Treat + vp_IPSC_noResid_percent$Ind_Treat > 0.5, na.rm=T) # 45458

# Total percent of genes for which treat + ind_treat explains half of variance after removing residual:
(54569 + 72120 + 45701) / (54569 + 72120 + 45701 + 32780 + 18607 + 45458)
[1] 0.6402957

######################################
######################################
# Now I have to repeat Alan's splicing within each cell type

# CARDIOMYOCYTES
CM <- fread("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/CM_juncs.txt.gz") %>% select(-V1) %>% column_to_rownames("Intron") %>% as.matrix()

# identify genes that pass expression cutoff.
isexpr <- rowSums(cpm(CM)>1) >= 0.5 * ncol(CM)

# create data structure with only expressed genes.
gExpr <- DGEList(counts=CM[isexpr,])

# Perform TMM normalization
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of uncertainty.
# Recommend including variables with a small number of categories that explain a substantial amount of variation.
meta <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/VarPart_CMmeta.txt", header = T)

# I need to make meta have names which match gExpr
meta$Filename <- paste0(meta$File, ".junc")
rownames(meta) <- meta$Filename

meta <- meta[order(match(meta$Filename, colnames(gExpr$counts))),] # Puts meta in same order as gExpr

design <- model.matrix( ~ Individual, meta)

# Estimate precision weights for each gene and sample.
# This models uncertainty in expression measurements.
vobjGenes <- voom(gExpr, design)

# Define formula.
form <- ~ (1|Individual) + (1|Plate.ID) + (1|Treatment.Name) + (1|Individual:Treatment.Name)

# variancePartition seamlessly deals with the result of voom() by default, it seamlessly models the precision weights.
# This can be turned off with useWeights=FALSE.
varPart <- fitExtractVarPartModel(vobjGenes, form, meta)

vp_splicing_CM <- varPart

colnames(vp_splicing_CM) <- c("Individual", "Ind:Treat", "Plate", "Treatment", "Residuals")
vp_splicing_CM <- vp_splicing_CM[ ,c("Individual", "Treatment", "Ind:Treat", "Plate", "Residuals")]

   Individual       Treatment          Ind:Treat            Plate
 Min.   :0.0000   Min.   :0.000000   Min.   :0.000000   Min.   :0.00000
 1st Qu.:0.2456   1st Qu.:0.005895   1st Qu.:0.000000   1st Qu.:0.02791
 Median :0.4473   Median :0.021706   Median :0.006858   Median :0.08135
 Mean   :0.4621   Mean   :0.038134   Mean   :0.023357   Mean   :0.11525
 3rd Qu.:0.6749   3rd Qu.:0.050222   3rd Qu.:0.033889   3rd Qu.:0.16764
 Max.   :0.9963   Max.   :0.870847   Max.   :0.454787   Max.   :0.99907
   Residuals
 Min.   :0.00075
 1st Qu.:0.19280
 Median :0.34047
 Mean   :0.36111
 3rd Qu.:0.50941
 Max.   :1.00000

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/final_splicing_CM_corrected.VarPlot.pdf")
plotVarPart(vp_splicing_CM, reorder=FALSE, col=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "grey85")) +
  theme(axis.text.x=element_text(size=18), axis.text.y = element_text(size = 16))
dev.off()

# -----------------------------------------------------------


# IPSC
IPSC <- fread("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/IPSC_juncs.txt.gz") %>% select(-V1) %>% column_to_rownames("Intron") %>% as.matrix()

# identify genes that pass expression cutoff.
isexpr <- rowSums(cpm(IPSC)>1) >= 0.5 * ncol(IPSC)

# create data structure with only expressed genes.
gExpr <- DGEList(counts=IPSC[isexpr,])

# Perform TMM normalization
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of uncertainty.
# Recommend including variables with a small number of categories that explain a substantial amount of variation.
meta <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/VarPart_IPSCmeta.txt", header = T)

# I need to make meta have names which match gExpr
meta$Filename <- paste0(meta$File, ".junc")
rownames(meta) <- meta$Filename

meta <- meta[order(match(meta$Filename, colnames(gExpr$counts))),] # Puts meta in same order as gExpr

design <- model.matrix( ~ Individual, meta)

# Estimate precision weights for each gene and sample.
# This models uncertainty in expression measurements.
vobjGenes <- voom(gExpr, design)

# Define formula.
form <- ~ (1|Individual) + (1|Plate.ID) + (1|Treatment.Name) + (1|Individual:Treatment.Name)

# variancePartition seamlessly deals with the result of voom() by default, it seamlessly models the precision weights.
# This can be turned off with useWeights=FALSE.
varPart <- fitExtractVarPartModel(vobjGenes, form, meta)

vp_splicing_IPSC <- varPart

colnames(vp_splicing_IPSC) <- c("Individual", "Ind:Treat", "Plate", "Treatment", "Residuals")
vp_splicing_IPSC <- vp_splicing_IPSC[ ,c("Individual", "Treatment", "Ind:Treat", "Plate", "Residuals")]


pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/final_splicing_IPSC_corrected.VarPlot.pdf")
plotVarPart(vp_splicing_IPSC, reorder=FALSE, col=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "grey85")) +
  theme(axis.text.x=element_text(size=18), axis.text.y = element_text(size = 16))
dev.off()

summary(vp_splicing_IPSC)
   Individual        Treatment         Ind:Treat             Plate
 Min.   :0.00000   Min.   :0.00000   Min.   :0.0000000   Min.   :0.00000
 1st Qu.:0.09698   1st Qu.:0.02296   1st Qu.:0.0000000   1st Qu.:0.01771
 Median :0.20272   Median :0.06085   Median :0.0000000   Median :0.06135
 Mean   :0.26305   Mean   :0.08744   Mean   :0.0092118   Mean   :0.09726
 3rd Qu.:0.38160   3rd Qu.:0.12018   3rd Qu.:0.0002183   3rd Qu.:0.14048
 Max.   :0.99630   Max.   :0.79057   Max.   :0.4963936   Max.   :0.94686
   Residuals
 Min.   :0.001547
 1st Qu.:0.389061
 Median :0.555948
 Mean   :0.543029
 3rd Qu.:0.705262
 Max.   :1.000000


# -----------------------------------------------------------


# LCL
LCL <- fread("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/LCL_juncs.txt.gz") %>% select(-V1) %>% column_to_rownames("Intron") %>% as.matrix()

# identify genes that pass expression cutoff.
isexpr <- rowSums(cpm(LCL)>1) >= 0.5 * ncol(LCL)

# create data structure with only expressed genes.
gExpr <- DGEList(counts=LCL[isexpr,])

# Perform TMM normalization
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of uncertainty.
# Recommend including variables with a small number of categories that explain a substantial amount of variation.
meta <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/VarPart_LCLmeta.txt", header = T)

# I need to make meta have names which match gExpr
meta$Filename <- paste0(meta$File, ".junc")
rownames(meta) <- meta$Filename

meta <- meta[order(match(meta$Filename, colnames(gExpr$counts))),] # Puts meta in same order as gExpr

design <- model.matrix( ~ Individual, meta)

# Estimate precision weights for each gene and sample.
# This models uncertainty in expression measurements.
vobjGenes <- voom(gExpr, design)

# Define formula.
form <- ~ (1|Individual) + (1|Plate.ID) + (1|Treatment.Name) + (1|Individual:Treatment.Name)

# variancePartition seamlessly deals with the result of voom() by default, it seamlessly models the precision weights.
# This can be turned off with useWeights=FALSE.
varPart <- fitExtractVarPartModel(vobjGenes, form, meta)

vp_splicing_LCL <- varPart

colnames(vp_splicing_LCL) <- c("Individual", "Ind:Treat", "Plate", "Treatment", "Residuals")
vp_splicing_LCL <- vp_splicing_LCL[ ,c("Individual", "Treatment", "Ind:Treat", "Plate", "Residuals")]


pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/final_splicing_LCL_corrected.VarPlot.pdf")
plotVarPart(vp_splicing_LCL, reorder=FALSE, col=c("#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "grey85")) +
  theme(axis.text.x=element_text(size=18), axis.text.y = element_text(size = 16))
dev.off()

summary(vp_splicing_LCL)
   Individual       Treatment         Ind:Treat           Plate
 Min.   :0.0000   Min.   :0.00000   Min.   :0.00000   Min.   :0.000000
 1st Qu.:0.1361   1st Qu.:0.05853   1st Qu.:0.00000   1st Qu.:0.008622
 Median :0.2811   Median :0.11770   Median :0.03126   Median :0.037539
 Mean   :0.3322   Mean   :0.14738   Mean   :0.04933   Mean   :0.066230
 3rd Qu.:0.4887   3rd Qu.:0.20452   3rd Qu.:0.07737   3rd Qu.:0.092816
 Max.   :0.9967   Max.   :0.88422   Max.   :0.44404   Max.   :0.769883
   Residuals
 Min.   :0.001445
 1st Qu.:0.239116
 Median :0.387449
 Mean   :0.404880
 3rd Qu.:0.556604
 Max.   :1.000000


# Save the splicing objects
save(list = c("vp_splicing_LCL", "vp_splicing_IPSC", "vp_splicing_CM"), file = "varPart_cellSep_Splice_corrected.RData")

