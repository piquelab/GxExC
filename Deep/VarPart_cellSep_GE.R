# Anthony Findley

# Run variance partition on CMs, IPSCs, and IPSCs separately (NovaSeq deep sequencing only). I'm only keeping transcripts in genes for which Alan could measure splicing. I've updated this so that the samples in the covariate file are in the same order as the gene expression object. This changes the results dramatically (much less residual).
# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/All/Variance_Partition

library(variancePartition)
library(doParallel)
library(limma)
library(edgeR)
library(tidyverse)

cl <- makeCluster(8)
registerDoParallel(cl)

# Annotation file
anno <- read.table("/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.tsv.gz", header=T, sep="\t", stringsAsFactors=F)

##### First run variance partition without removing potential confounders

# Read data in 
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM2R2/R2_R3/fpkm/GeneCounts/DCM2R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM1R2/R2_R3/fpkm/GeneCounts/DCM1R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM2R1/R2_R3/fpkm/GeneCounts/DCM2R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM1R1/R2_R3/fpkm/GeneCounts/DCM1R1.data.Rd")

data <- cbind(DCM1R1counts, DCM2R1counts, DCM1R2counts, DCM2R2counts)

#Add covariate information to dge$samples
cv5 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DCM1R1_covar.txt", as.is=T,header=T)
cv6 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DCM1R2_covar.txt", as.is=T,header=T)
cv7 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DCM2R1_covar.txt", as.is=T,header=T)
cv8 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DCM2R2_covar.txt", as.is=T,header=T)

cv5 <- cv5[,c(1:10,14,15)]
cv6 <- cv6[,c(1:10,14,15)]

# Remove water control from DCM1R2
cv6 <- cv6[-which(cv6$Treatment.ID == "CO1"),]
cv6$Control.ID <- "CO2"

cv <- rbind(cv5,cv6,cv7,cv8)
data <- data[,which(colnames(data) %in% cv$Filename)]
cv <- cv[which(cv$Filename %in% colnames(data)),]

# Load Alan's splice genes and select genes in the annotation file tested by Alan
Alan_CM <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/CM_tested_genes_unique.txt", stringsAsFactors=F, header=F)
anno_CM <- anno[which(anno$Gene.name %in% Alan_CM$V1),]
data_sub <- data[which(rownames(data) %in% anno_CM$Transcript.stable.ID), ]

# VariancePartition's reccomendation for expressed genes
isexpr <- rowSums(cpm(data_sub)>1) >= 0.5 * ncol(data_sub)

# Normalize counts
gExpr <- DGEList(counts=data_sub[isexpr,])
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of
# uncertainty.
# Recommend including variables with a small number of categories
# that explain a substantial amount of variation
cv <- cv[order(match(cv$Filename, colnames(gExpr$counts))),] # Puts cv in same order as gExpr

design <- model.matrix( ~ Plate.ID, cv)

# Estimate precision weights for each gene and sample
# This models uncertainty in expression measurements
vobjGenes <- voom(gExpr, design )

# Define formula
form <- ~ (1|Individual) + (1|Plate.ID) + (1|Treatment.Name) + (1|Individual:Treatment.Name)

# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
varPart <- fitExtractVarPartModel( vobjGenes, form, cv)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )

# Now save results so they don't get overwritten
vp_CM <- vp
data_CM <- data

# violin plot of contribution of each variable to total variance
vp_CM <- vp_CM[,c(1,4,2,3,5)] # Make order of bars same as Alan

summary(vp_CM)
   Individual     Individual:Treatment.Name    Plate.ID       Treatment.Name
 Min.   :0.0000   Min.   :0.00000           Min.   :0.00000   Min.   :0.00000
 1st Qu.:0.3504   1st Qu.:0.00000           1st Qu.:0.04526   1st Qu.:0.01528
 Median :0.5536   Median :0.01378           Median :0.12551   Median :0.03610
 Mean   :0.5406   Mean   :0.03053           Mean   :0.16338   Mean   :0.05718
 3rd Qu.:0.7425   3rd Qu.:0.04121           3rd Qu.:0.25332   3rd Qu.:0.07514
 Max.   :0.9941   Max.   :0.45771           Max.   :0.81022   Max.   :0.69760
   Residuals
 Min.   :0.004296
 1st Qu.:0.109834
 Median :0.186899
 Mean   :0.208270
 3rd Qu.:0.283168
 Max.   :0.956710

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/CMs_testedForSplice_corrected.VarPlot.pdf")
plotVarPart( vp_CM )
dev.off()

#################################################################################################
#################################################################################################

# Now do IPSCs

##### First run variance partition without removing potential confounders

# Read data in 
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC2R2/P1_R2/fpkm/GeneCounts/DIPSC2R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC1R2/P1_R2/fpkm/GeneCounts/DIPSC1R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC2R1/P1_R2/fpkm/GeneCounts/DIPSC2R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC1R1/P1_R2/fpkm/GeneCounts/DIPSC1R1.data.Rd")

data <- cbind(DIPSC1R1counts, DIPSC2R1counts, DIPSC1R2counts, DIPSC2R2counts)

#Add covariate information to dge$samples
cv1 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DIPSC1R1_covar.txt", as.is=T,header=T)
cv2 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DIPSC1R2_covar.txt", as.is=T,header=T)
cv3 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DIPSC2R1_covar.txt", as.is=T,header=T)
cv4 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DIPSC2R2_covar.txt", as.is=T,header=T)

cv <- rbind(cv1,cv2,cv3,cv4)
data <- data[,which(colnames(data) %in% cv$Filename)]
cv <- cv[which(cv$Filename %in% colnames(data)),]

# Load Alan's splice genes and select genes in the annotation file tested by Alan
Alan_IPSC <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/IPSC_tested_genes_unique.txt", stringsAsFactors=F, header=F)
anno_IPSC <- anno[which(anno$Gene.name %in% Alan_IPSC$V1),]
data_sub <- data[which(rownames(data) %in% anno_IPSC$Transcript.stable.ID), ]

# VariancePartition's reccomendation for expressed genes
isexpr <- rowSums(cpm(data_sub)>1) >= 0.5 * ncol(data_sub)

# Normalize counts
gExpr <- DGEList(counts=data_sub[isexpr,])
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of
# uncertainty.
# Recommend including variables with a small number of categories
# that explain a substantial amount of variation
cv <- cv[order(match(cv$Filename, colnames(gExpr$counts))),] # Puts cv in same order as gExpr

design <- model.matrix( ~ Plate.ID, cv)

# Estimate precision weights for each gene and sample
# This models uncertainty in expression measurements
vobjGenes <- voom(gExpr, design )

# Define formula
form <- ~ (1|Individual) + (1|Plate.ID) + (1|Treatment.Name) + (1|Individual:Treatment.Name)

# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
varPart <- fitExtractVarPartModel( vobjGenes, form, cv)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )

# Now save results so they don't get overwritten
vp_IPSC <- vp
data_IPSC <- data

summary(vp_IPSC)
   Individual     Treatment.Name       Plate.ID       Individual:Treatment.Name
 Min.   :0.0000   Min.   :0.00000   Min.   :0.00000   Min.   :0.000000
 1st Qu.:0.1423   1st Qu.:0.05056   1st Qu.:0.04109   1st Qu.:0.000000
 Median :0.2784   Median :0.11188   Median :0.11028   Median :0.000000
 Mean   :0.3287   Mean   :0.14612   Mean   :0.14908   Mean   :0.006128
 3rd Qu.:0.4776   3rd Qu.:0.20419   3rd Qu.:0.22370   3rd Qu.:0.000000
 Max.   :0.9957   Max.   :0.84356   Max.   :0.99808   Max.   :0.337417
   Residuals
 Min.   :0.001101
 1st Qu.:0.230421
 Median :0.350207
 Mean   :0.369999
 3rd Qu.:0.489057
 Max.   :1.000000

# violin plot of contribution of each variable to total variance
vp_IPSC <- vp_IPSC[,c(1,4,2,3,5)] # Make order of bars same as Alan

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/IPSCs_testedForSplice_corrected.VarPlot.pdf")
plotVarPart( vp_IPSC, reorder=FALSE )
dev.off()

#################################################################################################
#################################################################################################

# Now do LCLs

##### First run variance partition without removing potential confounders

# Read data in 
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL2R2/P1_R2/fpkm/GeneCounts/DLCL2R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL1R2/P1_R2/fpkm/GeneCounts/DLCL1R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL2R1/P1_R2/fpkm/GeneCounts/DLCL2R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL1R1/P1_R2/fpkm/GeneCounts/DLCL1R1.data.Rd")

data <- cbind(DLCL1R1counts, DLCL2R1counts, DLCL1R2counts, DLCL2R2counts)

#Add covariate information to dge$samples
cv1 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL1R1_covar.txt", as.is=T,header=T)
cv2 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL1R2_covar.txt", as.is=T,header=T)
cv3 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL2R1_covar.txt", as.is=T,header=T)
cv4 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL2R2_covar.txt", as.is=T,header=T)

cv <- rbind(cv1,cv2,cv3,cv4)
data <- data[,which(colnames(data) %in% cv$Filename)]
cv <- cv[which(cv$Filename %in% colnames(data)),]

# Load Alan's splice genes and select genes in the annotation file tested by Alan
Alan_LCL <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/LCL_tested_genes_unique.txt", stringsAsFactors=F, header=F)
anno_LCL <- anno[which(anno$Gene.name %in% Alan_LCL$V1),]
data_sub <- data[which(rownames(data) %in% anno_LCL$Transcript.stable.ID), ]

# VariancePartition's reccomendation for expressed genes
isexpr <- rowSums(cpm(data_sub)>1) >= 0.5 * ncol(data_sub)

# Normalize counts
gExpr <- DGEList(counts=data_sub[isexpr,])
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of
# uncertainty.
# Recommend including variables with a small number of categories
# that explain a substantial amount of variation
cv <- cv[order(match(cv$Filename, colnames(gExpr$counts))),] # Puts cv in same order as gExpr

design <- model.matrix( ~ Plate.ID, cv)

# Estimate precision weights for each gene and sample
# This models uncertainty in expression measurements
vobjGenes <- voom(gExpr, design )

# Define formula
form <- ~ (1|Individual) + (1|Plate.ID) + (1|Treatment.Name) + (1|Individual:Treatment.Name)

# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
varPart <- fitExtractVarPartModel( vobjGenes, form, cv)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )

# Now save results so they don't get overwritten
vp_LCL <- vp
data_LCL <- data

summary(vp_LCL)
   Individual     Treatment.Name   Individual:Treatment.Name    Plate.ID
 Min.   :0.0000   Min.   :0.0000   Min.   :0.000000          Min.   :0.000000
 1st Qu.:0.1939   1st Qu.:0.1196   1st Qu.:0.004943          1st Qu.:0.009764
 Median :0.3630   Median :0.2185   Median :0.040505          Median :0.033390
 Mean   :0.3985   Mean   :0.2542   Mean   :0.057279          Mean   :0.056759
 3rd Qu.:0.5816   3rd Qu.:0.3569   3rd Qu.:0.088950          3rd Qu.:0.078197
 Max.   :0.9937   Max.   :0.8980   Max.   :0.478050          Max.   :0.673480
   Residuals
 Min.   :0.002424
 1st Qu.:0.103873
 Median :0.191240
 Mean   :0.233267
 3rd Qu.:0.324349
 Max.   :0.972251


# violin plot of contribution of each variable to total variance
vp_LCL <- vp_LCL[,c(1,4,2,3,5)] # Make order of bars same as Alan

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/LCLs_testedForSplice_corrected.VarPlot.pdf")
plotVarPart( vp_LCL, reorder=FALSE )
dev.off()
