# Anthony Findley


# Run variance partition on IPSCs, IPSCs, and CMs (Novogene deep sequencing only; both runs for CMs). Updated to ensure covariate file is in correct order for gene expression count matrix
# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/All/Variance_Partition

library(variancePartition)
library(doParallel)
library(limma)
library(edgeR)
library(tidyverse)

##### First run variance partition without removing potential confounders

# Read data in 
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL1R1/P1_R2/fpkm/GeneCounts/DLCL1R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL2R1/P1_R2/fpkm/GeneCounts/DLCL2R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL1R2/P1_R2/fpkm/GeneCounts/DLCL1R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL2R2/P1_R2/fpkm/GeneCounts/DLCL2R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM2R2/R2_R3/fpkm/GeneCounts/DCM2R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM1R2/R2_R3/fpkm/GeneCounts/DCM1R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM2R1/R2_R3/fpkm/GeneCounts/DCM2R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM1R1/R2_R3/fpkm/GeneCounts/DCM1R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC1R1/P1_R2/fpkm/GeneCounts/DIPSC1R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC2R1/P1_R2/fpkm/GeneCounts/DIPSC2R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC1R2/P1_R2/fpkm/GeneCounts/DIPSC1R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC2R2/P1_R2/fpkm/GeneCounts/DIPSC2R2.data.Rd")

data <- cbind(DLCL1R1counts, DLCL2R1counts, DLCL2R2counts, DLCL1R2counts, DCM1R1counts, DCM2R1counts, DCM1R2counts, DCM2R2counts, DIPSC1R1counts, DIPSC2R1counts, DIPSC2R2counts, DIPSC1R2counts)

#Add covariate information to dge$samples
cv1 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL2R1_covar.txt", as.is=T,header=T)
cv2 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL1R1_covar.txt", as.is=T,header=T)
cv3 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL1R2_covar.txt", as.is=T,header=T)
cv4 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL2R2_covar.txt", as.is=T,header=T)
cv5 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DCM1R1_covar.txt", as.is=T,header=T)
cv6 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DCM1R2_covar.txt", as.is=T,header=T)
cv7 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DCM2R1_covar.txt", as.is=T,header=T)
cv8 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DCM2R2_covar.txt", as.is=T,header=T)
cv9 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DIPSC2R1_covar.txt", as.is=T,header=T)
cv10 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DIPSC1R1_covar.txt", as.is=T,header=T)
cv11 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DIPSC1R2_covar.txt", as.is=T,header=T)
cv12 <- read.table("~/rprdata/Anthony/GxE_iPSC/Deep/Docs/DIPSC2R2_covar.txt", as.is=T,header=T)

cv5 <- cv5[,c(1:10,14,15)]
cv6 <- cv6[,c(1:10,14,15)]

# Remove water control from DCM1R2
cv6 <- cv6[-which(cv6$Treatment.ID == "CO1"),]
cv6$Control.ID <- "CO2"

cv <- rbind(cv1,cv2,cv3,cv4,cv5,cv6,cv7,cv8,cv9,cv10,cv11,cv12)
data <- data[,which(colnames(data) %in% cv$Filename)]
cv <- cv[which(cv$Filename %in% colnames(data)),]

# VariancePartition's recommendation for expressed genes
isexpr <- rowSums(cpm(data)>1) >= 0.5 * ncol(data)

# Normalize counts
gExpr <- DGEList(counts=data[isexpr,])
gExpr <- calcNormFactors(gExpr)

# Specify variables to be included in the voom() estimates of
# uncertainty.
# Recommend including variables with a small number of categories
# that explain a substantial amount of variation
cv <- cv[order(match(cv$Filename, colnames(gExpr$counts))),] # Puts cv in same order as gExpr

design <- model.matrix( ~ CellType, cv)

# Estimate precision weights for each gene and sample
# This models uncertainty in expression measurements
vobjGenes <- voom(gExpr, design )

# Define formula
form <- ~ (1|Individual) + (1|CellType) + (1|Plate.ID) + (1|Treatment.Name)

# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
varPart <- fitExtractVarPartModel( vobjGenes, form, cv)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp <- sortCols( varPart )

# violin plot of contribution of each variable to total variance
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/VarPlot1_corrected.pdf")
plotVarPart( vp )
dev.off()

summary(vp)
   Cell Type        Individual        Treatment            Plate
 Min.   :0.0000   Min.   :0.00000   Min.   :0.000000   Min.   :0.00000
 1st Qu.:0.5238   1st Qu.:0.01701   1st Qu.:0.003329   1st Qu.:0.01513
 Median :0.7432   Median :0.04101   Median :0.011072   Median :0.03738
 Mean   :0.6736   Mean   :0.07221   Mean   :0.023468   Mean   :0.06119
 3rd Qu.:0.8763   3rd Qu.:0.09241   3rd Qu.:0.027886   3rd Qu.:0.08112
 Max.   :0.9994   Max.   :0.91470   Max.   :0.575272   Max.   :0.60201
   Residuals
 Min.   :0.0000478
 1st Qu.:0.0619223
 Median :0.1284727
 Mean   :0.1695756
 3rd Qu.:0.2436659
 Max.   :0.8375062

# For how many transcripts is individual component the greatest:
sum(vp$Individual > vp$CellType & vp$Individual > vp$Plate.ID & vp$Individual > vp$Treatment.Name & vp$Individual > vp$Residuals) # 1754

sum(vp$Treatment.Name > vp$CellType & vp$Treatment.Name > vp$Plate.ID & vp$Treatment.Name > vp$Individual & vp$Treatment.Name > vp$Residuals) # 80

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/final_GE_allCells_corrected_2.VarPlot.pdf")
plotVarPart(vp, reorder=FALSE, col=c("#e41a1c", "#377eb8", "#4daf4a", "#ff7f00", "grey85")) +
  theme(axis.text.x=element_text(size=18), axis.text.y = element_text(size = 16)) #+
 # stat_summary(geom = "text", vjust = 1.1, fontface = "bold", size = 7)
dev.off()


######## Now try adding the interaction effect to the formula
# Define formula
form2 <- ~ (1|Individual) + (1|CellType) + (1|Plate.ID) + (1|Treatment.Name) + (1|Treatment.Name:CellType)

# variancePartition seamlessly deals with the result of voom()
# by default, it seamlessly models the precision weights
# This can be turned off with useWeights=FALSE
varPart2 <- fitExtractVarPartModel( vobjGenes, form2, cv)

# sort variables (i.e. columns) by median fraction
# of variance explained
vp2 <- sortCols( varPart2 )

# violin plot of contribution of each variable to total variance
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/VarPlot2_corrected.pdf")
plotVarPart( vp2 )
dev.off()

# Save an image because it took a while to run
save.image("VarPart_R2_R3_corrected.Rd")

######################################

library(irlba)

p1 <- prcomp_irlba(vobjGenes$E, n=15)

pca <- t(p1$rotation)

pcavar <- (p1$sdev)^2

pcaprop <- pcavar/sum(pcavar)


varPart2PCA <- fitExtractVarPartModel( pca, form2, cv)

varPartSummary = t(as.matrix(varPart2PCA)) %*% (pcaprop)

varPartSummary

varPartSummary[-1]/sum(varPartSummary[-1])*100 
 

varPartSummary[-(1:3)]/sum(varPartSummary[-(1:3)])*100 


#######################################################

# Get transcripts expressed in all tissues

anno <- read.table("/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.tsv.gz", header=T, sep="\t", stringsAsFactors=F)

treatments <- c("T12C1", "T13C1", "T14C1", "T15C1", "T19C1", "T20C1", "T27C1", "T30C1", "T33C1", "T42C1", "T6C1", "T9C1")
i=1

CM <- read.table(paste0("/nfs/rprdata/Anthony/GxE_iPSC/Deep/CM_combined/P1_R2/out_data_DCM_P1R2_04.19.19/stats/DCM_P1R2_04.19.19_DEG_stats_", treatments[i], ".txt"), stringsAsFactors=F, header=T)
IPSC <- read.table(paste0("/nfs/rprdata/Anthony/GxE_iPSC/Deep/IPSC_combined/P1_R2/out_data_DIPSC_P1R2_04.19.19/stats/DIPSC_P1R2_04.19.19_DEG_stats_", treatments[i], ".txt"), stringsAsFactors=F, header=T)
LCL <- read.table(paste0("/nfs/rprdata/Anthony/GxE_iPSC/Deep/LCL_combined/P1_R2/out_data_DLCL_P1R2_04.19.19/stats/DLCL_P1R2_04.19.19_DEG_stats_", treatments[i], ".txt"), stringsAsFactors=F, header=T)

CM <- merge(CM, anno[, c("Transcript.stable.ID", "Gene.stable.ID", "Gene.name", "Gene.type")], by.x = "t.id", by.y = "Transcript.stable.ID")
colnames(CM) <- colnames(IPSC)

CM_DE <- unique(CM[which(CM$padj < 0.1 & (CM$logFC > 0.25 | CM$logFC < -0.25)), "g.id"])
IPSC_DE <- unique(IPSC[which(IPSC$padj < 0.1 & (IPSC$logFC > 0.25 | IPSC$logFC < -0.25)), "g.id"])
LCL_DE <- unique(LCL[which(LCL$padj < 0.1 & (LCL$logFC > 0.25 | LCL$logFC < -0.25)), "g.id"])

shared_t.id <- IPSC$t.id[which(IPSC$t.id %in% CM$t.id[which(CM$t.id %in% LCL$t.id)])]

# Normalize counts
gExpr_shared <- DGEList(counts=data[which(rownames(data) %in% shared_t.id),])
gExpr_shared <- calcNormFactors(gExpr_shared)

# Specify variables to be included in the voom() estimates of
# uncertainty.
# Recommend including variables with a small number of categories
# that explain a substantial amount of variation
cv <- cv[order(match(cv$Filename, colnames(gExpr_shared$counts))),] # Puts cv in same order as gExpr
design <- model.matrix( ~ CellType, cv)

# Estimate precision weights for each gene and sample
# This models uncertainty in expression measurements
vobjGenes_shared <- voom(gExpr_shared, design )

p1_shared <- prcomp_irlba(vobjGenes_shared$E, n=15)
pca_shared <- t(p1_shared$rotation)

pcavar_shared <- (p1_shared$sdev)^2
pcaprop_shared <- pcavar_shared/sum(pcavar_shared)


varPart2PCA_shared <- fitExtractVarPartModel( pca_shared, form2, cv)
varPartSummary_shared = t(as.matrix(varPart2PCA_shared)) %*% (pcaprop_shared)

varPartSummary_shared

#varPartSummary_shared[-1]/sum(varPartSummary_shared[-1])*100 


#############

# Add day to the covariate file

cv$Day <- rep(NA, length(cv$Filename))

cv[grep("1R1", cv$Filename),"Day"] <- 1
cv[grep("1R2", cv$Filename),"Day"] <- 2
cv[grep("2R1", cv$Filename),"Day"] <- 3
cv[grep("2R2", cv$Filename),"Day"] <- 4

cv$Day <- factor(cv$Day)

# Define formula
form3 <- ~ (1|Individual) + (1|CellType) + (1|Plate.ID) + (1|Treatment.Name) + (1|Treatment.Name:CellType) + (1|Day)

varPart3PCA <- fitExtractVarPartModel( pca, form3, cv)
varPartSummary3 = t(as.matrix(varPart3PCA)) %*% (pcaprop)

varPartSummary3

# Now do only expressed genes

varPart3PCA_shared <- fitExtractVarPartModel( pca_shared, form3, cv)
varPartSummary3_shared = t(as.matrix(varPart3PCA_shared)) %*% (pcaprop_shared)

varPartSummary3_shared

##############

# Add Individual:CellType, Individual:CellType:Treatment

form4 <- ~ (1|Individual) + (1|CellType) + (1|Plate.ID) + (1|Treatment.Name) + (1|Treatment.Name:CellType) + (1|Day) + (1|Individual:CellType) + (1|Individual:CellType:Treatment.Name)

varPart4PCA_shared <- fitExtractVarPartModel( pca_shared, form4, cv)
varPartSummary4_shared = t(as.matrix(varPart4PCA_shared)) %*% (pcaprop_shared)

varPartSummary4_shared

##############

# Add Individual: Individual:Treatment.Name

form5 <- ~ (1|Individual) + (1|CellType) + (1|Plate.ID) + (1|Treatment.Name) + (1|Treatment.Name:CellType) + (1|Day) + (1|Individual:CellType) + (1|Individual:CellType:Treatment.Name) + (1|Individual:Treatment.Name)

varPart5PCA_shared <- fitExtractVarPartModel( pca_shared, form5, cv)
# Roger might want to see the plot:
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/VarPlot5_corrected_PCA.violin.pdf")
plotVarPart( varPart5PCA_shared )
dev.off()


varPartSummary5_shared = t(as.matrix(varPart5PCA_shared)) %*% (pcaprop_shared)

varPartSummary5_shared

# Remove everything not including individual from formula5 to calculate genetic effect. 
(varPartSummary5_shared[-c(1,2,7:10)]/sum(varPartSummary5_shared[-c(1,2,7:10)]))*100 

# Make a bar plot

Ind_df <- data.frame(Group = c("Individual", "Individual:CellType", "Individual:CellType:Treatment", "Individual:Treatment"), value = varPartSummary5_shared[-c(1,2,7:10)]/sum(varPartSummary5_shared[-c(1,2,7:10)])*100)

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/IndVariance_corrected.box.pdf")
ggplot(data=Ind_df, aes(x=Group, y=value)) +
  geom_bar(aes(fill=Group), stat="identity", color="black") + scale_fill_manual(values = c("#beaed4", "#7fc97f", "#fdc086", "#fdc086")) +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "Percentage of Variance Explained")
dev.off()


Ind_df
                          Group       value
1                    Individual 17.58965339
2           Individual:CellType 81.99947209
3 Individual:CellType:Treatment  0.40536549
4          Individual:Treatment  0.00550903




form6 <- ~ (1|Individual) + (1|Individual:CellType) + (1|Individual:CellType:Treatment.Name) + (1|Individual:Treatment.Name)

varPart6PCA_shared <- fitExtractVarPartModel( pca_shared, form6, cv)

# Roger might want to see the plot:
#pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/Variance_Partition/R2_R3/VarPlot6_corrected_PCA.violin.pdf")
#plotVarPart( varPart6PCA_shared )
#dev.off()


varPartSummary6_shared = t(as.matrix(varPart6PCA_shared)) %*% (pcaprop_shared)

varPartSummary6_shared

# Remove everything not including individual from formula5 to calculate genetic effect. 
(varPartSummary5_shared[-c(1,2,7:10)]/sum(varPartSummary5_shared[-c(1,2,7:10)]))*100 
