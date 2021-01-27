# Anthony Findley
# 11/6/2020

# Purpose: Compute correlations of Quasar beta values between cell types, replicates, treatments, individual, etc. Use pseudocount and drop DCM1R2 H2O. Only include SNPs with ASE from ANOVA model. THIS DOESN'T INCLUDE SNPS ON THE X CHROMOSOME!
# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/Quasar_output

library(tidyverse)
library(Hmisc)
library(data.table)

# Get SNPs with ASE by ANOVA
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/fixed_effect/interaction_noReadCov_DF_ANOVA_add1_noX/anova_lm_vp.Rd")

ASE <- read.table("all_allOutput_noHeaders.sorted.bed", header=F, stringsAsFactors=F)
colnames(ASE) <- c("chr", "pos0", "pos", "ref", "alt", "rsID", "af", "cell.line", "treatment", "ref.reads", "alt.reads", "beta", "beta.se", "pval", "qval")

# Add 1 to every ref and alt read and calculate log ratio
ASE$ref.reads1 <- ASE$ref.reads + 1
ASE$alt.reads1 <- ASE$alt.reads + 1
ASE$beta1 <- log(ASE$ref.reads1 / ASE$alt.reads1)

# Select the SNPs with ASE by ANOVA
ANOVA_ASE <- myanova_all_5df_finite %>% filter(padj < 0.1) %>% pull(SNP_Individual)

ASE$SNP_Ind <- paste0(ASE$rsID, "_", ASE$cell.line)
ASE <- ASE %>% filter(SNP_Ind %in% ANOVA_ASE)

# First create a unique identifier for each library
# Add info from covariate file to ASE file
cv <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/all_covar.txt", header=T, stringsAsFactors=F)

ASE_cv <- merge(ASE, cv, by.x = "treatment", by.y = "Filename")
ASE_cv <- ASE_cv[!(ASE_cv$Plate.ID == "DCM1R2" & ASE_cv$Treatment.Name == "Water"),]
ASE_cv[(ASE_cv$Plate.ID == "DCM1R2" & ASE_cv$Control.ID == "CO1"), "Control.ID"] <- "CO2"
ASE_cv$Plate.ID2 <- gsub("R.*", "", ASE_cv$Plate.ID)
ASE_cv$Plate.Rep <- gsub(".*R", "", ASE_cv$Plate.ID)

ASE_cv$libraryID <- paste(ASE_cv$Individual, ASE_cv$CellType, ASE_cv$Treatment.ID, ASE_cv$Treatment.Name, ASE_cv$Plate.ID, ASE_cv$Plate.ID2, ASE_cv$Plate.Rep, sep="_")

ASE_wide <- spread(ASE_cv[,c("rsID", "libraryID", "beta1")], key = libraryID, value = beta1 ) 

rownames(ASE_wide) <- ASE_wide$rsID
ASE_wide <- ASE_wide[,2:495]

ASE_wide_cor <- rcorr(as.matrix(ASE_wide), type="spearman")

# ++++++++++++++++++++++++++++
# flattenCorrMatrix    # I took this code from:http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat, nmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut],
	n = nmat[ut] # I added this myself
    )
}

flat_ASE_wide_cor <- flattenCorrMatrix(ASE_wide_cor$r, ASE_wide_cor$P, ASE_wide_cor$n)

# Now "row" is one library, "column" is the other. All pairs are unique:
dim(flat_ASE_wide_cor)
[1] 121771      5

flat_ASE_wide_cor$row <- as.character(flat_ASE_wide_cor$row)
flat_ASE_wide_cor$column <- as.character(flat_ASE_wide_cor$column)

# Now I'll parse the library names to get treatment
flat_ASE_wide_cor %>% separate(row, c("ind.x", "cell.x", "treatID.x", "treat.x", "plate.x", "cell_plate.x", "rep.x")) %>% separate(column, c("ind.y", "cell.y", "treatID.y", "treat.y", "plate.y", "cell_plate.y", "rep.y")) -> master_cor

# Make plot of p-value distribution and number of SNPs for master_cor
master_cor <- master_cor[order(master_cor$p), ]

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/Quasar_output/ASE_correlations/pval_hist.ANOVA_ASE.pseudo.noX.pdf")
p <- ggplot(master_cor, aes(x=p))
p + geom_histogram(bins=20, color="black") + labs(title = "P-value distribution of ASE correlation pairs (SNPs with ANOVA ASE)", x = "p-value") + theme_bw() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(trans = 'log10')
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/Quasar_output/ASE_correlations/numSNPs_hist.ANOVA_ASE.pseudo.noX.pdf")
p <- ggplot(master_cor, aes(x=n))
p + geom_histogram(bins=30, color="black") + labs(title = "Number of SNPs shared between libraries for correlation (SNPs with ANOVA ASE)", x = "Number of SNPs") + theme_bw() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) #+ scale_y_continuous(trans = 'log10')
dev.off()


#########
#########
# The variables we want to look at are individual, cell type, treatment, and plate/replicate. I want every possible combination of whether these variables are shared between libraries (2^4), but it doesn't make sense to include replicate info across cell types, so there will be fewer than 16 groups.

# Start with correlation across replicates (cell, individual, and treatment are held constant)
master_cor %>% filter(ind.x == ind.y & cell.x == cell.y & treat.x == treat.y  & plate.x != plate.y) -> sameIndCellTreat_difRep

# pair with the same individual, cell-type, and replicate but different condition; i.e. on same plate within individual but different treatments
master_cor %>% filter(ind.x == ind.y & cell.x == cell.y & treat.x != treat.y  & plate.x == plate.y) -> sameIndCellRep_diffTreat

# same individual and treatment but different cell type (don't worry about replicate)
master_cor %>% filter(ind.x == ind.y & treat.x == treat.y & cell.x != cell.y) -> sameIndTreat_diffCell

# different individual, same treatment, cell type, plate
master_cor %>% filter(ind.x != ind.y & treat.x == treat.y & cell.x == cell.y & plate.x == plate.y) -> sameTreatCellPlate_diffInd

# Different individual and cell type, same treatment
master_cor %>% filter(ind.x != ind.y & treat.x == treat.y & cell.x != cell.y) -> sameTreat_diffIndCell

# Now do each variable by itself
master_cor %>% filter(treat.x == treat.y) -> sameTreat
master_cor %>% filter(cell.x == cell.y) -> sameCell
master_cor %>% filter(plate.x == plate.y) -> samePlate
master_cor %>% filter(ind.x == ind.y) -> sameInd

# There are a lot more combinations, but I want to start making some plots.

# First, broad overview of 4 variables
master_cor$group <- "All"
sameTreat$group <- "Same treatment"
sameCell$group <- "Same cell type"
samePlate$group <- "Same plate"
sameInd$group <- "Same individual"

for_plot <- rbind(master_cor, sameTreat, sameCell, samePlate, sameInd)

g <- ggplot(for_plot, aes(x=group, y=cor, fill=group), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE correlation within variables (SNPs with ANOVA ASE)") + ylab("Spearman correlation") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) #+ scale_fill_manual(values=c("#fde0dd", "#fa9fb5", "#c51b8a"))


pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/Quasar_output/ASE_correlations/big_picture.ANOVA_ASE.pseudo.box.noX.pdf")
g
dev.off()

# Make another plot with the other variables I made. (I don't know what Roger and Francesca will want to see.)

sameIndCellTreat_difRep$group <- "Across rep"
sameIndCellRep_diffTreat$group <- "Across treat"
sameIndTreat_diffCell$group <- "Across cell"
sameTreatCellPlate_diffInd$group <- "Across Ind"

for_plot <- rbind(master_cor, sameIndCellTreat_difRep, sameIndCellRep_diffTreat, sameIndTreat_diffCell, sameTreatCellPlate_diffInd)
for_plot$group <- factor(for_plot$group)
for_plot$group <- relevel(for_plot$group, ref="All")

g <- ggplot(for_plot, aes(x=group, y=cor, fill=group), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE correlation across variables (SNPs with ANOVA ASE)") + ylab("Spearman correlation") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))


pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/Quasar_output/ASE_correlations/across_big_picture.ANOVA_ASE.pseudo.box.noX.pdf")
g
dev.off()


# For across treatment ASE correlations, separate control vs control, treatment vs control, and treatment vs treatment

sameIndCellRep_diffTreat %>% filter(treat.x == "Water" & treat.y == "Ethanol") -> acrossTreat_ConCon
sameIndCellRep_diffTreat %>% filter((grepl("CO", sameIndCellRep_diffTreat$treatID.x) & !grepl("CO", sameIndCellRep_diffTreat$treatID.y)) | (!grepl("CO", sameIndCellRep_diffTreat$treatID.x) & grepl("CO", sameIndCellRep_diffTreat$treatID.y)) ) -> acrossTreat_TreatCon
sameIndCellRep_diffTreat %>% filter(!grepl("CO", sameIndCellRep_diffTreat$treatID.x) & !grepl("CO", sameIndCellRep_diffTreat$treatID.y)) -> acrossTreat_TreatTreat

acrossTreat_ConCon$group <- "Control_vs_Control"
acrossTreat_TreatCon$group <- "Treat_vs_Control"
acrossTreat_TreatTreat$group <- "Treat_vs_Treat"

for_plot <- rbind(master_cor, sameIndCellTreat_difRep, sameIndCellRep_diffTreat, acrossTreat_ConCon, acrossTreat_TreatCon, acrossTreat_TreatTreat)
for_plot$group <- factor(for_plot$group)
for_plot$group <- relevel(for_plot$group, ref="All")

g <- ggplot(for_plot, aes(x=group, y=cor, fill=group), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE correlation across variables (SNPs with ANOVA ASE)") + ylab("Spearman correlation") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))


pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/Quasar_output/ASE_correlations/across_big_picture.ConTreat.ANOVA_ASE.pseudo.box.noX.pdf")
g
dev.off()


# Roger and Francesca suggested a plot with “across cell-types”, “across individuals*(unphased)“, “across replicates”, “across treatments”, “control/control”.
for_paper <- rbind(sameIndTreat_diffCell, sameTreatCellPlate_diffInd, sameIndCellTreat_difRep, sameIndCellRep_diffTreat, acrossTreat_ConCon)

g <- ggplot(for_paper, aes(x=group, y=cor, fill=group), color="black") + geom_boxplot(notch=TRUE, outlier.shape=NA) + theme_bw() + labs(title="ASE correlation across variables (SNPs with ANOVA ASE)") + ylab("Spearman correlation") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/Quasar_output/ASE_correlations/across_big_picture.forPaper.ANOVA_ASE.pseudo.box.noX.pdf")
g
dev.off()
