# Anthony Findley
# 8/13/2020

# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/fixed_effect
# Purpose: Having 0 reads in either the reference or alternate messes up Quasar's ASE estimates. So add 1 to ref and alt for every ASE entry prior to running linear model. Now I'm changing DCM1R2 to use all ethanol controls because we thought there was a problem with water. This is the clean version and I remove sex chromosomes.

# https://r4stats.com/2017/04/18/group-by-modeling-in-r-made-easy/
# https://cran.r-project.org/web/packages/broom/vignettes/broom_and_dplyr.html

library(Hmisc)
library(data.table)
library(tidyverse)
library(broom)
library(reshape) # Not sure if I need this or not
library(VennDiagram)
library(RColorBrewer)
library(corrplot)
library(pheatmap)
library(UpSetR)
library(lme4)

# To re-create whole environment without needing to rerun linear models, load("interaction_noReadCov_DF_ANOVA_add1_noX/anova_lm_vp.Rd")
# Loads: "myanova_all_5df_finite", "mylm_all_5df", "vp_all"

# Function to implement linear model and include DF
mylm <- function(data){
	SNP_Individual=paste0(data$rsID[1],"_",data$Individual[1])
	obj <- lm(beta1 ~ Control.ID + CellType*Tr, data = data)
	mysum <- summary(obj)
	mycoeff <- as.data.frame(mysum$coefficients)  
	colnames(mycoeff) <- c("estimate","std.error","statistic","p.value")
	mycoeff	%>% rownames_to_column("term") %>%
		mutate(deg.f=obj$df.residual,
			rank=obj$qr$rank, 
			SNP_Individual=SNP_Individual)#,
			#kappa=kappa.qr(obj$qr))
}

myanova_reduced <- function(data){
	SNP_Individual=paste0(data$rsID[1],"_",data$Individual[1])
	full = lm(beta1 ~ Control.ID + CellType*Tr, data=data)
	reduced = lm(beta1 ~ 0, data=data)
	anov <- data.frame(anova(reduced, full))
	colnames(anov) <- c("Res.Df", "RSS", "Df", "Sum.Sq", "F.value", "p.value")
	anov %>% mutate(SNP_Individual=SNP_Individual)
}



ASE <- fread("../../Quasar_output/all_allOutput_noHeaders.sorted.bed", header=F, data.table=F)
colnames(ASE) <- c("chr", "pos0", "pos", "ref", "alt", "rsID", "af", "cell.line", "treatment", "ref.reads", "alt.reads", "beta", "beta.se", "pval", "qval")

# Remove sex chromsomes and add 1 to every ref and alt read and calculate log ratio
ASE %>% filter(chr != "X" & chr != "Y") %>% mutate(ref.reads1 = ref.reads + 1, alt.reads1 = alt.reads + 1, beta1 = log(ref.reads1 / alt.reads1)) -> ASE

cv <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/all_covar.txt", header=T, stringsAsFactors=F)

ASE <- merge(ASE, cv, by.x = "treatment", by.y = "Filename")

# For DCM1R2, change water control to ethanol
ASE <- ASE[!(ASE$Plate.ID == "DCM1R2" & ASE$Treatment.Name == "Water"),]
ASE[(ASE$Plate.ID == "DCM1R2" & ASE$Control.ID == "CO1"), "Control.ID"] <- "CO2"

ASE$SNP_Individual <- factor(paste0(ASE$rsID, "_", ASE$Individual))


# I think I can only run the above for SNP-individual pairs with multiple treatments and cell types
ASE.small <- ASE[,c("SNP_Individual", "CellType", "Treatment.Name")]

ASE.small.celltype <- dcast(ASE.small, SNP_Individual ~ CellType) # Creates df for every SNP_Individual, number of CM, IPSC, and IPSC treatments

# Only keep these SNPs to test
testable_SNP_Ind <- as.character(ASE.small.celltype[which(ASE.small.celltype$CM > 1 & ASE.small.celltype$IPSC > 1 & ASE.small.celltype$IPSC > 1), "SNP_Individual"])

ASE$SNP_Individual <- as.character(ASE$SNP_Individual) # I don't think we want these as factors yet

ASE_testable <- ASE[ASE$SNP_Individual %in% testable_SNP_Ind,]
ASE_testable$SNP_Individual <- factor(ASE_testable$SNP_Individual)

by_SNP.Ind <- group_by(ASE_testable, SNP_Individual)

by_SNP.Ind$Treatment.Name <- factor(by_SNP.Ind$Treatment.Name)
by_SNP.Ind$CellType <- factor(by_SNP.Ind$CellType)


############################

by_SNP.Ind$Tr <- as.character(by_SNP.Ind$Treatment.Name)
by_SNP.Ind$Tr[by_SNP.Ind$Tr=="Ethanol"] = "Control"
by_SNP.Ind$Tr[by_SNP.Ind$Tr=="Water"] = "Control"
by_SNP.Ind$Tr = factor(by_SNP.Ind$Tr)
by_SNP.Ind$Tr = relevel(by_SNP.Ind$Tr,ref="Control")
str(by_SNP.Ind$Tr)
 Factor w/ 13 levels "Control","Acetaminophin",..: 9 9 9 9 9 9 9 9 9 9 ...

str(by_SNP.Ind$CellType)
 Factor w/ 3 levels "CM","IPSC","LCL": 1 1 1 1 1 1 1 1 1 1 ...

by_SNP.Ind$CellType = relevel(by_SNP.Ind$CellType,ref="LCL")


by_SNP.Ind.small <- dcast(by_SNP.Ind, SNP_Individual ~ CellType) # Creates df for every SNP_Individual, number of CM, IPSC, and LCL treatments

testable_SNP_Ind <- as.character(by_SNP.Ind.small[which(by_SNP.Ind.small$CM > 4 & by_SNP.Ind.small$IPSC > 4 & by_SNP.Ind.small$LCL > 4), "SNP_Individual"])

by_SNP.Ind_testable <- by_SNP.Ind[by_SNP.Ind$SNP_Individual %in% testable_SNP_Ind,]

# Run the linear model
aux <- by_SNP.Ind_testable %>% nest()  
mylist <- aux$data
names(mylist) <- aux$SNP_Individual
mylm_all <- map_dfr(mylist,mylm)

# Run the ANOVA
myanova_all <- map_dfr(mylist,myanova_reduced)


# Remove values with NA for p.value
mylm_all <- mylm_all[which(!is.na(mylm_all$p.value)),]

# I think I want to remove SNPs with fewer than 5 degrees of freedom.
mylm_all_5df <- mylm_all[which(mylm_all$deg.f > 4),]
dim(mylm_all_5df)
[1] 2340107       8

# Only keep ANOVA SNPs with 5 degrees of freedom from linear model.
myanova_all_5df <- myanova_all[ myanova_all$SNP_Individual %in% unique(mylm_all_5df$SNP_Individual), ]

# There are only 2 lines per SNP, one for the reduced and one for the full. The reduced p.value is always NA. Get rid of all lines with p.value = NA
myanova_all_5df_finite <- myanova_all_5df[ is.finite(myanova_all_5df$p.value),]

# Now make plot
myanova_all_5df_finite <- myanova_all_5df_finite[order(myanova_all_5df_finite$p.value),]
myanova_all_5df_finite$padj <- p.adjust(myanova_all_5df_finite$p.value, method="BH")
sum(myanova_all_5df_finite$padj < 0.1) # How much ASE with FDR < 0.1
[1] 15497  # There was 15952 when including X chromosome

# How many genes?
gene_anno <- read.table("../../Quasar_output/ASE_SNPs.genes.uniq.txt", header=F, stringsAsFactors=F)
colnames(gene_anno) <- c("chr", "pos0", "pos1", "ref", "alt", "rsID", "ensg", "gene_type", "g.id")

myanova_all_5df_finite$rsID <- sapply(strsplit(myanova_all_5df_finite$SNP_Individual, "_"), "[", 1)
myanova_all_5df_finite_gene <- merge(myanova_all_5df_finite, gene_anno, by="rsID")
length(unique(myanova_all_5df_finite_gene$g.id)) # All genes tested
[1] 10142   # Was 10267 with X chromosome
myanova_all_5df_finite_gene %>% filter(padj < 0.1) %>% pull(g.id) %>% unique() %>% length()  # 5640; number of genes with significant ASE

myanova_all_5df_finite$exp <- ppoints(length(myanova_all_5df_finite$SNP_Individual))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/ANOVA_reduced.noDCM1R2_H2O.5df.noReadCov.noX.qq.pdf")
p <- ggplot(myanova_all_5df_finite, aes(-log10(exp), -log10(p.value)))
print(p + geom_point(size=2) + labs(title = "ANOVA reduced vs full model for all cell types together", x = "-log10(exp.pvalue)", y = "-log10(obs.pvalue)") + geom_abline(intercept=0, slope=1) + theme_bw() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")))
dev.off()

#######################################
#######################################

# Only keep SNPs with ANOVA ASE and then consider cell type, treatment, and cell:treatment together
all_anova_ASE <- unique(myanova_all_5df_finite[myanova_all_5df_finite$padj < 0.1, "SNP_Individual"])

mylm_all_5df_anovaASE <- mylm_all_5df[which(mylm_all_5df$SNP_Individual %in% all_anova_ASE),] # All ASE

mylm_all_5df_anovaASE_interact <- mylm_all_5df_anovaASE[grepl("T", mylm_all_5df_anovaASE$term),] # Picks up Treatment, CellType, and interaction

# Now split data into cell type, treatment, and cell type:treatment, and recalculate padj
mylm_all_5df_anovaASE_interact_cell <- mylm_all_5df_anovaASE_interact[which(grepl("^CellType", mylm_all_5df_anovaASE_interact$term) & !grepl("Tr", mylm_all_5df_anovaASE_interact$term)),]
mylm_all_5df_anovaASE_interact_cell <- mylm_all_5df_anovaASE_interact_cell[order(mylm_all_5df_anovaASE_interact_cell$p.value),]
mylm_all_5df_anovaASE_interact_cell$exp <- ppoints(length(mylm_all_5df_anovaASE_interact_cell$term))
mylm_all_5df_anovaASE_interact_cell$group <- rep("Cell Type Effect", length(mylm_all_5df_anovaASE_interact_cell$exp))
mylm_all_5df_anovaASE_interact_cell$padj <- p.adjust(mylm_all_5df_anovaASE_interact_cell$p.value, method="BH")
sum(mylm_all_5df_anovaASE_interact_cell$padj < 0.1)
[1] 9152 # Total cell type cASE; was 9743 with X
length(unique(mylm_all_5df_anovaASE_interact_cell[which(mylm_all_5df_anovaASE_interact_cell$padj < 0.1), "SNP_Individual"]))
[1] 6391 # SNP_Individuals with cell type cASE; was 6733 with X

mylm_all_5df_anovaASE_interact_Tr <- mylm_all_5df_anovaASE_interact[grep("^Tr", mylm_all_5df_anovaASE_interact$term),]
mylm_all_5df_anovaASE_interact_Tr <- mylm_all_5df_anovaASE_interact_Tr[order(mylm_all_5df_anovaASE_interact_Tr$p.value),]
mylm_all_5df_anovaASE_interact_Tr$exp <- ppoints(length(mylm_all_5df_anovaASE_interact_Tr$term))
mylm_all_5df_anovaASE_interact_Tr$group <- rep("Treatment Effect", length(mylm_all_5df_anovaASE_interact_Tr$exp))
mylm_all_5df_anovaASE_interact_Tr$padj <- p.adjust(mylm_all_5df_anovaASE_interact_Tr$p.value, method="BH")
sum(mylm_all_5df_anovaASE_interact_Tr$padj < 0.1)
[1] 2046 # Total Tr cASE; was 2076 with X chrom
length(unique(mylm_all_5df_anovaASE_interact_Tr[which(mylm_all_5df_anovaASE_interact_Tr$padj < 0.1), "SNP_Individual"]))
[1] 1482 # SNP_Individuals with Tr cASE; was 1504 with X

mylm_all_5df_anovaASE_interact_TrCell <- mylm_all_5df_anovaASE_interact[grepl(":", mylm_all_5df_anovaASE_interact$term),]
mylm_all_5df_anovaASE_interact_TrCell <- mylm_all_5df_anovaASE_interact_TrCell[order(mylm_all_5df_anovaASE_interact_TrCell$p.value),]
mylm_all_5df_anovaASE_interact_TrCell$exp <- ppoints(length(mylm_all_5df_anovaASE_interact_TrCell$term))
mylm_all_5df_anovaASE_interact_TrCell$group <- rep("Tr:Cell Effect", length(mylm_all_5df_anovaASE_interact_TrCell$exp))
mylm_all_5df_anovaASE_interact_TrCell$padj <- p.adjust(mylm_all_5df_anovaASE_interact_TrCell$p.value, method="BH")
sum(mylm_all_5df_anovaASE_interact_TrCell$padj < 0.1)
[1] 1457 # Total Tr:Cell cASE; was 1468 with X
length(unique(mylm_all_5df_anovaASE_interact_TrCell[which(mylm_all_5df_anovaASE_interact_TrCell$padj < 0.1), "SNP_Individual"]))
[1] 1011 # SNP_Individuals with Tr:Cell cASE; was 1020 with X

mylm_all_5df_anovaASE_interact_group_forPlot <- rbind(mylm_all_5df_anovaASE_interact_cell, mylm_all_5df_anovaASE_interact_Tr, mylm_all_5df_anovaASE_interact_TrCell)

# Plot control ID for ASE SNPs (should not be much significant)
mylm_all_5df_anovaASE_control <- mylm_all_5df_anovaASE[grepl("Control", mylm_all_5df_anovaASE$term),] # Picks up Treatment, CellType, and interaction
dim(mylm_all_5df_anovaASE_control)

mylm_all_5df_anovaASE_control$padj <- p.adjust(mylm_all_5df_anovaASE_control$p.value, method="BH")
sum(mylm_all_5df_anovaASE_control$padj < 0.1)
[1] 0  # Number of control ASE; was 0

# Set up qq plot
mylm_all_5df_anovaASE_control <- mylm_all_5df_anovaASE_control[order(mylm_all_5df_anovaASE_control$p.value),]
mylm_all_5df_anovaASE_control$exp <- ppoints(length(mylm_all_5df_anovaASE_control$term))

# Add control.ID to plot with treatment, cell type, and interaction
mylm_all_5df_anovaASE_control$group <- rep("Control ID", length(mylm_all_5df_anovaASE_control$p.value))
combined_for_plot <- rbind(mylm_all_5df_anovaASE_interact_group_forPlot, mylm_all_5df_anovaASE_control)

combined_for_plot$group <- factor(combined_for_plot$group, levels=unique(combined_for_plot$group))

p <- ggplot(combined_for_plot, aes(-log10(exp), -log10(p.value), color=group)) + geom_point(size=2) + labs(title = "Interaction Model by Group using ANOVA ASE", x = "-log10(exp.pvalue)", y = "-log10(obs.pvalue)") + geom_abline(intercept=0, slope=1) + theme_bw() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + scale_color_manual(values=c("#1b9e77", "#d95f02", "#7570b3", "#060606"))

p
ggsave("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/allCells_Interact_Group.controlID.no_DCM1R2_H2O.5df.noReadCov.noX.qq.png")

########################
########################
# I want to make a heatmap of counts for cASE between each pair of treatments, but I'm going to start with an upsetR plot because it is easier
mylm_all_5df_anovaASE_interact_Tr_sig <- mylm_all_5df_anovaASE_interact_Tr[ mylm_all_5df_anovaASE_interact_Tr$padj < 0.1, ]
cASE_upset <- cast(mylm_all_5df_anovaASE_interact_Tr_sig, SNP_Individual ~ term)
cASE_upset[] <- lapply(cASE_upset, as.character)
rownames(cASE_upset) <- cASE_upset$SNP_Individual
cASE_upset$SNP_Individual <- NULL
cASE_upset[!is.na(cASE_upset)] <- 1
cASE_upset[is.na(cASE_upset)] <- 0

cASE_upset[] <- lapply(cASE_upset, as.numeric)

# Make heatmap of significant cASE sharing including all singificant SNPs
shared_cASE <- t(as.matrix(cASE_upset)) %*% as.matrix(cASE_upset)
shared_cASE_melt <- data.frame(melt(shared_cASE))

# Now select only SNPs which have cASE in at least 2 conditions.
cASE_upset_double <- data.frame(cASE_upset[rowSums(cASE_upset) > 1,])

shared_cASE_double <- t(as.matrix(cASE_upset_double)) %*% as.matrix(cASE_upset_double)
shared_cASE_double_melt <- data.frame(melt(shared_cASE_double))

### Rather than plotting the number of SNPs which overlap between treatments, Roger and Francesca asked me to plot (# shared SNPs)/(# SNPs in union between conditions). I need to find the union between all pairs of treatments.
# I'm going to create 2 identical DFs with different column names which will contain the total number of cASE for a treatment. Then I'll merge these with "shared_cASE_melt"

X1_df <- data.frame(X1 = names(colSums(cASE_upset)), X1_total = colSums(cASE_upset))
X2_df <- data.frame(X2 = names(colSums(cASE_upset)), X2_total = colSums(cASE_upset))

shared_cASE_melt_union <- merge(shared_cASE_melt, X1_df)
shared_cASE_melt_union <- merge(shared_cASE_melt_union, X2_df)

# Calculate union of cASE SNPs. "value" is the overlap between the two
shared_cASE_melt_union %>% mutate(uni = X1_total + X2_total - value, prop = value/uni) -> shared_cASE_melt_union

# Put into wide format so I can plot with pheatmap
shared_cASE_melt_union %>% select(X1, X2, prop) %>% spread(key = X2, value = prop) %>% column_to_rownames(var = "X1") -> shared_cASE_union

# Make the diagonal 0 and try again.

shared_cASE_union2 <- shared_cASE_union
shared_cASE_union2[shared_cASE_union2 == 1] <- NA

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/allCells_cASE_shared.union_noDiag..noX.pheat.pdf")
pheatmap(shared_cASE_union2)
dev.off()


######################
######################
# Compare distribution of beta values to what we found before. EDIT: SKIP THIS.
#mylm_all_5df_anovaASE_interact_group_forPlot$ASE <- "ANOVA"
#mylm_all_5df_ASE_interact_group_forPlot$ASE <- "Linear_Model"

#beta_bothModel <- rbind(mylm_all_5df_anovaASE_interact_group_forPlot, mylm_all_5df_ASE_interact_group_forPlot)

#pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/allCells_Interact_Beta.5df.noReadCov.dens.pdf")
#p <- ggplot(beta_bothModel, aes(x=estimate))
#print(p + geom_density(aes(fill=ASE), color="black", alpha=0.5) + facet_grid(group ~ .) + labs(title = "Beta distribution of coefficients for SNPs with ASE by linear or ANOVA model", x = "beta") + theme_bw() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + xlim(-5,5))
#dev.off()

####################
####################
# Plot some of the top SNPs
forest_plot <- function(snpID){

	for_plot <- ASE[which(ASE$SNP_Individual == snpID),]

    for_plot$CI_hi <- for_plot$beta1 + (1.96 * for_plot$beta.se)
    for_plot$CI_lo <- for_plot$beta1 - (1.96 * for_plot$beta.se)
    for_plot$Treatment.ID <- factor(for_plot$Treatment.ID)

	for_plot <- for_plot[order(for_plot$Treatment.Name, for_plot$beta1),]
    for_plot$id <- paste0(for_plot$Plate.ID, "_", for_plot$Treatment.Name)
    for_plot$id <- factor(for_plot$id, levels=for_plot$id)

    fp <- ggplot(data=for_plot, aes(x=id, y=beta1, ymin=CI_lo, ymax=CI_hi)) +
    geom_pointrange(aes(color=Treatment.Name, pch=CellType)) +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Plate and Treatment") + ylab("ASE (95% CI)") + theme_bw() + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=Treatment.Name), width=0.5) + ggtitle(paste0(snpID, " in model with all cell types"))

    pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/forest/", snpID, ".fp.pdf"), height=7)
        print(fp)
    dev.off()
		
	# Now plot the beta's from the Tr cASE model
	beta_plot <- mylm_all_5df_anovaASE_interact_Tr[which(mylm_all_5df_anovaASE_interact_Tr$SNP_Individual == snpID),]
    beta_plot$CI_hi <- beta_plot$estimate + (1.96 * beta_plot$std.error)
    beta_plot$CI_lo <- beta_plot$estimate - (1.96 * beta_plot$std.error)

	beta_plot <- beta_plot[order(beta_plot$estimate),]
    beta_plot$term <- factor(beta_plot$term, levels=beta_plot$term)

    fp <- ggplot(data=beta_plot, aes(x=term, y=estimate, ymin=CI_lo, ymax=CI_hi)) +
        geom_pointrange(aes(color=padj < 0.1)) +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Plate and Treatment") + ylab("ASE (95% CI)") + theme_bw() + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=padj < 0.1), width=0.5) + ggtitle(paste0(snpID, " for Treatment cASE (model with all cell types)"))

    pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/forest/", snpID, ".Tr.model_betas.fp.pdf"), height=3.5)
        print(fp)
    dev.off()
	
	# Now plot the beta's from the Tr:Cell cASE model
	beta_plot <- mylm_all_5df_anovaASE_interact_TrCell[which(mylm_all_5df_anovaASE_interact_TrCell$SNP_Individual == snpID),]
    beta_plot$CI_hi <- beta_plot$estimate + (1.96 * beta_plot$std.error)
    beta_plot$CI_lo <- beta_plot$estimate - (1.96 * beta_plot$std.error)

	beta_plot <- beta_plot[order(beta_plot$estimate),]
    beta_plot$term <- factor(beta_plot$term, levels=beta_plot$term)

    fp <- ggplot(data=beta_plot, aes(x=term, y=estimate, ymin=CI_lo, ymax=CI_hi)) +
        geom_pointrange(aes(color=padj < 0.1)) +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Plate and Treatment") + ylab("ASE (95% CI)") + theme_bw() + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=padj < 0.1), width=0.5) + ggtitle(paste0(snpID, " for Treatment:Cell cASE (model with all cell types)"))

    pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/forest/", snpID, ".TrCell.model_betas.fp.pdf"), height=3.5)
        print(fp)
    dev.off()

	# Now plot the beta's from the Tr:Cell cASE model
	beta_plot <- mylm_all_5df_anovaASE_interact_cell[which(mylm_all_5df_anovaASE_interact_cell$SNP_Individual == snpID),]
    beta_plot$CI_hi <- beta_plot$estimate + (1.96 * beta_plot$std.error)
    beta_plot$CI_lo <- beta_plot$estimate - (1.96 * beta_plot$std.error)

	beta_plot <- beta_plot[order(beta_plot$estimate),]
    beta_plot$term <- factor(beta_plot$term, levels=beta_plot$term)

    fp <- ggplot(data=beta_plot, aes(x=term, y=estimate, ymin=CI_lo, ymax=CI_hi)) +
        geom_pointrange(aes(color=padj < 0.1)) +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Plate and Treatment") + ylab("ASE (95% CI)") + theme_bw() + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=padj < 0.1), width=0.5) + ggtitle(paste0(snpID, " for Cell cASE (model with all cell types)"))

    pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/forest/", snpID, ".Cell.model_betas.fp.pdf"), height=1.5)
        print(fp)
    dev.off()

}


###########################################
###########################################
# Only keep SNPs which have at least one water and one ethanol control for each cell type (except DCM1R2, which needs just one ethanol)

ASE_etoh <- ASE[ASE$Treatment.Name == "Ethanol",]
ASE_etoh_cast <- dcast(ASE_etoh, SNP_Individual ~ CellType) # Creates df of # of etoh per cell type
good_ASE_etoh <- ASE_etoh_cast[ASE_etoh_cast$CM > 0 & ASE_etoh_cast$IPSC > 0 & ASE_etoh_cast$LCL > 0, "SNP_Individual"]
length(good_ASE_etoh)
[1] 68599  # Was 69133 with X chrom


ASE_h2o <- ASE[ASE$Treatment.Name == "Water",]
ASE_h2o_cast <- dcast(ASE_h2o, SNP_Individual ~ CellType) # Creates df of # of h2o per cell type
good_ASE_h2o <- ASE_h2o_cast[ASE_h2o_cast$CM > 0 & ASE_h2o_cast$IPSC > 0 & ASE_h2o_cast$LCL > 0, "SNP_Individual"]
length(good_ASE_h2o)
[1] 66643 # was 67163 with X chrom

# How many SNPs can we test with at least 1 of each control:
mylm_all_5df_anovaASE_interact_Tr %>% filter(SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) %>% pull(SNP_Individual) %>% unique() %>% length() # 13009; 13400 with X

mylm_all_5df_anovaASE_interact_Tr %>% filter(padj < 0.1 & SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) -> Tr_cASE_etoh_h2o

# I need to make this again with another name because it changes below
mylm_all_5df_anovaASE_interact_Tr %>% filter(padj < 0.1 & SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) -> Tr_cASE_etoh_h2o_2

dim(Tr_cASE_etoh_h2o)
[1] 1409   11  # 1102 unique SNP_Ind; was 1438 and 1077 with X chrom

forest_plot("rs1046165_GM18858")
forest_plot("rs12063729_GM18855")
forest_plot("rs1046164_GM18858")
forest_plot("rs77878473_GM18912")

mylm_all_5df_anovaASE_interact_TrCell %>% filter(padj < 0.1 & SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) -> TrCell_cASE_etoh_h2o

dim(TrCell_cASE_etoh_h2o)
[1] 929 # 715 unique SNP_Ind; was 943 and 725 with X

# How many genes?
TrCell_cASE_etoh_h2o$rsID <- sapply(strsplit(TrCell_cASE_etoh_h2o$SNP_Individual, "_"), "[", 1)
TrCell_cASE_etoh_h2o_gene <- merge(TrCell_cASE_etoh_h2o, gene_anno, by="rsID")
length(unique(TrCell_cASE_etoh_h2o_gene$g.id))
[1] 689 # was 695 with X

mylm_all_5df_anovaASE_interact_cell %>% filter(padj < 0.1 & SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) -> Cell_cASE_etoh_h2o

dim(Cell_cASE_etoh_h2o)
[1] 7866 # 5452 unique SNP_Ind; was 8401 and 5760 with X

# How many genes? (Later on I add the gene annotation)
Cell_cASE_etoh_h2o$rsID <- sapply(strsplit(Cell_cASE_etoh_h2o$SNP_Individual, "_"), "[", 1)
Cell_cASE_etoh_h2o_gene <- merge(Cell_cASE_etoh_h2o, gene_anno, by="rsID")
length(unique(Cell_cASE_etoh_h2o_gene$g.id))
[1] 2822 # was 2943 with X

# How many total genes?
length(unique(c(unique(Tr_cASE_etoh_h2o$g.id), unique(TrCell_cASE_etoh_h2o_gene$g.id), unique(Cell_cASE_etoh_h2o_gene$g.id))))
[1] 3198  # 3314 with X

# How many SNPs with significant GxE?
length(unique(c(unique(Tr_cASE_etoh_h2o$SNP_Individual), unique(TrCell_cASE_etoh_h2o_gene$SNP_Individual), unique(Cell_cASE_etoh_h2o_gene$SNP_Individual)))) # 5984; 6274 with X



################################################
################################################
################################################
# Now go ahead and do the mixed model to calculate ASE variance

calcVarPart <- function(data){
	data$Tr <- factor(data$Tr)
	data$Control.ID <- factor(data$Control.ID)
	data$CellType <- factor(data$CellType)
	obj <- lmer(beta1 ~ 1 + (1|Control.ID) + (1|CellType) + (1|Tr) + (1|CellType:Tr), data=data, control=lmerControl(check.conv.hess="stop", check.conv.grad="stop"))
	conv <- obj@optinfo$conv$opt # I think "0" means that it converged
	aux <- as.data.frame(VarCorr(obj))
	vp <- c(aux$vcov, conv)
	names(vp) <- c(aux$grp, "convergence")
	vp
}

by_SNP.Ind_testable %>% nest() %>% mutate(vp=map(data, possibly(calcVarPart, otherwise = NULL))) -> aux 

# Now only keep the SNP_Individuals which have 5 df and are tested for ASE
aux_5df <- aux[as.character(aux$SNP_Individual) %in% mylm_all_5df[, "SNP_Individual"], ]

# Get rid of the NULL rows
aux_5df_noNull <- aux_5df[-which(lengths(aux_5df$vp) == 0),]

vp_all <- do.call(rbind,aux_5df_noNull$vp)
rownames(vp_all) <- aux_5df_noNull$SNP_Individual

# Select just the SNPs with ASE by ANOVA model which have at least one ethanol and water control per cell type
vp_all_ASE <- vp_all[rownames(vp_all) %in% all_anova_ASE & rownames(vp_all) %in% good_ASE_h2o & rownames(vp_all) %in% good_ASE_etoh, ] # 12948 SNPs; was 15,952 with X

# Drop convergence (all are 0)
vp_all_ASE <- vp_all_ASE[,c(1:5)]

# Now calculate percent variance explained
vp_percent <- vp_all_ASE/rowSums(vp_all_ASE)
vp_percent_melt <- melt(vp_percent, id.vars=0)
vp_percent_melt$X2 <- factor(vp_percent_melt$X2, levels=c("CellType", "Tr", "CellType:Tr", "Control.ID", "Residual"))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/varPart_allCells.ASE.noX.box.pdf")
ggplot(vp_percent_melt, aes(x=X2, y=value, fill=X2)) +  geom_boxplot() + theme_bw() + ylab("Percent Variance Explained")
dev.off()

# Remove outliers
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/varPart_allCells.ASE.noOut.noX.box.pdf")
ggplot(vp_percent_melt, aes(x=X2, y=value, fill=X2)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + theme_bw() + ylab("Percent Variance Explained")
dev.off()

### 
# See if I can add plate
calcVarPart_plate <- function(data){
	data$Tr <- factor(data$Tr)
	data$Control.ID <- factor(data$Control.ID)
	data$CellType <- factor(data$CellType)
	data$Plate.ID <- factor(data$Plate.ID)
	obj <- lmer(beta1 ~ 1 + (1|Control.ID) + (1|CellType) + (1|Tr) + (1|CellType:Tr) + (1|Plate.ID), data=data, control=lmerControl(check.conv.hess="stop", check.conv.grad="stop"))
	conv <- obj@optinfo$conv$opt # I think "0" means that it converged
	aux <- as.data.frame(VarCorr(obj))
	vp <- c(aux$vcov, conv)
	names(vp) <- c(aux$grp, "convergence")
	vp
}

by_SNP.Ind_testable %>% nest() %>% mutate(vp=map(data, possibly(calcVarPart_plate, otherwise = NULL))) -> aux_plate

# Now only keep the SNP_Individuals which have 5 df and are tested for ASE
aux_plate_5df <- aux_plate[as.character(aux_plate$SNP_Individual) %in% mylm_all_5df[, "SNP_Individual"], ]

# Get rid of the NULL rows
aux_plate_5df_noNull <- aux_plate_5df[-which(lengths(aux_plate_5df$vp) == 0),]

vp_plate_all <- do.call(rbind,aux_plate_5df_noNull$vp)
rownames(vp_plate_all) <- aux_plate_5df_noNull$SNP_Individual

# Select just the SNPs with ASE by ANOVA model which have at least one ethanol and water control per cell type
vp_plate_all_ASE <- vp_plate_all[rownames(vp_plate_all) %in% all_anova_ASE & rownames(vp_plate_all) %in% good_ASE_h2o & rownames(vp_plate_all) %in% good_ASE_etoh, ] # 12,983 SNPs; was 12,948 without plate

# Drop convergence (all are 0)
vp_plate_all_ASE <- vp_plate_all_ASE[,c(1:6)]

# Now calculate percent variance vp_plate_all_ASE
vp_plate_percent <- vp_plate_all_ASE/rowSums(vp_plate_all_ASE)
vp_plate_percent_melt <- melt(vp_plate_percent, id.vars=0)
vp_plate_percent_melt$X2 <- factor(vp_plate_percent_melt$X2, levels=c("CellType", "Tr", "CellType:Tr", "Plate.ID", "Control.ID", "Residual"))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/varPart_allCells_plate.ASE.noX.box.pdf")
ggplot(vp_plate_percent_melt, aes(x=X2, y=value, fill=X2)) +  geom_boxplot() + theme_bw() + ylab("Percent Variance Explained")
dev.off()

# Remove outliers
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/varPart_allCells_plate.ASE.noOut.noX.box.pdf")
ggplot(vp_plate_percent_melt, aes(x=X2, y=value, fill=X2)) +  geom_boxplot(outlier.shape=NA, notch=TRUE) + theme_bw() + ylab("Percent Variance Explained")
dev.off()

# Make it a violin plot 
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/varPart_allCells_plate.ASE.noX.violin.pdf")
ggplot(vp_plate_percent_melt, aes(x=X2, y=value)) +  geom_violin(aes(fill=X2), scale="width") + geom_boxplot(width=0.1, fill="grey") + theme_bw() + theme(legend.position = "none") + ylab("Percent Variance Explained")
dev.off()

#####################
#####################
# Make barplot for number of cASE
# Get cASE SNPs per treatment
Tr_cASE_etoh_h2o$Treatment <- gsub("^Tr", "", Tr_cASE_etoh_h2o$term)
cASEbyTreat <- aggregate(cbind(count = SNP_Individual) ~ Treatment, data=unique(Tr_cASE_etoh_h2o[,c("SNP_Individual", "Treatment")]), FUN = function(x){NROW(x)})

cASEbyTreat <- cASEbyTreat[order(cASEbyTreat$count),]
cASEbyTreat$Treatment <- factor(cASEbyTreat$Treatment, levels=cASEbyTreat$Treatment) # We have slighlty more cASE per treatment without X chromosome

# Now do one without the colors
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/allCells_Interact_cASEbyTreat.noColor.noX.bar.pdf")
ggplot(data=cASEbyTreat, aes(x=Treatment, y=count)) +
  geom_bar(stat="identity", color="black", fill="palevioletred1") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "Number of cASE SNPs")
dev.off()

# Francesca wants it by gene instead of SNP.
# How many genes?
# Unique based on the gene
# Get cASE SNPs per treatment
Tr_cASE_etoh_h2o$rsID <- sapply(strsplit(Tr_cASE_etoh_h2o$SNP_Individual, "_"), "[", 1)
Tr_cASE_etoh_h2o <- merge(Tr_cASE_etoh_h2o, gene_anno[, c("rsID", "ensg", "gene_type", "g.id")], by="rsID")

cASEbyTreat_gene <- aggregate(cbind(count = g.id) ~ Treatment, data=unique(Tr_cASE_etoh_h2o[,c("g.id", "Treatment")]), FUN = function(x){NROW(x)})

cASEbyTreat_gene <- cASEbyTreat_gene[order(cASEbyTreat_gene$count),]
cASEbyTreat_gene$Treatment <- factor(cASEbyTreat_gene$Treatment, levels=cASEbyTreat_gene$Treatment)


# Instead of separating treatment x cell type between IPSCs and CMs, plot the unique number of cell x treat total (i.e. if SNP has treat x cell for both IPSC and CM, it's only counted once and is probably a IPSC-specific signal).
TrCell_cASE_etoh_h2o_gene$Treatment <- gsub("^Tr", "", sapply(strsplit(TrCell_cASE_etoh_h2o_gene$term, ":"), "[", 2))

# Unique on SNP_Individual and Treatment
cASEbyTreatCell_oneCell <- aggregate(cbind(count = SNP_Individual) ~ Treatment, data=unique(TrCell_cASE_etoh_h2o_gene[,c("SNP_Individual", "Treatment")]), FUN = function(x){NROW(x)})

cASEbyTreatCell_oneCell <- cASEbyTreatCell_oneCell[ order(cASEbyTreatCell_oneCell$count), ]
cASEbyTreatCell_oneCell$Treatment <- factor(cASEbyTreatCell_oneCell$Treatment, levels=cASEbyTreatCell_oneCell$Treatment)

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/allCells_Interact_cASEbyCellTreat.noColor.noX.bar.pdf")
ggplot(data=cASEbyTreatCell_oneCell, aes(x=Treatment, y=count)) +
  geom_bar(stat="identity", color="black", fill="palevioletred1") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "Number of cASE SNPs")
dev.off()

# Combine Treatment cASE with Treatment x Cell Type cASE into one plot
#bar_comb <- rbind(mylm_all_5df_anovaASE_interact_Tr_sig[,c("SNP_Individual", "Treatment")], mylm_all_5df_anovaASE_interact_TrCell_sig_uniq[,c("SNP_Individual", "Treatment")])

#bar_comb <- unique(bar_comb)

#cASEbyComb <- aggregate(cbind(count = SNP_Individual) ~ Treatment, data=bar_comb, FUN = function(x){NROW(x)})

#cASEbyComb <- cASEbyComb[ order(cASEbyComb$count), ]
#cASEbyComb$Treatment <- factor(cASEbyComb$Treatment, levels=cASEbyComb$Treatment)

#pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/allCells_Interact_cASEbyComb.noColor.bar.pdf")
#ggplot(data=cASEbyComb, aes(x=Treatment, y=count)) +
#  geom_bar(stat="identity", color="black", fill="palevioletred1") +
#  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "Number of cASE SNPs")
#dev.off()

# Make stacked bar plot
cASEbyTreat$cASE <- "Treatment"
cASEbyTreatCell_oneCell$cASE <- "Treatment x Cell Type"
cASE_stacked <- rbind(cASEbyTreat, cASEbyTreatCell_oneCell)

cASE_stacked$Treatment <- as.character(cASE_stacked$Treatment)
cASE_stacked <- cASE_stacked[order(cASE_stacked$Treatment),]
cASE_stacked$Treatment <- factor(cASE_stacked$Treatment, levels=unique(cASE_stacked$Treatment))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/allCells_Interact_cASEbyComb.noX.stacked.bar.pdf")
ggplot(data=cASE_stacked, aes(x=Treatment, y=count, fill=cASE)) +
  geom_bar(stat="identity", color="black") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "Number of cASE SNPs")
dev.off()


save(list=c("myanova_all_5df_finite", "mylm_all_5df", "vp_all"), file="interaction_noReadCov_DF_ANOVA_add1_noX/anova_lm_vp.Rd")



###
# Make barplot of number of cASE per SNP
TreatPerSNP <- aggregate(cbind(count = Treatment) ~ SNP_Individual, data=Tr_cASE_etoh_h2o %>% select(Treatment, SNP_Individual) %>% unique(), FUN = function(x){NROW(x)})

sum(TreatPerSNP$count > 1)
[1] 209  # was 217 with X
sum(TreatPerSNP$count <= 1)
[1] 848  # was 907 with X
dim(TreatPerSNP)
[1] 1057    2  # was 1124 with X


pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/allCells_Interact_TreatPerSNP.noX.hist.pdf")
ggplot(data=TreatPerSNP, aes(x=count)) +
  geom_histogram(breaks=seq(0.5, max(TreatPerSNP$count) + 0.5, by=1), color="darkblue", fill="lightblue") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(y = "Number of Treatment cASE per SNP") + 
  scale_x_continuous(breaks=seq(1, max(TreatPerSNP$count), by=1))
dev.off()

# Treatment-by-CellType
TreatCellPerSNP <- aggregate(cbind(count = term) ~ SNP_Individual, data=TrCell_cASE_etoh_h2o %>% select(term, SNP_Individual) %>% unique(), FUN = function(x){NROW(x)})

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/allCells_Interact_TreatCellPerSNP.noX.hist.pdf")
ggplot(data=TreatCellPerSNP, aes(x=count)) +
  geom_histogram(breaks=seq(0.5, max(TreatCellPerSNP$count) + 0.5, by=1), color="darkblue", fill="lightblue") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(y = "Number of Treatment x Cell cASE per SNP") + 
  scale_x_continuous(breaks=seq(1, max(TreatCellPerSNP$count), by=1))
dev.off()






























####################
####################
# Make TWAS heatmap

### Consider significant genes in TWAS. This file only has significant TWAS genes from many studies.
TWAS <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/FastQTL/GxE_mapping/normalized_metas/TWAS_overlap/media-4(2).txt", header=T, stringsAsFactors=)
TWAS$ensg <- gsub("\\..*", "", TWAS$Gene)


# Add TWAS info to cASE
Tr_cASE_TWAS <- merge(Tr_cASE_etoh_h2o, TWAS, by = "ensg")

# I want to make sure everything is unique based on the columns I'm interested in. Those are the term (treatment), Trait (from TWAS), and I'm going to use the SNP_Individual (from ASE). If one gene has cASE in the same treatment from multiple SNPs, then it will be counted multiple times.
Tr_cASE_TWAS_uniq <- unique(Tr_cASE_TWAS[,c("SNP_Individual", "term", "Trait")])

# This has a lot of traits for which we don't have much cASE. Subset for only TWAS which have at least 20 cASE across all conditions.
cASE_TWAS_table <- table(Tr_cASE_TWAS[,c("term", "Trait")])

# How many unique genes with treatment cASE are in CVD TWAS?
unique(Tr_cASE_TWAS[Tr_cASE_TWAS$Trait %in% c("CARDIoGRAM_C4D_CAD", "GLGC_Mc_HDL", "GLGC_Mc_LDL", "GLGC_Mc_TG", "UKB_20002_1065_self_reported_hypertension", "UKB_20002_1473_self_reported_high_cholesterol", "HRGene_HeartRate", "ICBP_DiastolicPressure", "ICBP_SystolicPressure", "MAGNETIC_HDL.C", "MAGNETIC_IDL.TG", "MAGNETIC_LDL.C", "UKB_20002_1094_self_reported_deep_venous_thrombosis_dvt", "UKB_6150_1_Vascular_or_heart_problems_diagnosed_by_doctor_Heart_attack"), "ensg"]) # 169 genes; was 171


length(unique(TWAS$Trait)) # 103 total TWAS
sum(colSums(cASE_TWAS_table) > 0) # 85 with at least 1 cASE; was 85 with X
sum(colSums(cASE_TWAS_table) > 19) # 44 with at least 20 cASE; was 45 with X

# Here are the TWAS with at least 20 cASE
cASE_TWAS_20 <- names(which(colSums(cASE_TWAS_table) > 19))

for_heatmap <- Tr_cASE_TWAS_uniq %>% count(term, Trait)

## Now add just ASE overlap with TWAS
sig_ASE <- myanova_all_5df_finite[myanova_all_5df_finite$padj < 0.1, ]
sig_ASE$rsID <- sapply(strsplit(sig_ASE$SNP_Individual, "_"), "[", 1)

# I need to add the gene annotations to the SNPs
sig_ASE <- merge(sig_ASE, gene_anno[, c("rsID", "ensg", "gene_type", "g.id")], by="rsID")

# Add TWAS info to ASE
sig_ASE_TWAS <- merge(sig_ASE, TWAS, by = "ensg")

# I want to make sure everything is unique based on the columns I'm interested in. Those are ther term Trait (from TWAS), and I'm going to use the SNP_Individual (from ASE). If one gene has cASE in the same treatment from multiple SNPs, then it will be counted multiple times.
sig_ASE_TWAS_uniq <- unique(sig_ASE_TWAS[,c("SNP_Individual", "Trait")])

# I'm going to add this to the cASE heatmap
for_heatmap_ASE <- sig_ASE_TWAS_uniq %>% count(Trait) %>% add_column(term = "ASE")

for_heatmap <- bind_rows(for_heatmap, for_heatmap_ASE)

for_heatmap_20 <- for_heatmap[ for_heatmap$Trait %in% cASE_TWAS_20, ]

# We have so much ASE that it distorts differences between # of cASE overlaps, so now remove ASE
for_heatmap_20_noASE <- for_heatmap_20 %>% filter(term != "ASE")

p <- ggplot(for_heatmap_20_noASE, aes(Trait, term, fill= log(n))) + 
  geom_tile(colour="black") + scale_fill_gradient(low="white", high="blue") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_20cASE.noASE.noX.pdf", width=12)
p
dev.off()

### I'm going to divide the # of cASE by the # of ASE per trait to see if it makes the plot look better
for_heatmap_20_noASE_spread <- for_heatmap_20_noASE %>% spread(term, n)

ASE_toDivide <- for_heatmap_20 %>% filter(term == "ASE") %>% select(Trait, n)

ASE_toDivide_mat<- as.matrix(ASE_toDivide[,"n"])
for_heatmap_20_noASE_spread_mat <- as.matrix(for_heatmap_20_noASE_spread[,2:13])

# Divide all cASE by ASE
cASE_div_ASE <- data.frame(for_heatmap_20_noASE_spread_mat) / do.call("cbind", replicate(12, data.frame(ASE_toDivide_mat), simplify = FALSE))

cASE_div_ASE <- cbind(for_heatmap_20_noASE_spread[,1], cASE_div_ASE) # Add rownames back
colnames(cASE_div_ASE)[1] <- "Trait"

for_heatmap_cASE_div_ASE <- melt(cASE_div_ASE, id="Trait")

p <- ggplot(for_heatmap_cASE_div_ASE, aes(Trait, variable, fill= value)) + 
  geom_tile(colour="black") + scale_fill_gradient(low="white", high="blue") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_20cASE.noASE.divcASE.noX.pdf", width=12)
p
dev.off()

# Remove blood traits because they aren't very interesting and change NA to 0
no_blood <- for_heatmap_cASE_div_ASE[!grepl("Astle", for_heatmap_cASE_div_ASE$Trait),]
no_blood[is.na(no_blood)] <- 0

p <- ggplot(no_blood, aes(Trait, variable, fill= value)) + 
  geom_tile(colour="black") + scale_fill_gradient(low="white", high="blue") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_20cASE.noASE.divcASE.noBlood.noX.pdf", width=12)
p
dev.off()

# Now plot CARDIoGRAM as a barplot per FL's request
no_blood %>% filter(grepl("CARDIoGRAM", Trait)) -> cardiogram_overlap

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_cardiogram.divcASE.noX.box.pdf")
ggplot(data=cardiogram_overlap, aes(x=variable, y=value)) +
  geom_bar(stat="identity", color="black", fill="palevioletred1") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "cASE - CARDIoGRAM overlap")
dev.off()

# I think it might be better just plotting the total numbers of overlapping genes, not any type of proportion
for_heatmap_20 %>% filter(grepl("CARDIoGRAM", Trait), term != "ASE") %>% mutate(term2 = gsub("^Tr", "", term)) -> cardiogram_overlap_num

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_cardiogram.num.noX.bar.pdf")
ggplot(data=cardiogram_overlap_num, aes(x=term2, y=n)) +
  geom_bar(stat="identity", color="black", fill="palevioletred1") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "# of genes")
dev.off()

# Roger wants me to add confidence intervals using binomial distribution. I have percentages in cardiogram_overlap and the total number of ASE is in ASE_toDivide.

# Total number of ASE in CARDIoGRAM:
cardiogram_overlap$num_ASE <- rep(ASE_toDivide[grep("CARDIo", ASE_toDivide$Trait), "n"], 12)

# Total number of cASE per treatment
cardiogram_overlap$num_cASE <- round(as.numeric(cardiogram_overlap$value) * as.numeric(cardiogram_overlap$num_ASE))

# Add confidence interval
cardiogram_overlap <- cbind(cardiogram_overlap, binconf(as.numeric(cardiogram_overlap$num_cASE), as.numeric(cardiogram_overlap$num_ASE)))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_cardiogram.divcASE.confInt.noX.box.pdf")
ggplot(data=cardiogram_overlap, aes(x=variable, y=value, ymin=Lower, ymax=Upper)) +
  geom_bar(stat="identity", color="black", fill="palevioletred1") + geom_errorbar() +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "cASE - CARDIoGRAM overlap")
dev.off()

## Do the same for BMI
for_heatmap_20 %>% filter(grepl("Body_mass_index", Trait), term != "ASE") %>% mutate(term2 = gsub("^Tr", "", term)) -> BMI_overlap_num

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_BMI.num.noX.bar.pdf")
ggplot(data=BMI_overlap_num, aes(x=term2, y=n)) +
  geom_bar(stat="identity", color="black", fill="palevioletred1") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "# of genes")
dev.off()

# Add 6 heart rate genes on top of cardiogram plot
for_heatmap %>% filter(grepl("CARDIoGRAM", Trait) | grepl("HRGene", Trait), term != "ASE") %>% mutate(term2 = gsub("^Tr", "", term)) -> cardiogram_HRGene_overlap_num

cardiogram_HRGene_overlap_num$Trait <- factor(cardiogram_HRGene_overlap_num$Trait, levels=c("HRGene_HeartRate", "CARDIoGRAM_C4D_CAD"))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_cardiogram_HRGene.num.bar.noX.pdf", width=9)
ggplot(data=cardiogram_HRGene_overlap_num, aes(x=term2, y=n, fill=Trait)) +
  geom_bar(stat="identity", color="black") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "# of genes")
dev.off()

# From the big AFib GWAS (Nielsen et al 2018, Nature Genetics), see if any of the priortized genes are cASE genes.
Afib <- read.csv("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/annotation/Nielsen2018_AFib_GWAS_SuppTable11.csv", header=T, stringsAsFactors=F, comment.char="#")
colnames(Afib)[8] <- "Gene_set"

Afib_TrcASE <- merge(Tr_cASE_etoh_h2o, Afib, by.x = "g.id", by.y = "Gene")

# Insulin, Nicotine, Triclosan, and zinc have no cASE overlap, but I still want to show them in the plot
Afib_TrcASE %>% count(term) %>% bind_rows(data.frame(term = c("TrInsulin", "TrNicotine", "TrTriclosan", "TrZinc"), n = c(0,0,0,0))) %>% mutate(term2 = gsub("^Tr", "", term)) %>% arrange(term2) -> Afib_TrcASE_count

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_AFib.num.noX.bar.pdf")
ggplot(data=Afib_TrcASE_count, aes(x=term2, y=n)) +
  geom_bar(stat="identity", color="black", fill="palevioletred1") +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=1, hjust=1, color="black"), text = element_text(size=16), axis.text.y = element_text(color="black")) + labs(x = element_blank(), y = "# of genes")
dev.off()



# prop.test

prop.test(for_heatmap_20_noASE_spread_mat[1,], rep(372, 12))

for_heatmap_20_noASE_spread_mat_noNA <- for_heatmap_20_noASE_spread_mat

for_heatmap_20_noASE_spread_mat_noNA[is.na(for_heatmap_20_noASE_spread_mat_noNA)] <- 0


p_dist <- data.frame(p.val = rep(NA, 44))  # This was 45 before

for (i in c(1:44)){

	test <- prop.test(for_heatmap_20_noASE_spread_mat_noNA[i,], rep(ASE_toDivide_mat[i,], 12))
	
	p_dist[i,1] <- test$p.val
}

rownames(p_dist) <- for_heatmap_20_noASE_spread$Trait

p_dist$padj <- p.adjust(p_dist$p.val, method="BH")

# Roger asked me to plot it a different way: "Let say O(t,c) is the number of overlap of cell-type c cASE with a trait t TWAS. Then N(t)=sum(O(t,…)) sum of all overlaps over all cell-types for a given trait. We can caluclate the poprotion P(t,c) = O(t,c) / N(t), and compare that proportion to P0(c)=mean(  P(…,c)). The proprotion test based on Z-score would be Z(t,c)=P(t,c)-P0(c)/sqrt(1/N(t) * P0(c) * (1-P0(c))).  you can plot the sqrt(1/N(t) * P0(c) * (1-P0(c))) and is the standard error."

# Add number of ASE and cASE
for_heatmap_cASE_div_ASE[is.na(for_heatmap_cASE_div_ASE)] <- 0
for_heatmap_cASE_div_ASE <- merge(for_heatmap_cASE_div_ASE, ASE_toDivide, by="Trait")
colnames(for_heatmap_cASE_div_ASE)[4] <- "num_ASE"

for_heatmap_cASE_div_ASE$num_cASE <- round(as.numeric(for_heatmap_cASE_div_ASE$value) * as.numeric(for_heatmap_cASE_div_ASE$num_ASE))

# First sum of all overlaps over all cell-types for a given trait
N_t_Tr <- for_heatmap_cASE_div_ASE %>% select(Trait, num_cASE) %>% group_by(Trait) %>% summarise(N_t = sum(num_cASE))
forest_plot_Tr <- merge(for_heatmap_cASE_div_ASE, N_t_Tr, by="Trait")
forest_plot_Tr$P_tc <- forest_plot_Tr$num_cASE / forest_plot_Tr$N_t
P0_c_Tr <- forest_plot_Tr %>% select(P_tc, variable) %>% group_by(variable) %>% summarise(P0_c = mean(P_tc))

forest_plot_Tr <- merge(forest_plot_Tr, P0_c_Tr, by="variable")
forest_plot_Tr %>% mutate(Z_tc = (P_tc - P0_c)/sqrt(1/(N_t) * P0_c * (1 - P0_c))) -> forest_plot_Tr
forest_plot_Tr %>% mutate(prop_plot = P_tc - P0_c, prop_se = sqrt(1/(N_t) * P0_c * (1 - P0_c))) -> forest_plot_Tr
forest_plot_Tr %>% mutate(se_above = prop_plot + prop_se, se_below = prop_plot - prop_se) -> forest_plot_Tr

forest_plot_Tr$dummyName <- paste0(forest_plot_Tr$Trait, "_", forest_plot_Tr$variable)

# Plot Roger wanted
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Tr.RogerPropTest.confInt.noX.fp.pdf", height=18, width=10)
ggplot(data=forest_plot_Tr, aes(x=dummyName, y=prop_plot, ymin=se_below, ymax=se_above)) +
    geom_pointrange(aes(color=variable)) +
    coord_flip() + 
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()
dev.off()

# Plot heatmap of Z-scores
p <- ggplot(forest_plot_Tr, aes(Trait, variable, fill=Z_tc)) + 
  geom_tile(colour="black") + scale_fill_gradient2(low="red", mid="white", high="blue") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.RogerPadj.Zscore.heat.noX.pdf", height=5, width=10)
p
dev.off()

# Rename TWAS (in excel) and load in new names
TWAS_names <- read.delim("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/annotation/TWAS_key.tab")
colnames(TWAS_names) <- c("Trait", "Short_Trait")
forest_plot_Tr_shortName <- merge(forest_plot_Tr, TWAS_names, by="Trait")

# Remove "Tr" from treatment names
forest_plot_Tr_shortName$Tr_name <- gsub("^Tr", "", forest_plot_Tr_shortName$variable)
forest_plot_Tr_shortName[forest_plot_Tr_shortName$Tr_name == "Acetaminophin", "Tr_name"] <- "Acetaminophen"

# Plot same heatmap as above
p <- ggplot(forest_plot_Tr_shortName, aes(Short_Trait, Tr_name, fill=Z_tc)) + 
  geom_tile(colour="black") + scale_fill_gradient2(low="blue", mid="white", high="red") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1), axis.text = element_text(size=13))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.RogerPadj.Zscore.final.noX.heat.pdf", height=5, width=12)
p
dev.off()


# Make anything with Z-score < 2 white
forest_plot_Tr_white <- forest_plot_Tr
forest_plot_Tr_white[abs(forest_plot_Tr_white$Z_tc) < 2 , "Z_tc"] <- 0

p <- ggplot(forest_plot_Tr_white, aes(Trait, variable, fill=Z_tc)) + 
  geom_tile(colour="black") + scale_fill_gradient2(low="blue", mid="white", high="red") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.RogerPadj.Zscore.Sig.noX.heat.pdf", height=5, width=10)
p
dev.off()


# Just plot traits with at least 1 significant value by padj
forest_plot_Tr$Roger_pval <- 2*pnorm(-abs(forest_plot_Tr$Z_tc))
forest_plot_Tr$Roger_padj <- p.adjust(forest_plot_Tr$Roger_pval, method="BH")
sig_TWAS_RogerPadj_Tr <-unique(forest_plot_Tr[forest_plot_Tr$Roger_padj < 0.1, "Trait"])
forest_plot_RogerPadj_sigTrait_Tr <- forest_plot_Tr %>% filter(Trait %in% sig_TWAS_RogerPadj_Tr)
dim(forest_plot_RogerPadj_sigTrait_Tr)
[1] 0 17

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.RogerPadj.SigTraits.confInt.noX.fp.pdf", height=5, width=10)
ggplot(data=forest_plot_RogerPadj_sigTrait_Tr, aes(x=dummyName, y=prop_plot, ymin=se_below, ymax=se_above, pch=forest_plot_RogerPadj_sigTrait_Tr$Roger_padj < 0.1)) +
    geom_pointrange(aes(color=variable)) +
    coord_flip() + 
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()
dev.off()


# Set significance threshold to Z-score > 2
sig_TWAS_RogerPropTest_Tr <-unique(forest_plot_Tr[abs(forest_plot_Tr$Z_tc) > 2, "Trait"])
forest_plot_RogerPropTest_sigTrait_Tr <- forest_plot_Tr %>% filter(Trait %in% sig_TWAS_RogerPropTest_Tr)
dim(forest_plot_RogerPropTest_sigTrait_Tr)
[1] 108  16  # Was 120 with X

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Tr.RogerPropTest.SigTraits.confInt.noX.fp.pdf", height=7, width=20)
ggplot(data=forest_plot_RogerPropTest_sigTrait_Tr, aes(x=variable, y=prop_plot, ymin=se_below, ymax=se_above, pch=abs(forest_plot_RogerPropTest_sigTrait_Tr$Z_tc) > 2)) +
    geom_pointrange(aes(color=Trait)) +
	facet_grid(. ~ Trait) +
    coord_flip() + 
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()
dev.off()

# Don't facet
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Tr.RogerPropTest.SigTraits.confInt.noFacet.noX.fp.pdf", height=15, width=10)
ggplot(data=forest_plot_RogerPropTest_sigTrait_Tr, aes(x=dummyName, y=prop_plot, ymin=se_below, ymax=se_above, pch=abs(forest_plot_RogerPropTest_sigTrait_Tr$Z_tc) > 2)) +
    geom_pointrange(aes(color=Trait)) +
    coord_flip() + 
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()
dev.off()

# Just plot sig trait-cASE pairs
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Tr.RogerPropTest.SigTraits.confInt.noFacet.sigOnly.noX.fp.pdf", height=4, width=10)
ggplot(data=forest_plot_RogerPropTest_sigTrait_Tr %>% filter(abs(Z_tc) > 2) , aes(x=dummyName, y=prop_plot, ymin=se_below, ymax=se_above)) +
    geom_pointrange() +
	geom_hline(yintercept=0, linetype="dashed") + 
    coord_flip() + 
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()
dev.off()



# Do TWAS overlap with cell type cASE
# Add TWAS info to cASE
Cell_cASE_TWAS <- merge(Cell_cASE_etoh_h2o_gene, TWAS, by = "ensg")

# I want to make sure everything is unique based on the columns I'm interested in. Those are the term (treatment), Trait (from TWAS), and I'm going to use the SNP_Individual (from ASE). If one gene has cASE in the same treatment from multiple SNPs, then it will be counted multiple times.
Cell_cASE_TWAS_uniq <- unique(Cell_cASE_TWAS[,c("SNP_Individual", "term", "Trait")])

# If there is cell specific cASE in both CMs and IPSCs, call it LCL cASE.

Cell_cASE_LCL <- Cell_cASE_TWAS_uniq %>% mutate(yesno = 1) %>% spread(term, yesno, fill=0) %>% filter(CellTypeCM == "1" & CellTypeIPSC == "1") %>% pull(SNP_Individual) %>% unique()

Cell_cASE_TWAS_uniq[Cell_cASE_TWAS_uniq$SNP_Individual %in% Cell_cASE_LCL,"term"] <- "CellTypeLCL"
Cell_cASE_TWAS_uniq <- unique(Cell_cASE_TWAS_uniq)


# This has a lot of traits for which we don't have much cASE. Subset for only TWAS which have at least 20 cASE across all conditions.
Cell_cASE_TWAS_table <- table(Cell_cASE_TWAS_uniq[,c("term", "Trait")])

sum(colSums(Cell_cASE_TWAS_table) > 0) # 90 with at least 1 cASE; same as X
sum(colSums(Cell_cASE_TWAS_table) > 19) # 59 with at least 20 cASE; same as X

# Here are the TWAS with at least 20 cASE
Cell_cASE_TWAS_20 <- names(which(colSums(Cell_cASE_TWAS_table) > 19))

Cell_for_heatmap <- Cell_cASE_TWAS_uniq %>% count(term, Trait)

# Add ASE values so I can divide cASE as I did for treatment
Cell_for_heatmap_ASE <- Cell_cASE_TWAS_uniq %>% count(Trait) %>% add_column(term = "ASE")

Cell_for_heatmap <- bind_rows(Cell_for_heatmap, Cell_for_heatmap_ASE)

Cell_for_heatmap_20 <- Cell_for_heatmap[ Cell_for_heatmap$Trait %in% Cell_cASE_TWAS_20, ]

# We have so much ASE that it distorts differences between # of cASE overlaps, so now remove ASE
Cell_for_heatmap_20_noASE <- Cell_for_heatmap_20 %>% filter(term != "ASE")

### I'm going to divide the # of cASE by the # of ASE per trait to see if it makes the plot look better
Cell_for_heatmap_20_noASE_spread <- Cell_for_heatmap_20_noASE %>% spread(term, n)

Cell_ASE_toDivide <- Cell_for_heatmap_20 %>% filter(term == "ASE") %>% select(Trait, n)

Cell_ASE_toDivide_mat<- as.matrix(Cell_ASE_toDivide[,"n"])
Cell_for_heatmap_20_noASE_spread_mat <- as.matrix(Cell_for_heatmap_20_noASE_spread[,2:4])

# Divide all cASE by ASE
Cell_cASE_div_ASE <- data.frame(Cell_for_heatmap_20_noASE_spread_mat) / do.call("cbind", replicate(3, data.frame(Cell_ASE_toDivide_mat), simplify = FALSE))

Cell_cASE_div_ASE <- cbind(Cell_for_heatmap_20_noASE_spread[,1], Cell_cASE_div_ASE) # Add rownames back
colnames(Cell_cASE_div_ASE)[1] <- "Trait"

Cell_for_heatmap_cASE_div_ASE <- melt(Cell_cASE_div_ASE, id="Trait")

p <- ggplot(Cell_for_heatmap_cASE_div_ASE, aes(Trait, variable, fill= value)) + 
  geom_tile(colour="black") + scale_fill_gradient(low="white", high="blue") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell_20cASE.noASE.divcASE.noX.pdf", width=12)
p
dev.off()

# Remove blood traits because they aren't very interesting and change NA to 0
Cell_no_blood <- Cell_for_heatmap_cASE_div_ASE[!grepl("Astle", Cell_for_heatmap_cASE_div_ASE$Trait),]
Cell_no_blood[is.na(Cell_no_blood)] <- 0

p <- ggplot(Cell_no_blood, aes(Trait, variable, fill= value)) + 
  geom_tile(colour="black") + scale_fill_gradient(low="white", high="blue") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell_20cASE.noASE.divcASE.noBlood.noX.pdf", width=12)
p
dev.off()

### Turn the heatmap into a forest plot with error bars?
# Roger wants me to add confidence intervals using binomial distribution. I have percentages in Cell_for_heatmap_cASE_div_ASE and the total number of ASE is in ASE_toDivide.

forest_plot <- Cell_for_heatmap_cASE_div_ASE

# Total number of ASE in CARDIoGRAM:
Cell_ASE_toDivide$Trait <- as.character(Cell_ASE_toDivide$Trait)
forest_plot$Trait <- as.character(forest_plot$Trait)
forest_plot <- merge(forest_plot, Cell_ASE_toDivide, by = "Trait")

# Total number of cASE per treatment
forest_plot$num_cASE <- round(as.numeric(forest_plot$value) * as.numeric(forest_plot$n))

# Add confidence interval
forest_plot <- cbind(forest_plot, binconf(as.numeric(forest_plot$num_cASE), as.numeric(forest_plot$n)))

forest_plot$dummyName <- paste0(forest_plot$Trait, "_", forest_plot$variable)

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.divcASE.confInt.noX.fp.pdf", height=20, width=10)
ggplot(data=forest_plot, aes(x=dummyName, y=PointEst, ymin=Lower, ymax=Upper)) +
    geom_pointrange(aes(color=variable)) +
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
	#facet_grid(Trait ~ .) +
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()# + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=for_plot$Plate.ID), width=0.5)
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.divcASE.confInt.facetTrait.noX.fp.pdf", height=50, width=5)
ggplot(data=forest_plot, aes(x=variable, y=PointEst, ymin=Lower, ymax=Upper)) +
    geom_pointrange(aes(color=variable)) +
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
	facet_grid(Trait ~ .) +
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()# + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=for_plot$Plate.ID), width=0.5)
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.divcASE.confInt.facetCell.noX.fp.pdf", height=15, width=15)
ggplot(data=forest_plot, aes(x=Trait, y=PointEst, ymin=Lower, ymax=Upper)) +
    geom_pointrange(aes(color=variable)) +
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
	facet_grid(. ~ variable) +
    xlab("Trait") + ylab("Proportion (95% CI)") + theme_bw()# + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=for_plot$Plate.ID), width=0.5)
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.divcASE.confInt.facetWrap.noX.fp.pdf", height=10, width=20)
ggplot(data=forest_plot, aes(x=variable, y=PointEst, ymin=Lower, ymax=Upper)) +
    geom_pointrange(aes(color=variable)) +
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
	facet_wrap(vars(Trait)) +
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()# + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=for_plot$Plate.ID), width=0.5)
dev.off()


p_dist_cell <- data.frame(p.val = rep(NA, 59))

for (i in c(1:59)){

	test <- prop.test(Cell_for_heatmap_20_noASE_spread_mat[i,], rep(Cell_ASE_toDivide_mat[i,], 3))
	
	p_dist_cell[i,1] <- test$p.val
}

p_dist_cell$Trait <- Cell_for_heatmap_20_noASE_spread$Trait

p_dist_cell$padj <- p.adjust(p_dist_cell$p.val, method="BH")

### 
# Refine forest plots. Roger likes the one with no facets, colored by cell type, with trait-cell type on y axis. He wants me to only keep interesting TWAS. I'm going to remove blood TWAS and any that are not significant.

forest_plot_p <- merge(forest_plot, p_dist_cell, by="Trait")
forest_plot_p_noBlood <- forest_plot_p[!grepl("Astle",forest_plot_p$Trait), ]

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.noBloodSig.divcASE.confInt.noX.fp.pdf", height=20, width=10)
ggplot(data=forest_plot_p_noBlood %>% filter(padj < 0.1), aes(x=dummyName, y=PointEst, ymin=Lower, ymax=Upper)) +
    geom_pointrange(aes(color=variable)) +
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
	#facet_grid(Trait ~ .) +
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()# + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=for_plot$Plate.ID), width=0.5)
dev.off()

# Now try further removing baldness, gout, sleep, Neuroticism, intelligence, Education, Inflammatory Bowel (similar to UC), breast cancer
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.noBloodSig_Select.divcASE.confInt.noX.fp.pdf", height=12, width=10)
ggplot(data=forest_plot_p_noBlood %>% filter(padj < 0.1) %>% filter(!grepl("Breast|Bowel|Education|Sleep|Morning|Gout|intelligence|Neuroticism|Hair", Trait)), aes(x=dummyName, y=PointEst, ymin=Lower, ymax=Upper)) +
    geom_pointrange(aes(color=variable)) +
    coord_flip() + 
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()
dev.off()

# Francesca preferred the one facetted by cell type
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.noBloodSig_Select.divcASE.confInt.facetCell.noX.fp.pdf", height=8, width=15)
ggplot(data=forest_plot_p_noBlood %>% filter(padj < 0.1) %>% filter(!grepl("Breast|Bowel|Education|Sleep|Morning|Gout|intelligence|Neuroticism|Hair", Trait)), aes(x=Trait, y=PointEst, ymin=Lower, ymax=Upper)) +
    geom_pointrange(aes(color=variable)) +
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
	facet_grid(. ~ variable) +
    xlab("Trait") + ylab("Proportion (95% CI)") + theme_bw()# + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=for_plot$Plate.ID), width=0.5)
dev.off()

###
# Roger asked me to plot it a different way: "Let say O(t,c) is the number of overlap of cell-type c cASE with a trait t TWAS. Then N(t)=sum(O(t,…)) sum of all overlaps over all cell-types for a given trait. We can caluclate the poprotion P(t,c) = O(t,c) / N(t), and compare that proportion to P0(c)=mean(  P(…,c)). The proprotion test based on Z-score would be Z(t,c)=P(t,c)-P0(c)/sqrt(1/N(t) * P0(c) * (1-P0(c))).  you can plot the sqrt(1/N(t) * P0(c) * (1-P0(c))) and is the standard error."

# First sum of all overlaps over all cell-types for a given trait
N_t <- forest_plot %>% select(Trait, num_cASE) %>% group_by(Trait) %>% summarise(N_t = sum(num_cASE))
forest_plot <- merge(forest_plot, N_t, by="Trait")
forest_plot$P_tc <- forest_plot$num_cASE / forest_plot$N_t
P0_c <- forest_plot %>% select(P_tc, variable) %>% group_by(variable) %>% summarise(P0_c = mean(P_tc))

forest_plot <- merge(forest_plot, P0_c, by="variable")
forest_plot %>% mutate(Z_tc = (P_tc - P0_c)/sqrt(1/(N_t) * P0_c * (1 - P0_c))) -> forest_plot
forest_plot %>% mutate(prop_plot = P_tc - P0_c, prop_se = sqrt(1/(N_t) * P0_c * (1 - P0_c))) -> forest_plot
forest_plot %>% mutate(se_above = prop_plot + prop_se, se_below = prop_plot - prop_se) -> forest_plot

# Plot Roger wanted
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.RogerPropTest.confInt.noX.fp.pdf", height=18, width=10)
ggplot(data=forest_plot, aes(x=dummyName, y=prop_plot, ymin=se_below, ymax=se_above)) +
    geom_pointrange(aes(color=variable)) +
    coord_flip() + 
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()
dev.off()

# Just plot traits with at least 1 significant value by padj
forest_plot$Roger_pval <- 2*pnorm(-abs(forest_plot$Z_tc))
forest_plot$Roger_padj <- p.adjust(forest_plot$Roger_pval, method="BH")
sig_TWAS_RogerPadj <-unique(forest_plot[forest_plot$Roger_padj < 0.1, "Trait"])
forest_plot_RogerPadj_sigTrait <- forest_plot %>% filter(Trait %in% sig_TWAS_RogerPadj)
dim(forest_plot_RogerPadj_sigTrait)
[1] 15 19  # Same as with X

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.RogerPadj.SigTraits.confInt.noX.fp.pdf", height=5, width=10)
ggplot(data=forest_plot_RogerPadj_sigTrait, aes(x=dummyName, y=prop_plot, ymin=se_below, ymax=se_above, pch=forest_plot_RogerPadj_sigTrait$Roger_padj < 0.1)) +
    geom_pointrange(aes(color=variable)) +
    coord_flip() + 
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()
dev.off()


# Set significance threshold to Z-score > 2
sig_TWAS_RogerPropTest <-unique(forest_plot[abs(forest_plot$Z_tc) > 2, "Trait"])
forest_plot_RogerPropTest_sigTrait <- forest_plot %>% filter(Trait %in% sig_TWAS_RogerPropTest)
dim(forest_plot_RogerPropTest_sigTrait)
[1] 45 19  # Same as with X

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.RogerPropTest.SigTraits.confInt.noX.fp.pdf", height=10, width=10)
ggplot(data=forest_plot_RogerPropTest_sigTrait, aes(x=dummyName, y=prop_plot, ymin=se_below, ymax=se_above, pch=abs(forest_plot_RogerPropTest_sigTrait$Z_tc) > 2)) +
    geom_pointrange(aes(color=variable)) +
    coord_flip() + 
    xlab("Cell Type") + ylab("Proportion (95% CI)") + theme_bw()
dev.off()


# Francesca preferred the one facetted by cell type
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.RogerPropTest.confInt.facetCell.noX.fp.pdf", height=12, width=15)
ggplot(data=forest_plot, aes(x=Trait, y=prop_plot, ymin=se_below, ymax=se_above)) +
    geom_pointrange(aes(color=variable)) +
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
	facet_grid(. ~ variable) +
    xlab("Trait") + ylab("Proportion (95% CI)") + theme_bw()# + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=for_plot$Plate.ID), width=0.5)
dev.off()

### Plot just CM cASE and only traits related to CVD
forest_plot %>% filter(variable == "CellTypeCM",
	Trait %in% c("CARDIoGRAM_C4D_CAD", "MAGNETIC_LDL.C", "GLGC_Mc_HDL", "UKB_20002_1065_self_reported_hypertension", "GLGC_Mc_LDL", "GLGC_Mc_TG", "UKB_20002_1473_self_reported_high_cholesterol", "HRGene_HeartRate")) -> forest_plot_CM

# Correct p-values just for CM cASE and CVD
forest_plot_CM$Roger_padj_CM <- p.adjust(forest_plot_CM$Roger_pval, method="BH")
forest_plot_CM <- forest_plot_CM[order(forest_plot_CM$prop_plot),]
forest_plot_CM$Trait <- factor(forest_plot_CM$Trait, levels=forest_plot_CM$Trait)

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.CM.RogerPadj.SigTraits.confInt.noX.fp.pdf", height=3, width=8)
ggplot(data=forest_plot_CM, aes(x=Trait, y=prop_plot, ymin=se_below, ymax=se_above)) +
    geom_pointrange(aes(color=variable, pch=forest_plot_CM$Roger_padj_CM < 0.1)) +
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Trait") + ylab("Proportion (95% CI)") + theme_bw()# + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=for_plot$Plate.ID), width=0.5)
dev.off()

# Use color to signify significance
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell.CM.RogerPadj.SigTraits.confInt.color.noX.fp.pdf", height=3, width=7)
ggplot(data=forest_plot_CM, aes(x=Trait, y=prop_plot, ymin=se_below, ymax=se_above)) +
    geom_pointrange(aes(color=Roger_padj_CM < 0.1)) +
    #geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
    coord_flip() +  # flip coordinates (puts labels on y axis)
    xlab("Trait") + ylab("Proportion (95% CI)") + theme_bw()# + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=for_plot$Plate.ID), width=0.5)
dev.off()


###############################################################
###############################################################

# Now I'm going to try for an enrichment.

# First get all genes which had ASE in the linear model.
# Tr_cASE <- mylm_all_5df_ASE_interact_Tr[ mylm_all_5df_ASE_interact_Tr$padj < 0.1, ]
mylm_all_5df_anovaASE_interact_Tr$rsID <- sapply(strsplit(mylm_all_5df_anovaASE_interact_Tr$SNP_Individual, "_"), "[", 1)

mylm_all_5df_anovaASE_interact_Tr_gene <- merge(mylm_all_5df_anovaASE_interact_Tr, gene_anno[, c("rsID", "ensg", "gene_type", "g.id")], by="rsID")

# Keep just the entries which have one water and ethanol per cell type
mylm_all_5df_anovaASE_interact_Tr_gene %>% filter(SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) -> Tr_cASE__allTested_etoh_h2o

genesTestedForCASE <- unique(Tr_cASE__allTested_etoh_h2o$ensg)

# I have the TWAS(+)cASE(+) entry from the overlap. If I get the total TWAS(+) for each trait, then I can easily get the TWAS(+)cASE(-) entry.

TWAS_uniq <- unique(TWAS[,c("Trait", "ensg")])

# Subset to only genes tested for cASE
TWAS_uniq <- TWAS_uniq[TWAS_uniq$ensg %in% genesTestedForCASE, ]

TWAS_uniq_count <- TWAS_uniq %>% count(Trait)

# Get number of total cASE genes per treatment
Tr_cASE_uniq <- unique(Tr_cASE_etoh_h2o[,c("term", "ensg")])
Tr_cASE_count <- Tr_cASE_uniq %>% count(term)


## Try an example
# Treat <- "TrAcetaminophin"
# TWAS_trait <- "Astle_et_al_2016_Lymphocyte_counts"

# Create an output file. Begin with a header
sink(file="TWAS_cASE_enrichment.noX.tab")
cat(paste0("Treatment", "\t", "Trait", "\t", "Odds_Ratio", "\t", "Low", "\t", "Hi", "\t", "p.value", "\t", "cASE_P_TWAS_P", "\t", "cASE_N_TWAS_P", "\t", "cASE_P_TWAS_N", "\t", "cASE_N_TWAS_N", "\n"))
sink()

for (i in c(1:519)){  # Was 530
Treat <- for_heatmap_20_noASE[i, "term"] %>% unlist(use.names=F)
TWAS_trait <- for_heatmap_20_noASE[i, "Trait"] %>% unlist(use.names=F)


# cASE(+)TWAS(+)
cASE_P_TWAS_P <- for_heatmap_20_noASE[ which(for_heatmap_20_noASE$term == Treat & as.character(for_heatmap_20_noASE$Trait) == TWAS_trait), "n"] %>% unlist(use.names=F)

# cASE(-)TWAS(+): Total in TWAS - cASE(+)TWAS(+)
cASE_N_TWAS_P <- TWAS_uniq_count[TWAS_uniq_count$Trait == TWAS_trait, "n"] %>% unlist(use.names=F) - cASE_P_TWAS_P

# cASE(+)TWAS(-): Total cASE - cASE_P_TWAS_P
cASE_P_TWAS_N <- Tr_cASE_count[Tr_cASE_count$term == Treat, "n"] %>% unlist(use.names=F) - cASE_P_TWAS_P

# cASE(-)TWAS(-): Total genes tested - other 3 values
cASE_N_TWAS_N <- length(genesTestedForCASE) - cASE_P_TWAS_P - cASE_N_TWAS_P - cASE_P_TWAS_N

### Do Fisher's test to see if significantly different
enrich <-
matrix(c(cASE_P_TWAS_P, cASE_N_TWAS_P, cASE_P_TWAS_N, cASE_N_TWAS_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       TWAS = c("TWAS+", "TWAS-")))
fish_test <- fisher.test(enrich, alternative = "two.sided")
names(fish_test$estimate) <- "Odds_Ratio"

sink(file="TWAS_cASE_enrichment.noX.tab", append=TRUE)
cat(paste0(Treat, "\t", TWAS_trait, "\t", fish_test$estimate, "\t", fish_test$conf.int[1], "\t", fish_test$conf.int[2], "\t", fish_test$p.value, "\t", cASE_P_TWAS_P, "\t", cASE_N_TWAS_P, "\t", cASE_P_TWAS_N, "\t", cASE_N_TWAS_N, "\n"))
sink()

print(i)

}

# Now load the data in
fishtest <- read.table("TWAS_cASE_enrichment.noX.tab", header=T, stringsAsFactors=F, sep="\t")

# BH correct p.value
fishtest$padj <- p.adjust(fishtest$p.value, method="BH")

write.table(fishtest, file="TWAS_cASE_enrichment.padj.tab", col.names=T, row.names=F, quote=F, sep="\t")

# Make a plot where color is odds ratio, but everything with p > 0.05 should be blank
fishtest$adj_OR <- fishtest$Odds_Ratio
fishtest[fishtest$p.value > 0.05, "adj_OR"] <- 0

p <- ggplot(fishtest, aes(Trait, Treatment, fill= adj_OR)) + 
  geom_tile(colour="black") + scale_fill_gradient(low="white", high="blue") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_20cASE.noASE.odds_ratio.noX.pdf", width=12)
p
dev.off()

##############################################################
##############################################################
##############################################################
# Try cell-type cASE
# I think I need to go from long to wide. I want to have the SNP_Indiviual, trait, and then 0/1 for CelltypeIPSC column and CellTypeCM column. Add dummy column with 1's for spread.
Cell_cASE_TWAS_uniq$value <- 1

cell_cASE_TWAS_spread <- spread(Cell_cASE_TWAS_uniq, key=term, value=value)

cell_cASE_TWAS_spread[is.na(cell_cASE_TWAS_spread)] <- 0

# Create column which will indicate ifcASE is in CMs, IPSCs, or both
cell_cASE_TWAS_spread$Cell <- "Test"
cell_cASE_TWAS_spread[cell_cASE_TWAS_spread$CellTypeCM == 1 & cell_cASE_TWAS_spread$CellTypeIPSC == 1, "Cell"] <- "Both"
cell_cASE_TWAS_spread[cell_cASE_TWAS_spread$CellTypeCM == 0 & cell_cASE_TWAS_spread$CellTypeIPSC == 1, "Cell"] <- "IPSC"
cell_cASE_TWAS_spread[cell_cASE_TWAS_spread$CellTypeCM == 1 & cell_cASE_TWAS_spread$CellTypeIPSC == 0, "Cell"] <- "CM"

# See how many of each cell type goes with each trait
for_heatmap_cell <- cell_cASE_TWAS_spread %>% count(Cell, Trait)

p <- ggplot(for_heatmap_cell, aes(Trait, Cell, fill= n)) + 
  geom_tile(colour="black") + scale_fill_gradient(low="white", high="blue") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/cASE_TWAS_Cell_20cASE.noASE.noX.pdf", width=12)
p
dev.off()


# Now I'm going to try for an enrichment.
# Get number of total cell cASE genes per treatment
cell_cASE_uniq <- unique(Cell_cASE_etoh_h2o[,c("term", "ensg")])
cell_cASE_count <- cell_cASE_uniq %>% count(term)


## Try an example
# Treat <- "TrAcetaminophin"
# TWAS_trait <- "Astle_et_al_2016_Lymphocyte_counts"

##### I THINK THIS WAS AN ERROR TO INCLUDE HERE, SO I'VE REMOVED IT
# Create an output file. Begin with a header
#sink(file="TWAS_cASE_enrichment.noX.tab")
#cat(paste0("Treatment", "\t", "Trait", "\t", "Odds_Ratio", "\t", "Low", "\t", "Hi", "\t", "p.value", "\t", "cASE_P_TWAS_P", "\t", "cASE_N_TWAS_P", "\t", "cASE_P_TWAS_N", "\t", "cASE_N_TWAS_N", "\n"))
#sink()

#for (i in c(1:519)){  # Was 661
#Treat <- for_heatmap_20_noASE[i, "term"] %>% unlist(use.names=F)
#TWAS_trait <- for_heatmap_20_noASE[i, "Trait"] %>% unlist(use.names=F)


# cASE(+)TWAS(+)
#cASE_P_TWAS_P <- for_heatmap_20_noASE[ which(for_heatmap_20_noASE$term == Treat & as.character(for_heatmap_20_noASE$Trait) == TWAS_trait), "n"] %>% unlist(use.names=F)

# cASE(-)TWAS(+): Total in TWAS - cASE(+)TWAS(+)
#cASE_N_TWAS_P <- TWAS_uniq_count[TWAS_uniq_count$Trait == TWAS_trait, "n"] %>% unlist(use.names=F) - cASE_P_TWAS_P

# cASE(+)TWAS(-): Total cASE - cASE_P_TWAS_P
#cASE_P_TWAS_N <- Tr_cASE_count[Tr_cASE_count$term == Treat, "n"] %>% unlist(use.names=F) - cASE_P_TWAS_P

# cASE(-)TWAS(-): Total genes tested (3353) - other 3 values
#cASE_N_TWAS_N <- 3353 - cASE_P_TWAS_P - cASE_N_TWAS_P - cASE_P_TWAS_N

### Do Fisher's test to see if significantly different
#enrich <-
#matrix(c(cASE_P_TWAS_P, cASE_N_TWAS_P, cASE_P_TWAS_N, cASE_N_TWAS_N),
#       nrow = 2,
#       dimnames = list(cASE = c("cASE+", "cASE-"),
#                       TWAS = c("TWAS+", "TWAS-")))
#fish_test <- fisher.test(enrich, alternative = "two.sided")
#names(fish_test$estimate) <- "Odds_Ratio"

#sink(file="TWAS_cASE_enrichment.noX.tab", append=TRUE)
#cat(paste0(Treat, "\t", TWAS_trait, "\t", fish_test$estimate, "\t", fish_test$conf.int[1], "\t", fish_test$conf.int[2], "\t", fish_test$p.value, "\t", cASE_P_TWAS_P, "\t", cASE_N_TWAS_P, "\t", cASE_P_TWAS_N, "\t", cASE_N_TWAS_N, "\n"))
#sink()

#print(i)

#}

# Now load the data in
#fishtest <- read.table("TWAS_cASE_enrichment.noX.tab", header=T, stringsAsFactors=F, sep="\t")

# BH correct p.value
#fishtest$padj <- p.adjust(fishtest$p.value, method="BH")

#write.table(fishtest, file="TWAS_cASE_enrichment.noX.padj.tab", col.names=T, row.names=F, quote=F, sep="\t")



######################################
######################################
# Output cASE and ASE SNPs for Alan. Only give him SNPs which are not on X chromosome
all_anova_ASE_forWrite <- data.frame(ASE_SNPs = unique(gsub("_.*", "", all_anova_ASE)))
all_anova_ASE_forWrite <- merge(all_anova_ASE_forWrite, gene_anno, by.x="ASE_SNPs", by.y="rsID")
all_anova_ASE_forWrite %>% filter(chr != "X") ->  all_anova_ASE_forWrite
all_anova_ASE_forWrite %>% select(ASE_SNPs) -> all_anova_ASE_forWrite
all_anova_ASE_forWrite <- unique(all_anova_ASE_forWrite)
write.table(all_anova_ASE_forWrite, file="interaction_noReadCov_DF_ANOVA_add1_noX/ASE_SNPs.txt", quote=F, sep="\t", row.names=F, col.names=F)

Tr_cASE_forWrite <- data.frame(cASE_SNPs = unique(gsub("_.*", "", Tr_cASE_etoh_h2o_2$SNP_Individual)))
Tr_cASE_forWrite <- merge(Tr_cASE_forWrite, gene_anno, by.x="cASE_SNPs", by.y="rsID")
Tr_cASE_forWrite %>% filter(chr != "X") ->  Tr_cASE_forWrite
Tr_cASE_forWrite %>% select(cASE_SNPs) -> Tr_cASE_forWrite
Tr_cASE_forWrite <- unique(Tr_cASE_forWrite)
write.table(Tr_cASE_forWrite, file="interaction_noReadCov_DF_ANOVA_add1_noX/Tr_cASE_SNPs.txt", quote=F, sep="\t", row.names=F, col.names=F)


#################################
#################################
######## Overlap with GxE database
gxe_db <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/FastQTL/GxE_mapping/normalized_metas/external_validation/datasets/GxE_full_db_build1/GxE_db_full.txt", header=T, sep="\t", stringsAsFactors=F)

gxe_db %>% filter(significant | pvalue < 0.05) -> gxe_db_05

sum( (unique(Tr_cASE_etoh_h2o$g.id) %in% gxe_db_05$g.id) | (unique(Tr_cASE_etoh_h2o$ensg) %in% gxe_db_05$ensg) )
[1] 850  # was 855 with X

length(unique(Tr_cASE_etoh_h2o$g.id))
[1] 979  # was 997 with X

# Make a table of the results
# The GxE table has g.id's for some entries, and ensg's for others. Using my genes with cASE, add missing names.
Tr_cASE_gene_anno <- unique(Tr_cASE_etoh_h2o[,c("g.id", "ensg")])

gxe_db_05 %>% filter(is.na(ensg)) %>% select(-ensg) %>% inner_join(Tr_cASE_gene_anno) -> gxe_db_05_gid
gxe_db_05 %>% filter(is.na(g.id)) %>% select(-g.id) %>% inner_join(Tr_cASE_gene_anno) -> gxe_db_05_ensg
gxe_db_05 %>% filter(!is.na(ensg) & !is.na(g.id) & ensg %in% Tr_cASE_gene_anno$ensg & g.id %in% Tr_cASE_gene_anno$g.id) -> gxe_db_05_both

gxe_db_05_TrcASE <- rbind(gxe_db_05_both, gxe_db_05_ensg, gxe_db_05_gid)


Tr_cASE_etoh_h2o %>% inner_join(gxe_db_05_TrcASE, by=c("ensg", "g.id")) -> Tr_cASE_gxe

Tr_cASE_gxe %>% select(g.id, Treatment, condition) %>% distinct() %>% group_by(Treatment, condition) %>% summarise(num = n())



####################
# How many genes with ASE have multiple SNPs for which we can measure ASE within them?

ASE_genes <- data.frame(SNP_Individual = all_anova_ASE, rsID = sapply(strsplit(all_anova_ASE, "_"), "[", 1), Individual = sapply(strsplit(all_anova_ASE, "_"), "[", 2))
ASE_genes  <- merge(ASE_genes, gene_anno[, c("rsID", "ensg", "gene_type", "g.id")], by="rsID")

# Add gene annotation to all SNPs tested for ASE
all_genes_tested <- data.frame(SNP_Individual = unique(as.character(by_SNP.Ind_testable$SNP_Individual)), rsID = sapply(strsplit(unique(as.character(by_SNP.Ind_testable$SNP_Individual)), "_"), "[", 1), Individual = sapply(strsplit(unique(as.character(by_SNP.Ind_testable$SNP_Individual)), "_"), "[", 2))
all_genes_tested  <- merge(all_genes_tested, gene_anno[, c("rsID", "ensg", "gene_type", "g.id")], by="rsID")

# Only keep the genes which had ASE. We want to do this for each gene in an individual, so create a new variable for that
ASE_genes$gene_ind <- paste0(ASE_genes$g.id, "_", ASE_genes$Individual)
all_genes_tested$gene_ind <- paste0(all_genes_tested$g.id, "_", all_genes_tested$Individual)

ASE_genes_tested <- all_genes_tested %>% filter(gene_ind %in% unique(ASE_genes$gene_ind))

# Make sure we only have unique rsID - gene_ind pairs
ASE_genes_tested <- ASE_genes_tested[!duplicated(ASE_genes_tested[,c("gene_ind", "rsID")]),]

SNPs_per_gene <- aggregate(rsID ~ gene_ind, ASE_genes_tested, FUN = function(x){NROW(x)})






##########################
##########################
# Write files for output

# Generate file with treatment cASE intersected with nominally significant GxE
Tr_cASE_gxe_write <- Tr_cASE_gxe %>% select(ensg, g.id, SNP_Individual, Treatment, estimate, std.error, statistic, p.value, padj, deg.f, rsID.y, condition, cell.type, effect.size, Bayes.factor, pvalue, p.adjusted, significant, study, DOI, Suppl.Table)
write.table(Tr_cASE_gxe_write, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/Tr_cASE_nomSigGxE.noX.tab", col.names=T, row.names=F, quote=F, sep="\t")

# Linear model results
mylm_all_5df_anovaASE_write <- mylm_all_5df_anovaASE %>% select(SNP_Individual, term, estimate, std.error, statistic, p.value, deg.f)
write.table(mylm_all_5df_anovaASE, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/LinearModel_allCells.noX.tab", col.names=T, row.names=F, quote=F, sep="\t")

# ANOVA results
myanova_all_5df_finite_write <- myanova_all_5df_finite %>% select(SNP_Individual, Res.Df, RSS, Df, Sum.Sq, F.value, p.value, padj)
write.table(myanova_all_5df_finite_write, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/AnovaModel_allCells.noX.tab", col.names=T, row.names=F, quote=F, sep="\t")

# 44 TWAS I need to rename
#write.table(data.frame(x = unique(forest_plot_Tr$Trait)), file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/annotation/TWAS_key.noX.tab", quote=F, row.names=F, col.names=F, sep="\t")

# Afib TWAS overlap
Afib_TrcASE_write <- Afib_TrcASE %>% mutate(treatment = gsub("^Tr", "", term)) %>% select(ensg, g.id, SNP_Individual, treatment, estimate, std.error, statistic, p.value, padj, deg.f)
write.table(Afib_TrcASE_write, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/Tr_cASE_Afib.noX.tab", col.names=T, row.names=F, quote=F, sep="\t")

# Mixed effect model
vp_all_ASE_df <- data.frame(vp_all_ASE)
vp_all_ASE_df$SNP_Individual <- rownames(vp_all_ASE_df)
vp_all_ASE_df <- vp_all_ASE_df %>% select(SNP_Individual, Tr, CellType, CellType.Tr, Control.ID, Residual)
write.table(vp_all_ASE_df, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/MixedEffect_allCell.tab", row.names=F, col.names=T, quote=F, sep="\t")



#########################
#########################
# Make table that Roger wants with every term from linear model, number of tests, and number with padj < 10%

## Control
# How many SNPs tested for control
mylm_all_5df_anovaASE_control %>% filter(SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) %>% pull(SNP_Individual) %>% unique() %>% length() # 13009
# How many sig?
mylm_all_5df_anovaASE_control %>% filter(SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh & padj < 0.1) %>% pull(SNP_Individual) %>% unique() %>% length() # 0

## Intercept
mylm_all_5df_anovaASE_intercept <- mylm_all_5df_anovaASE[grepl("(Intercept)", mylm_all_5df_anovaASE$term),] 
mylm_all_5df_anovaASE_intercept$padj <- p.adjust(mylm_all_5df_anovaASE_intercept$p.value, method="BH")
# How many SNPs tested for intercept
mylm_all_5df_anovaASE_intercept %>% filter(SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) %>% pull(SNP_Individual) %>% unique() %>% length() # 13009
# How many sig?
mylm_all_5df_anovaASE_intercept %>% filter(SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh & padj < 0.1) %>% pull(SNP_Individual) %>% unique() %>% length() # 0

## Treatment
mylm_all_5df_anovaASE_interact_Tr %>% filter(SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) %>% group_by(term) %>% summarise(n = n())

mylm_all_5df_anovaASE_interact_Tr %>% filter(SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh & padj < 0.1) %>% group_by(term) %>% summarise(n = n())
