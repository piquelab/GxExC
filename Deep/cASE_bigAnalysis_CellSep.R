#Anthony Findley
# 5/8/2020

# Purpose: Call ASE based on ANOVA model. This is for the model where cell types are analyzed separately.
# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/fixed_effect

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

# Function to implement linear model and include DF
mylm <- function(data){
	SNP_Individual=paste0(data$rsID[1],"_",data$Individual[1])
	obj <- lm(beta1 ~ Control.ID + Tr + PlateVar, data = data)
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
	full = lm(beta1 ~ Control.ID + Tr + PlateVar, data=data)
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

####################################
####################################
### Do CMs first
ASE_CM <- ASE %>% filter(CellType == "CM")
ASE_CM$Tr <- as.character(ASE_CM$Treatment.Name)
ASE_CM$Tr[ASE_CM$Tr=="Ethanol"] = "Control"
ASE_CM$Tr[ASE_CM$Tr=="Water"] = "Control"
ASE_CM$Tr = factor(ASE_CM$Tr)
ASE_CM$Tr = relevel(ASE_CM$Tr,ref="Control")
ASE.small <- ASE_CM[,c("SNP_Individual", "Tr")]
ASE.small.Tr <- dcast(ASE.small, SNP_Individual ~ Tr) # Creates df for every SNP_Individual, number of CM, IPSC, and LCL treatments
rownames(ASE.small.Tr) <- ASE.small.Tr$SNP_Individual
ASE.small.Tr <- ASE.small.Tr[,c(2:14)]
ASE.small.Tr[ASE.small.Tr > 0] <-1

# Only keep these SNPs to test
testable_SNP_Ind <- rownames(ASE.small.Tr[which(rowSums(ASE.small.Tr) > 4), ])
ASE_CM$SNP_Individual <- as.character(ASE_CM$SNP_Individual)
ASE_testable <- ASE_CM[ASE_CM$SNP_Individual %in% testable_SNP_Ind,]
ASE_testable$SNP_Individual <- factor(ASE_testable$SNP_Individual)
by_SNP.Ind <- group_by(ASE_testable, SNP_Individual)

# Now make sure every SNP has at least one of each control
ASE.small.Ctrl <- dcast(ASE_testable, SNP_Individual ~ Control.ID)
testable_SNP_Ind <- ASE.small.Ctrl[which(ASE.small.Ctrl$CO1 > 0 & ASE.small.Ctrl$CO2 > 0), "SNP_Individual"]
ASE_testable$SNP_Individual <- as.character(ASE_testable$SNP_Individual)
ASE_testable <- ASE_testable[ASE_testable$SNP_Individual %in% testable_SNP_Ind,]
ASE_testable$SNP_Individual <- factor(ASE_testable$SNP_Individual)
by_SNP.Ind <- group_by(ASE_testable, SNP_Individual)

### Now try adding plate variable
by_SNP.Ind$PlateVar <- factor(gsub(".*R", "R", by_SNP.Ind$Plate.ID))
ASE.small.Plate <- dcast(by_SNP.Ind, SNP_Individual ~ PlateVar)
testable_SNP_Ind <- ASE.small.Plate[which(ASE.small.Plate$R1 > 0 & ASE.small.Plate$R2 > 0), "SNP_Individual"]
by_SNP.Ind <- by_SNP.Ind[which(by_SNP.Ind$SNP_Individual %in% testable_SNP_Ind),]

# Run the linear model
aux <- by_SNP.Ind %>% nest()  
mylist <- aux$data
names(mylist) <- aux$SNP_Individual
mylm_CM <- map_dfr(mylist,mylm)

# Run the ANOVA
myanova_reduced_CM <- map_dfr(mylist,myanova_reduced)

#########
# Now do LCLs
ASE_LCL <- ASE %>% filter(CellType == "LCL")
ASE_LCL$Tr <- as.character(ASE_LCL$Treatment.Name)
ASE_LCL$Tr[ASE_LCL$Tr=="Ethanol"] = "Control"
ASE_LCL$Tr[ASE_LCL$Tr=="Water"] = "Control"
ASE_LCL$Tr = factor(ASE_LCL$Tr)
ASE_LCL$Tr = relevel(ASE_LCL$Tr,ref="Control")
ASE.small <- ASE_LCL[,c("SNP_Individual", "Tr")]
ASE.small.Tr <- dcast(ASE.small, SNP_Individual ~ Tr) # Creates df for every SNP_Individual, number of LCL, IPSC, and LCL treatments
rownames(ASE.small.Tr) <- ASE.small.Tr$SNP_Individual
ASE.small.Tr <- ASE.small.Tr[,c(2:14)]
ASE.small.Tr[ASE.small.Tr > 0] <-1

# Only keep these SNPs to test
testable_SNP_Ind <- rownames(ASE.small.Tr[which(rowSums(ASE.small.Tr) > 4), ])
ASE_LCL$SNP_Individual <- as.character(ASE_LCL$SNP_Individual)
ASE_testable <- ASE_LCL[ASE_LCL$SNP_Individual %in% testable_SNP_Ind,]
ASE_testable$SNP_Individual <- factor(ASE_testable$SNP_Individual)
by_SNP.Ind <- group_by(ASE_testable, SNP_Individual)

# Now make sure every SNP has at least one of each control
ASE.small.Ctrl <- dcast(ASE_testable, SNP_Individual ~ Control.ID)
testable_SNP_Ind <- ASE.small.Ctrl[which(ASE.small.Ctrl$CO1 > 0 & ASE.small.Ctrl$CO2 > 0), "SNP_Individual"]
ASE_testable$SNP_Individual <- as.character(ASE_testable$SNP_Individual)
ASE_testable <- ASE_testable[ASE_testable$SNP_Individual %in% testable_SNP_Ind,]
ASE_testable$SNP_Individual <- factor(ASE_testable$SNP_Individual)
by_SNP.Ind <- group_by(ASE_testable, SNP_Individual)

### Now try adding plate variable
by_SNP.Ind$PlateVar <- factor(gsub(".*R", "R", by_SNP.Ind$Plate.ID))
ASE.small.Plate <- dcast(by_SNP.Ind, SNP_Individual ~ PlateVar)
testable_SNP_Ind <- ASE.small.Plate[which(ASE.small.Plate$R1 > 0 & ASE.small.Plate$R2 > 0), "SNP_Individual"]
by_SNP.Ind <- by_SNP.Ind[which(by_SNP.Ind$SNP_Individual %in% testable_SNP_Ind),]

# Run the linear model
aux <- by_SNP.Ind %>% nest()  
mylist <- aux$data
names(mylist) <- aux$SNP_Individual
mylm_LCL <- map_dfr(mylist,mylm)

# Run the ANOVA
myanova_reduced_LCL <- map_dfr(mylist,myanova_reduced)

#########
# Now do IPSCs
ASE_IPSC <- ASE[which(ASE$CellType == "IPSC"),]
ASE_IPSC$Tr <- as.character(ASE_IPSC$Treatment.Name)
ASE_IPSC$Tr[ASE_IPSC$Tr=="Ethanol"] = "Control"
ASE_IPSC$Tr[ASE_IPSC$Tr=="Water"] = "Control"
ASE_IPSC$Tr = factor(ASE_IPSC$Tr)
ASE_IPSC$Tr = relevel(ASE_IPSC$Tr,ref="Control")
ASE.small <- ASE_IPSC[,c("SNP_Individual", "Tr")]
ASE.small.Tr <- dcast(ASE.small, SNP_Individual ~ Tr) # Creates df for every SNP_Individual, number of IPSC, IPSC, and IPSC treatments
rownames(ASE.small.Tr) <- ASE.small.Tr$SNP_Individual
ASE.small.Tr <- ASE.small.Tr[,c(2:14)]
ASE.small.Tr[ASE.small.Tr > 0] <-1

# Only keep these SNPs to test
testable_SNP_Ind <- rownames(ASE.small.Tr[which(rowSums(ASE.small.Tr) > 4), ])
ASE_IPSC$SNP_Individual <- as.character(ASE_IPSC$SNP_Individual)
ASE_testable <- ASE_IPSC[ASE_IPSC$SNP_Individual %in% testable_SNP_Ind,]
ASE_testable$SNP_Individual <- factor(ASE_testable$SNP_Individual)
by_SNP.Ind <- group_by(ASE_testable, SNP_Individual)

# Now make sure every SNP has at least one of each control
ASE.small.Ctrl <- dcast(ASE_testable, SNP_Individual ~ Control.ID)
testable_SNP_Ind <- ASE.small.Ctrl[which(ASE.small.Ctrl$CO1 > 0 & ASE.small.Ctrl$CO2 > 0), "SNP_Individual"]
ASE_testable$SNP_Individual <- as.character(ASE_testable$SNP_Individual)
ASE_testable <- ASE_testable[ASE_testable$SNP_Individual %in% testable_SNP_Ind,]
ASE_testable$SNP_Individual <- factor(ASE_testable$SNP_Individual)
by_SNP.Ind <- group_by(ASE_testable, SNP_Individual)

### Now try adding plate variable
by_SNP.Ind$PlateVar <- factor(gsub(".*R", "R", by_SNP.Ind$Plate.ID))
ASE.small.Plate <- dcast(by_SNP.Ind, SNP_Individual ~ PlateVar)
testable_SNP_Ind <- ASE.small.Plate[which(ASE.small.Plate$R1 > 0 & ASE.small.Plate$R2 > 0), "SNP_Individual"]
by_SNP.Ind <- by_SNP.Ind[which(by_SNP.Ind$SNP_Individual %in% testable_SNP_Ind),]

# Run the linear model
aux <- by_SNP.Ind %>% nest()  
mylist <- aux$data
names(mylist) <- aux$SNP_Individual
mylm_IPSC <- map_dfr(mylist,mylm)

# Run the ANOVA
aux <- by_SNP.Ind %>% nest()  
mylist <- aux$data
names(mylist) <- aux$SNP_Individual
myanova_reduced_IPSC <- map_dfr(mylist,myanova_reduced)

#####################################
#####################################
# Analyze results from linear model. Start with CMs

# Remove values with NA for p.value
mylm_CM <- mylm_CM[which(!is.na(mylm_CM$p.value)),]
mylm_CM_5df <- mylm_CM[which(mylm_CM$deg.f > 4),]
mylm_CM_5df_intercept <- mylm_CM_5df[which(mylm_CM_5df$term == "(Intercept)"),]
dim(mylm_CM_5df_intercept)
[1]  67496    8 # 68160 with X
mylm_CM_5df_intercept <- mylm_CM_5df_intercept[order(mylm_CM_5df_intercept$p.value),]
mylm_CM_5df_intercept$padj <- p.adjust(mylm_CM_5df_intercept$p.value, method="BH")
sum(mylm_CM_5df_intercept$padj < 0.1) # How much ASE with FDR < 0.1
[1] 2522 # was 2653 with X
mylm_CM_5df_intercept$exp <- ppoints(length(mylm_CM_5df_intercept$term))
### Make qq-plots for treatments
mylm_CM_5df_Tr <- mylm_CM_5df[grep("Tr", mylm_CM_5df$term),]
mylm_CM_5df_Tr <- mylm_CM_5df_Tr[order(mylm_CM_5df_Tr$p.value),]
# Subset cASE for SNPs which display ASE
CM_ASE <- mylm_CM_5df_intercept[which(mylm_CM_5df_intercept$padj < 0.1),"SNP_Individual"]
mylm_CM_5df_Tr_ASE <- mylm_CM_5df_Tr[which(mylm_CM_5df_Tr$SNP_Individual %in% CM_ASE),]
mylm_CM_5df_Tr_ASE$padj <- p.adjust(mylm_CM_5df_Tr_ASE$p.value, method="BH")
mylm_CM_5df_Tr_ASE <- mylm_CM_5df_Tr_ASE[order(mylm_CM_5df_Tr_ASE$p.value),]
mylm_CM_5df_Tr_ASE$exp <- ppoints(length(mylm_CM_5df_Tr_ASE$term))
### Make qq-plots for plate
mylm_CM_5df_plate <- mylm_CM_5df[grep("Plate", mylm_CM_5df$term),]
# Subset plate for SNPs which display ASE
mylm_CM_5df_plate_ASE <- mylm_CM_5df_plate[which(mylm_CM_5df_plate$SNP_Individual %in% as.character(CM_ASE)),]
mylm_CM_5df_plate_ASE <- mylm_CM_5df_plate_ASE[order(mylm_CM_5df_plate_ASE$p.value),]
mylm_CM_5df_plate_ASE$padj <- p.adjust(mylm_CM_5df_plate_ASE$p.value, method="BH")
mylm_CM_5df_plate_ASE$exp <- ppoints(length(mylm_CM_5df_plate_ASE$p.value))

# LCLs
mylm_LCL <- mylm_LCL[which(!is.na(mylm_LCL$p.value)),]
mylm_LCL_5df <- mylm_LCL[which(mylm_LCL$deg.f > 4),]
mylm_LCL_5df_intercept <- mylm_LCL_5df[which(mylm_LCL_5df$term == "(Intercept)"),]

# IPSCs
mylm_IPSC <- mylm_IPSC[which(!is.na(mylm_IPSC$p.value)),]
mylm_IPSC_5df <- mylm_IPSC[which(mylm_IPSC$deg.f > 4),]
mylm_IPSC_5df_intercept <- mylm_IPSC_5df[which(mylm_IPSC_5df$term == "(Intercept)"),]



######################################
######################################
# Forget the intercept model for now. Analyze the ANOVA since that's what we're interested in.
# Start with CMs

# Only keep SNPs with 5 degrees of freedom from linear model.
myanova_reduced_CM_5df <- myanova_reduced_CM[ myanova_reduced_CM$SNP_Individual %in% unique(mylm_CM_5df$SNP_Individual), ]
myanova_reduced_CM_5df_finite <- myanova_reduced_CM_5df[ is.finite(myanova_reduced_CM_5df$p.value),]

# Now make plot
myanova_reduced_CM_5df_finite <- myanova_reduced_CM_5df_finite[order(myanova_reduced_CM_5df_finite$p.value),]
myanova_reduced_CM_5df_finite$padj <- p.adjust(myanova_reduced_CM_5df_finite$p.value, method="BH")
sum(myanova_reduced_CM_5df_finite$padj < 0.1) # How much ASE with FDR < 0.1
[1] 6811 # There was 6980 with X

myanova_reduced_CM_5df_finite$exp <- ppoints(length(myanova_reduced_CM_5df_finite$SNP_Individual))

#########################################################
# LCLs 

# Only keep SNPs with 5 degrees of freedom from linear model.
myanova_reduced_LCL_5df <- myanova_reduced_LCL[ myanova_reduced_LCL$SNP_Individual %in% unique(mylm_LCL_5df$SNP_Individual), ]
myanova_reduced_LCL_5df_finite <- myanova_reduced_LCL_5df[ is.finite(myanova_reduced_LCL_5df$p.value),]

# Now make plot
myanova_reduced_LCL_5df_finite <- myanova_reduced_LCL_5df_finite[order(myanova_reduced_LCL_5df_finite$p.value),]
myanova_reduced_LCL_5df_finite$padj <- p.adjust(myanova_reduced_LCL_5df_finite$p.value, method="BH")
sum(myanova_reduced_LCL_5df_finite$padj < 0.1) # How much ASE with FDR < 0.1
[1] 7127 # There was 7877 with X

myanova_reduced_LCL_5df_finite$exp <- ppoints(length(myanova_reduced_LCL_5df_finite$SNP_Individual))

#########################################################
# IPSCs 

# Only keep SNPs with 5 degrees of freedom from linear model.
myanova_reduced_IPSC_5df <- myanova_reduced_IPSC[ myanova_reduced_IPSC$SNP_Individual %in% unique(mylm_IPSC_5df$SNP_Individual), ]
myanova_reduced_IPSC_5df_finite <- myanova_reduced_IPSC_5df[ is.finite(myanova_reduced_IPSC_5df$p.value),]

# Now make plot
myanova_reduced_IPSC_5df_finite <- myanova_reduced_IPSC_5df_finite[order(myanova_reduced_IPSC_5df_finite$p.value),]
myanova_reduced_IPSC_5df_finite$padj <- p.adjust(myanova_reduced_IPSC_5df_finite$p.value, method="BH")
sum(myanova_reduced_IPSC_5df_finite$padj < 0.1) # How much ASE with FDR < 0.1
[1] 6794 # There was 7048 with X

myanova_reduced_IPSC_5df_finite$exp <- ppoints(length(myanova_reduced_IPSC_5df_finite$SNP_Individual))

######
# Plot 3 cell types together
myanova_reduced_CM_5df_finite$cell <- "CM"
myanova_reduced_IPSC_5df_finite$cell <- "IPSC"
myanova_reduced_LCL_5df_finite$cell <- "LCL"

myanova_reduced_allCell <- rbind(myanova_reduced_LCL_5df_finite, myanova_reduced_IPSC_5df_finite, myanova_reduced_CM_5df_finite)

p <- ggplot(myanova_reduced_allCell, aes(-log10(exp), -log10(p.value), color=cell))
p + geom_point(size=2) + labs(title = "ANOVA reduced vs full model", x = "-log10(exp.pvalue)", y = "-log10(obs.pvalue)") + geom_abline(intercept=0, slope=1) + theme_bw() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggsave("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/AllSeparate_ANOVA_reduced.5df.noReadCov.qq.noX.png")


########################################
########################################
# Previously for the fixed effect model, we subset cASE based on ASE determined by the intercept, which we now know is not a good idea. So now use the ANOVA ASE to subset cASE and see how the numbers change (5df cutoff is still in effect).

### Make qq-plots for treatments
CM_anova_ASE <- unique(myanova_reduced_CM_5df_finite[myanova_reduced_CM_5df_finite$padj < 0.1, "SNP_Individual"])
mylm_CM_5df_Tr <- mylm_CM_5df[grep("Tr", mylm_CM_5df$term),]
mylm_CM_5df_Tr <- mylm_CM_5df_Tr[order(mylm_CM_5df_Tr$p.value),]

# Subset cASE for SNPs which display ASE
mylm_CM_5df_Tr_anovaASE <- mylm_CM_5df_Tr[which(mylm_CM_5df_Tr$SNP_Individual %in% CM_anova_ASE),]
mylm_CM_5df_Tr_anovaASE$padj <- p.adjust(mylm_CM_5df_Tr_anovaASE$p.value, method="BH")
mylm_CM_5df_Tr_anovaASE <- mylm_CM_5df_Tr_anovaASE[order(mylm_CM_5df_Tr_anovaASE$p.value),]
mylm_CM_5df_Tr_anovaASE$exp <- ppoints(length(mylm_CM_5df_Tr_anovaASE$term))
mylm_CM_5df_Tr_anovaASE$CellType <- "CM"


# LCLs
LCL_anova_ASE <- unique(myanova_reduced_LCL_5df_finite[myanova_reduced_LCL_5df_finite$padj < 0.1, "SNP_Individual"])
mylm_LCL_5df_Tr <- mylm_LCL_5df[grep("Tr", mylm_LCL_5df$term),]
mylm_LCL_5df_Tr <- mylm_LCL_5df_Tr[order(mylm_LCL_5df_Tr$p.value),]
mylm_LCL_5df_Tr_anovaASE <- mylm_LCL_5df_Tr[which(mylm_LCL_5df_Tr$SNP_Individual %in% LCL_anova_ASE),]
mylm_LCL_5df_Tr_anovaASE$padj <- p.adjust(mylm_LCL_5df_Tr_anovaASE$p.value, method="BH")
mylm_LCL_5df_Tr_anovaASE <- mylm_LCL_5df_Tr_anovaASE[order(mylm_LCL_5df_Tr_anovaASE$p.value),]
mylm_LCL_5df_Tr_anovaASE$exp <- ppoints(length(mylm_LCL_5df_Tr_anovaASE$term))
mylm_LCL_5df_Tr_anovaASE$CellType <- "LCL"

# IPSCs
IPSC_anova_ASE <- unique(myanova_reduced_IPSC_5df_finite[myanova_reduced_IPSC_5df_finite$padj < 0.1, "SNP_Individual"])
mylm_IPSC_5df_Tr <- mylm_IPSC_5df[grep("Tr", mylm_IPSC_5df$term),]
mylm_IPSC_5df_Tr <- mylm_IPSC_5df_Tr[order(mylm_IPSC_5df_Tr$p.value),]
mylm_IPSC_5df_Tr_anovaASE <- mylm_IPSC_5df_Tr[which(mylm_IPSC_5df_Tr$SNP_Individual %in% IPSC_anova_ASE),]
mylm_IPSC_5df_Tr_anovaASE$padj <- p.adjust(mylm_IPSC_5df_Tr_anovaASE$p.value, method="BH")
mylm_IPSC_5df_Tr_anovaASE <- mylm_IPSC_5df_Tr_anovaASE[order(mylm_IPSC_5df_Tr_anovaASE$p.value),]
mylm_IPSC_5df_Tr_anovaASE$exp <- ppoints(length(mylm_IPSC_5df_Tr_anovaASE$term))
mylm_IPSC_5df_Tr_anovaASE$CellType <- "IPSC"

# See how many significant interactions and SNP_Individuals
sum(mylm_CM_5df_Tr_anovaASE$padj < 0.1) # 492; was 441 with X
length(unique(mylm_CM_5df_Tr_anovaASE[ mylm_CM_5df_Tr_anovaASE$padj < 0.1, "SNP_Individual"])) # 352; was 320 with X

sum(mylm_IPSC_5df_Tr_anovaASE$padj < 0.1) # 511; was 462 with X
length(unique(mylm_IPSC_5df_Tr_anovaASE[ mylm_IPSC_5df_Tr_anovaASE$padj < 0.1, "SNP_Individual"])) # 350; was was 321 with X

sum(mylm_LCL_5df_Tr_anovaASE$padj < 0.1) # 230; was 141 with X
length(unique(mylm_LCL_5df_Tr_anovaASE[ mylm_LCL_5df_Tr_anovaASE$padj < 0.1, "SNP_Individual"])) # 168; was 104 with X

# Make QQ-plot
allCellsSeparate_Tr_anovaASE <- rbind(mylm_LCL_5df_Tr_anovaASE, mylm_IPSC_5df_Tr_anovaASE, mylm_CM_5df_Tr_anovaASE)

p <- ggplot(allCellsSeparate_Tr_anovaASE, aes(-log10(exp), -log10(p.value), color=CellType))
p + geom_point(size=2) + labs(title = "Treatment cASE using ANOVA ASE SNPs", x = "-log10(exp.pvalue)", y = "-log10(obs.pvalue)") + geom_abline(intercept=0, slope=1) + theme_bw() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggsave("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/AllSeparate_Tr.anovaASE.5df.noReadCov.qq.noX.png")

### See if the plate effect is significant for any
### Make qq-plots for treatments
mylm_CM_5df_PlateVar <- mylm_CM_5df[grep("PlateVar", mylm_CM_5df$term),]
mylm_CM_5df_PlateVar <- mylm_CM_5df_PlateVar[order(mylm_CM_5df_PlateVar$p.value),]

# Subset cASE for SNPs which display ASE
mylm_CM_5df_PlateVar_anovaASE <- mylm_CM_5df_PlateVar[which(mylm_CM_5df_PlateVar$SNP_Individual %in% CM_anova_ASE),]
mylm_CM_5df_PlateVar_anovaASE$padj <- p.adjust(mylm_CM_5df_PlateVar_anovaASE$p.value, method="BH")
mylm_CM_5df_PlateVar_anovaASE <- mylm_CM_5df_PlateVar_anovaASE[order(mylm_CM_5df_PlateVar_anovaASE$p.value),]
mylm_CM_5df_PlateVar_anovaASE$exp <- ppoints(length(mylm_CM_5df_PlateVar_anovaASE$term))

# LCLs
mylm_LCL_5df_PlateVar <- mylm_LCL_5df[grep("PlateVar", mylm_LCL_5df$term),]
mylm_LCL_5df_PlateVar <- mylm_LCL_5df_PlateVar[order(mylm_LCL_5df_PlateVar$p.value),]
mylm_LCL_5df_PlateVar_anovaASE <- mylm_LCL_5df_PlateVar[which(mylm_LCL_5df_PlateVar$SNP_Individual %in% LCL_anova_ASE),]
mylm_LCL_5df_PlateVar_anovaASE$padj <- p.adjust(mylm_LCL_5df_PlateVar_anovaASE$p.value, method="BH")
mylm_LCL_5df_PlateVar_anovaASE <- mylm_LCL_5df_PlateVar_anovaASE[order(mylm_LCL_5df_PlateVar_anovaASE$p.value),]
mylm_LCL_5df_PlateVar_anovaASE$exp <- ppoints(length(mylm_LCL_5df_PlateVar_anovaASE$term))
mylm_LCL_5df_PlateVar_anovaASE$CellType <- "LCL"

# IPSCs
mylm_IPSC_5df_PlateVar <- mylm_IPSC_5df[grep("PlateVar", mylm_IPSC_5df$term),]
mylm_IPSC_5df_PlateVar <- mylm_IPSC_5df_PlateVar[order(mylm_IPSC_5df_PlateVar$p.value),]
mylm_IPSC_5df_PlateVar_anovaASE <- mylm_IPSC_5df_PlateVar[which(mylm_IPSC_5df_PlateVar$SNP_Individual %in% IPSC_anova_ASE),]
mylm_IPSC_5df_PlateVar_anovaASE$padj <- p.adjust(mylm_IPSC_5df_PlateVar_anovaASE$p.value, method="BH")
mylm_IPSC_5df_PlateVar_anovaASE <- mylm_IPSC_5df_PlateVar_anovaASE[order(mylm_IPSC_5df_PlateVar_anovaASE$p.value),]
mylm_IPSC_5df_PlateVar_anovaASE$exp <- ppoints(length(mylm_IPSC_5df_PlateVar_anovaASE$term))
mylm_IPSC_5df_PlateVar_anovaASE$CellType <- "IPSC"

# See how many significant interactions and SNP_Individuals
sum(mylm_CM_5df_PlateVar_anovaASE$padj < 0.1) # 68; was 72 with X
sum(mylm_IPSC_5df_PlateVar_anovaASE$padj < 0.1) # 98; was 96 with X
sum(mylm_LCL_5df_PlateVar_anovaASE$padj < 0.1) # 0; was 0 with X

# Make QQ-plot
mylm_CM_5df_PlateVar_anovaASE$CellType <- "CM"
allCellsSeparate_PlateVar_anovaASE <- rbind(mylm_LCL_5df_PlateVar_anovaASE, mylm_IPSC_5df_PlateVar_anovaASE, mylm_CM_5df_PlateVar_anovaASE)

p <- ggplot(allCellsSeparate_PlateVar_anovaASE, aes(-log10(exp), -log10(p.value), color=CellType))
p + geom_point(size=2) + labs(title = "Plate replicate cASE using ANOVA ASE SNPs", x = "-log10(exp.pvalue)", y = "-log10(obs.pvalue)") + geom_abline(intercept=0, slope=1) + theme_bw() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggsave("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/AllSeparate_PlateVar.anovaASE.5df.noReadCov.qq.noX.png")


######################################

CM_anova_Tr_cASE <- mylm_CM_5df_Tr_anovaASE[mylm_CM_5df_Tr_anovaASE$padj < 0.1,]
CM_anova_Tr_cASE <- unique(paste0(CM_anova_Tr_cASE$term, "_", CM_anova_Tr_cASE$SNP_Individual))

LCL_anova_Tr_cASE <- mylm_LCL_5df_Tr_anovaASE[mylm_LCL_5df_Tr_anovaASE$padj < 0.1,]
LCL_anova_Tr_cASE <- unique(paste0(LCL_anova_Tr_cASE$term, "_", LCL_anova_Tr_cASE$SNP_Individual))

IPSC_anova_Tr_cASE <- mylm_IPSC_5df_Tr_anovaASE[mylm_IPSC_5df_Tr_anovaASE$padj < 0.1,]
IPSC_anova_Tr_cASE <- unique(paste0(IPSC_anova_Tr_cASE$term, "_", IPSC_anova_Tr_cASE$SNP_Individual))

CM_anova_cASE <- unique(mylm_CM_5df_Tr_anovaASE[mylm_CM_5df_Tr_anovaASE$padj < 0.1, "SNP_Individual"])
LCL_anova_cASE <- unique(mylm_LCL_5df_Tr_anovaASE[mylm_LCL_5df_Tr_anovaASE$padj < 0.1, "SNP_Individual"])
IPSC_anova_cASE <- unique(mylm_IPSC_5df_Tr_anovaASE[mylm_IPSC_5df_Tr_anovaASE$padj < 0.1, "SNP_Individual"])

######
# Plot 


forest_plot_CM <- function(snpID){

	for_plot <- ASE_CM[which(ASE_CM$SNP_Individual == snpID),]

    for_plot$CI_hi <- for_plot$beta1 + (1.96 * for_plot$beta.se)
    for_plot$CI_lo <- for_plot$beta1 - (1.96 * for_plot$beta.se)
    for_plot$Treatment.ID <- factor(for_plot$Treatment.ID)

	for_plot <- for_plot[order(for_plot$Treatment.Name, for_plot$beta1),]
    for_plot$id <- paste0(for_plot$Plate.ID, "_", for_plot$Treatment.Name)
    for_plot$id <- factor(for_plot$id, levels=for_plot$id)

    fp <- ggplot(data=for_plot, aes(x=id, y=beta1, ymin=CI_lo, ymax=CI_hi)) +
    geom_pointrange(aes(color=Treatment.Name)) +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Plate and Treatment") + ylab("ASE (95% CI)") + theme_bw() + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=Treatment.Name), width=0.5) + ggtitle(paste0(snpID, " in CMs"))

    pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/forest/", snpID, ".CM.fp.noX.pdf"))
        print(fp)
    dev.off()
		
	# Now plot the beta's from the model
	beta_plot <- mylm_CM_5df_Tr_anovaASE[which(mylm_CM_5df_Tr_anovaASE$SNP_Individual == snpID),]
    beta_plot$CI_hi <- beta_plot$estimate + (1.96 * beta_plot$std.error)
    beta_plot$CI_lo <- beta_plot$estimate - (1.96 * beta_plot$std.error)

	beta_plot <- beta_plot[order(beta_plot$estimate),]
    beta_plot$term <- factor(beta_plot$term, levels=beta_plot$term)

    fp <- ggplot(data=beta_plot, aes(x=term, y=estimate, ymin=CI_lo, ymax=CI_hi)) +
        geom_pointrange(aes(color=padj < 0.1)) +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Plate and Treatment") + ylab("ASE (95% CI)") + theme_bw() + theme(axis.text=element_text(size=14)) + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=padj < 0.1), width=0.5) + ggtitle(paste0(snpID, " in CMs"))

    pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/forest/", snpID, ".CM.model_betas.fp.noX.pdf"),height=3.5)
        print(fp)
    dev.off()
}
	
#### Make plots for LCLs
forest_plot_LCL <- function(snpID){

	for_plot <- ASE_LCL[which(ASE_LCL$SNP_Individual == snpID),]

    for_plot$CI_hi <- for_plot$beta1 + (1.96 * for_plot$beta.se)
    for_plot$CI_lo <- for_plot$beta1 - (1.96 * for_plot$beta.se)
    for_plot$Treatment.ID <- factor(for_plot$Treatment.ID)

	for_plot <- for_plot[order(for_plot$Treatment.Name, for_plot$beta1),]
    for_plot$id <- paste0(for_plot$Plate.ID, "_", for_plot$Treatment.Name)
    for_plot$id <- factor(for_plot$id, levels=for_plot$id)

    fp <- ggplot(data=for_plot, aes(x=id, y=beta1, ymin=CI_lo, ymax=CI_hi)) +
    geom_pointrange(aes(color=Treatment.Name)) +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Plate and Treatment") + ylab("ASE (95% CI)") + theme_bw() + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=Treatment.Name), width=0.5) + ggtitle(paste0(snpID, " in LCLs"))

    pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/forest/", snpID, ".LCL.fp.noX.pdf"))
        print(fp)
    dev.off()
		
	# Now plot the beta's from the model
	beta_plot <- mylm_LCL_5df_Tr_anovaASE[which(mylm_LCL_5df_Tr_anovaASE$SNP_Individual == snpID),]
    beta_plot$CI_hi <- beta_plot$estimate + (1.96 * beta_plot$std.error)
    beta_plot$CI_lo <- beta_plot$estimate - (1.96 * beta_plot$std.error)

	beta_plot <- beta_plot[order(beta_plot$estimate),]
    beta_plot$term <- factor(beta_plot$term, levels=beta_plot$term)

    fp <- ggplot(data=beta_plot, aes(x=term, y=estimate, ymin=CI_lo, ymax=CI_hi)) +
        geom_pointrange(aes(color=padj < 0.1)) +
        geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Plate and Treatment") + ylab("ASE (95% CI)") + theme_bw() + theme(axis.text=element_text(size=14)) + geom_errorbar(aes(ymin=CI_lo, ymax=CI_hi, col=padj < 0.1), width=0.5) + ggtitle(paste0(snpID, " in LCLs"))

    pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/forest/", snpID, ".LCL.model_betas.fp.noX.pdf"), height=3.5)
        print(fp)
    dev.off()
}

forest_plot_CM("rs1142851_GM19209")
forest_plot_CM("rs182840333_GM18858")
forest_plot_CM("rs112436800_GM18855")

forest_plot_LCL("rs1046165_GM18858")
forest_plot_LCL("rs2280349_GM18855")
forest_plot_LCL("rs116444657_GM18912")


# Now plot SNPs which are significant, but are not the top cASE
forest_plot_CM("rs3889_GM18855")


#####################
######################
# I'm going to run the mixed effect model, but only for SNPs with ASE

calcVarPart <- function(data){
	data$Tr <- factor(data$Tr)
	data$Control.ID <- factor(data$Control.ID)
	data$PlateVar <- factor(data$PlateVar)
	obj <- lmer(beta1 ~ 1 + (1|Control.ID) + (1|Tr) + (1|PlateVar), data=data, control=lmerControl(check.conv.hess="stop", check.conv.grad="stop"))
	# browser()
	conv <- obj@optinfo$conv$opt # I think "0" means that it converged
	aux <- as.data.frame(VarCorr(obj))
	vp <- c(aux$vcov, conv)
	names(vp) <- c(aux$grp, "convergence")
	vp
}


ASE_testable_CM <- ASE_CM[ASE_CM$SNP_Individual %in% CM_anova_ASE,]
ASE_testable_CM$SNP_Individual <- factor(ASE_testable_CM$SNP_Individual)
ASE_testable_CM$PlateVar <- factor(gsub(".*R", "R", ASE_testable_CM$Plate.ID))
by_SNP.Ind_CM <- group_by(ASE_testable_CM, SNP_Individual)

options(nwarnings = 10000)
# Now run model on SNPs tested previously
by_SNP.Ind_CM %>% nest() %>% mutate(vp=map(data, possibly(calcVarPart, otherwise = NULL))) -> aux_CM

names(aux_CM$vp) <- aux_CM$SNP_Individual
vp_CM <- do.call(rbind,aux_CM$vp)
vp_CM <- data.frame(vp_CM)
vp_CM$SNP_Individual <- rownames(vp_CM)

# Now do LCLs
ASE_testable_LCL <- ASE_LCL[ASE_LCL$SNP_Individual %in% LCL_anova_ASE,]
ASE_testable_LCL$SNP_Individual <- factor(ASE_testable_LCL$SNP_Individual)
ASE_testable_LCL$PlateVar <- factor(gsub(".*R", "R", ASE_testable_LCL$Plate.ID))
by_SNP.Ind_LCL <- group_by(ASE_testable_LCL, SNP_Individual)
by_SNP.Ind_LCL %>% nest() %>% mutate(vp=map(data, possibly(calcVarPart, otherwise = NULL))) -> aux_LCL
names(aux_LCL$vp) <- aux_LCL$SNP_Individual
vp_LCL <- do.call(rbind,aux_LCL$vp)
vp_LCL <- data.frame(vp_LCL)
vp_LCL$SNP_Individual <- rownames(vp_LCL)

# Now do IPSCs
ASE_testable_IPSC <- ASE_IPSC[ASE_IPSC$SNP_Individual %in% IPSC_anova_ASE,]
ASE_testable_IPSC$SNP_Individual <- factor(ASE_testable_IPSC$SNP_Individual)
ASE_testable_IPSC$PlateVar <- factor(gsub(".*R", "R", ASE_testable_IPSC$Plate.ID))
by_SNP.Ind_IPSC <- group_by(ASE_testable_IPSC, SNP_Individual)
by_SNP.Ind_IPSC %>% nest() %>% mutate(vp=map(data, possibly(calcVarPart, otherwise = NULL))) -> aux_IPSC
names(aux_IPSC$vp) <- aux_IPSC$SNP_Individual
vp_IPSC <- do.call(rbind,aux_IPSC$vp)
vp_IPSC <- data.frame(vp_IPSC)
vp_IPSC$SNP_Individual <- rownames(vp_IPSC)

### Put everything together and make plots
vp_LCL$CellType <- rep("LCL", length(vp_LCL$Tr))
vp_IPSC$CellType <- rep("IPSC", length(vp_IPSC$Tr))
vp_CM$CellType <- rep("CM", length(vp_CM$Tr))

vp_all <- rbind(vp_LCL, vp_IPSC, vp_CM)
vp_all <- vp_all[vp_all$convergence == 0,] # Remove things which didn't converge

# Now calculate percent variance explained
vp_all_percent <- vp_all
vp_all_percent$Tr <- vp_all_percent$Tr / rowSums(vp_all[,c("Tr", "PlateVar", "Control.ID", "Residual")])
vp_all_percent$PlateVar <- vp_all_percent$PlateVar / rowSums(vp_all[,c("Tr", "PlateVar", "Control.ID", "Residual")])
vp_all_percent$Control.ID <- vp_all_percent$Control.ID / rowSums(vp_all[,c("Tr", "PlateVar", "Control.ID", "Residual")])
vp_all_percent$Residual <- vp_all_percent$Residual / rowSums(vp_all[,c("Tr", "PlateVar", "Control.ID", "Residual")])

# Make plot
vp_all_percent_melt <- melt(vp_all_percent, id.vars=c("SNP_Individual", "CellType", "convergence"), measure.vars = c("Tr", "PlateVar", "Control.ID", "Residual"))

vp_all_percent_melt$CellType <- factor(vp_all_percent_melt$CellType, levels=c("LCL", "IPSC", "CM"))
vp_all_percent_melt$variable <- factor(vp_all_percent_melt$variable, levels=c("Tr", "PlateVar", "Control.ID", "Residual"))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart_eachCell.anovaASE.5df.noDCM1R2H2O.box.noX.pdf", width=10, height=7)
	ggplot(vp_all_percent_melt, aes(x=variable, y=value, fill=CellType)) +  geom_boxplot(notch=T) + scale_fill_manual(values=c("#C77CFF", "#619CFF", "#00BA38")) + theme_bw() + ylab("Percent variance explained")
dev.off()

# Remove outliers
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart_eachCell.anovaASE.5df.noDCM1R2H2O.noOut.box.noX.pdf", width=10, height=7)
	ggplot(vp_all_percent_melt, aes(x=variable, y=value, fill=CellType)) +  geom_boxplot(outlier.shape=NA, notch=T) + scale_fill_manual(values=c("#C77CFF", "#619CFF", "#00BA38")) + theme_bw() + ylab("Percent variance explained")
dev.off()

# How much variance attributed to each category per cell type
vp_all_percent_melt %>% group_by(CellType, variable) %>% summarise(med = median(value))

###################
###################
# Now separate DE vs non DE for ASE variance
# Get DE genes
Deep_genes <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/All/DE_genes.txt", header=F, stringsAsFactors=F)
colnames(Deep_genes) <- c("enst", "DE_padj", "DE_pval", "DE_logFC", "ensg", "g.id", "type", "file.path")

# Put Deep_genes in a format which can be used to make plot
Deep_genes$file <- sapply(strsplit(Deep_genes$file.path, "/"), "[", 6) # Useful function that's similar to cut in unix
Deep_genes$CellType <- sapply(strsplit(Deep_genes$file, "_"), "[", 1) # Useful function that's similar to cut in unix
Deep_genes$Treatment.ID <- sapply(strsplit(Deep_genes$file, "_"), "[", 6) # Useful function that's similar to cut in unix
Deep_genes$Treatment.ID <- gsub(".txt", "", Deep_genes$Treatment.ID)
Deep_genes$CellType <- gsub("D", "", Deep_genes$CellType)

# Get gene info for SNPs
gene_info <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/Quasar_output/ASE_SNPs.genes.uniq.txt", header=F, stringsAsFactors=F)
colnames(gene_info) <- c("chr", "pos0", "pos1", "ref", "alt", "rsID", "ensg", "gene_class", "g.id")

# Get TATA box info
promoter_info <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/mixed_effect/hg38_tata_annotation_v3.txt", header=T, stringsAsFactors=F, sep="\t", fill=T)
TATA_box_genes <- unique(promoter_info[which(!promoter_info$TATA.Box.TBP..Promoter.Homer.Distance.From.Peak.sequence.strand.conservation. == ""),"Nearest.Ensembl"])

# Get CpG. I'm going to keep entries which have a "Nearest.Ensembl" listed and the smallest "Distance.to.TSS"
prom_CpG <- promoter_info[ promoter_info$Nearest.Ensembl != "", ]
prom_CpG <- prom_CpG[ order(prom_CpG$Nearest.Ensembl, abs(prom_CpG$Distance.to.TSS)), ]

prom_CpG_one <- prom_CpG[!duplicated(prom_CpG$Nearest.Ensembl),] # keeps nearest promoter

###
# Add gene info to SNPs (start with CMs)
vp_CM$rsID <- sapply(strsplit(rownames(vp_CM), "_"), "[", 1)
vp_CM_gene <- merge(vp_CM, gene_info, by="rsID")

# Try DE vs non-DE genes
Deep_genes_CM <- Deep_genes[which(Deep_genes$CellType == "CM"),]
vp_CM_gene$DE <- rep("Not_DE", length(vp_CM_gene$rsID))
vp_CM_gene$DE[vp_CM_gene$g.id %in% unique(Deep_genes_CM$g.id)] <- "DE"

# Add even more DE (logFC > 1)
Deep_genes_CM_1 <- Deep_genes_CM[ abs(Deep_genes_CM$DE_logFC) > 0.75, ]
vp_CM_gene$DE[vp_CM_gene$g.id %in% unique(Deep_genes_CM_1$g.id)] <- "More DE"

sum(vp_CM_gene$DE == "DE")
[1] 3598 # 3703 with X
sum(vp_CM_gene$DE == "More DE")
[1] 387 # Was 392 with X
sum(vp_CM_gene$DE == "Not_DE")
[1] 2896  # 3029 with X

vp_CM_gene$DE <- factor(vp_CM_gene$DE, levels=c("Not_DE", "DE", "More DE"))

# Add TATA info
vp_CM_gene$TATA <- rep("No_TATA", length(vp_CM_gene$ensg))
vp_CM_gene$TATA[which(vp_CM_gene$ensg %in% TATA_box_genes)] <- "TATA"

# Try CpG
vp_CM_gene <- merge(vp_CM_gene, prom_CpG_one[, c("Distance.to.TSS", "Nearest.Ensembl", "CpG.", "GC.")], by.x="ensg", by.y="Nearest.Ensembl", all.x=T, all.y=F)

###
# Now do LCLs
vp_LCL$rsID <- sapply(strsplit(rownames(vp_LCL), "_"), "[", 1)
vp_LCL_gene <- merge(vp_LCL, gene_info, by="rsID")

# Try DE vs non-DE genes
Deep_genes_LCL <- Deep_genes[which(Deep_genes$CellType == "LCL"),]
vp_LCL_gene$DE <- rep("Not_DE", length(vp_LCL_gene$rsID))
vp_LCL_gene$DE[vp_LCL_gene$g.id %in% unique(Deep_genes_LCL$g.id)] <- "DE"

# Add even more DE (logFC > 1)
Deep_genes_LCL_1 <- Deep_genes_LCL[ abs(Deep_genes_LCL$DE_logFC) > 0.75, ]
vp_LCL_gene$DE[vp_LCL_gene$g.id %in% unique(Deep_genes_LCL_1$g.id)] <- "More DE"

sum(vp_LCL_gene$DE == "DE")
[1] 5152 # 5800 with X
sum(vp_LCL_gene$DE == "More DE")
[1] 1978  # 2136 with X
sum(vp_LCL_gene$DE == "Not_DE")
[1] 382  # 434 with X

vp_LCL_gene$DE <- factor(vp_LCL_gene$DE, levels=c("Not_DE", "DE", "More DE"))

# Add TATA info
vp_LCL_gene$TATA <- rep("No_TATA", length(vp_LCL_gene$ensg))
vp_LCL_gene$TATA[which(vp_LCL_gene$ensg %in% TATA_box_genes)] <- "TATA"

# Try CpG
vp_LCL_gene <- merge(vp_LCL_gene, prom_CpG_one[, c("Distance.to.TSS", "Nearest.Ensembl", "CpG.", "GC.")], by.x="ensg", by.y="Nearest.Ensembl", all.x=T, all.y=F)

###
# Now do IPSCs
vp_IPSC$rsID <- sapply(strsplit(rownames(vp_IPSC), "_"), "[", 1)
vp_IPSC_gene <- merge(vp_IPSC, gene_info, by="rsID")

# Try DE vs non-DE genes
Deep_genes_IPSC <- Deep_genes[which(Deep_genes$CellType == "IPSC"),]
vp_IPSC_gene$DE <- rep("Not_DE", length(vp_IPSC_gene$rsID))
vp_IPSC_gene$DE[vp_IPSC_gene$g.id %in% unique(Deep_genes_IPSC$g.id)] <- "DE"

# Add even more DE (logFC > 1)
Deep_genes_IPSC_1 <- Deep_genes_IPSC[ abs(Deep_genes_IPSC$DE_logFC) > 0.75, ]
vp_IPSC_gene$DE[vp_IPSC_gene$g.id %in% unique(Deep_genes_IPSC_1$g.id)] <- "More DE"

sum(vp_IPSC_gene$DE == "DE")
[1] 4380  # was 4620
sum(vp_IPSC_gene$DE == "More DE")
[1] 1074  # was 1110 with X
sum(vp_IPSC_gene$DE == "Not_DE")
[1] 1575  # was 1634 with X

vp_IPSC_gene$DE <- factor(vp_IPSC_gene$DE, levels=c("Not_DE", "DE", "More DE"))

# Add TATA info
vp_IPSC_gene$TATA <- rep("No_TATA", length(vp_IPSC_gene$ensg))
vp_IPSC_gene$TATA[which(vp_IPSC_gene$ensg %in% TATA_box_genes)] <- "TATA"

# Try CpG
vp_IPSC_gene <- merge(vp_IPSC_gene, prom_CpG_one[, c("Distance.to.TSS", "Nearest.Ensembl", "CpG.", "GC.")], by.x="ensg", by.y="Nearest.Ensembl", all.x=T, all.y=F)

###################################
# Put all cell types in one graph. 

vp_allCellTypes_gene <- rbind(vp_LCL_gene, vp_IPSC_gene, vp_CM_gene)
vp_allCellTypes_gene_melt <- melt(vp_allCellTypes_gene, id.vars=c("SNP_Individual", "rsID", "CellType", "ensg", "g.id", "TATA", "DE", "CpG.", "GC."), measure.vars=c("Tr", "PlateVar", "Control.ID", "Residual"))

# TATA cells separate
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes.TATA.allVar.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt, aes(x=CellType, y=value, fill=TATA), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by TATA box") + xlab("Cell Type") + ylab("ASE Variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=13)) + facet_grid(. ~ variable) #+ ylim(-1,0)
dev.off()

# TATA cells together
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.TATA.allVar.box.noX.pdf", width=12)
ggplot(vp_allCellTypes_gene_melt, aes(x=TATA, y=value, fill=TATA), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by TATA box") + xlab("Group") + ylab("ASE variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + facet_grid(. ~ variable) #+ ylim(0,0.7)
dev.off()

# TATA cells together, only plotting residual variance
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.TATA.ResidualVar.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt %>% filter(variable == "Residual"), aes(x=TATA, y=value, fill=TATA), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by TATA box") + xlab("Group") + ylab("ASE variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) #+ ylim(0,0.7)
dev.off()


# DE cells separate
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes.DE.allVar.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt, aes(x=CellType, y=value, fill=DE), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by DE status") + xlab("Cell Type") + ylab("ASE Variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=13)) + facet_grid(. ~ variable) #+ ylim(-1,0)
dev.off()

# DE cells together
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.DE.allVar.box.noX.pdf", width=12)
ggplot(vp_allCellTypes_gene_melt, aes(x=DE, y=value, fill=DE), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by DE status") + xlab("Group") + ylab("ASE variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + facet_grid(. ~ variable) #+ ylim(0,0.7)
dev.off()

# For DE, plot just Residual and treatment variance since those 2 are significant
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.DE.Resid_Tr.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt %>% filter(variable == "Tr" | variable == "Residual"), aes(x=DE, y=value, fill=DE), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by DE status") + xlab("Group") + ylab("ASE variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + facet_grid(. ~ variable) + ylim(0,0.6)
dev.off()


# GC and CpG graphs
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.GC.allVar.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_melt, aes(x=as.numeric(GC.), y=log(value + 0.0001), color=CellType)) + geom_point() + theme_bw() + labs(title="ASE variance vs GC Content") + xlab("GC") + ylab("log(Variance)") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + facet_grid(variable ~ .) + scale_color_manual(values=c("#00BA38", "#619CFF", "#C77CFF")) #+ ylim(-12,3)
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.CpG.allVar.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_melt, aes(x=as.numeric(CpG.), y=log(value + 0.0001), color=CellType)) + geom_point() + theme_bw() + labs(title="ASE variance vs CpG") + xlab("CpG") + ylab("log(Variance)") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + facet_grid(variable ~ .) + scale_color_manual(values=c("#00BA38", "#619CFF", "#C77CFF")) #+ ylim(-12,3)
dev.off()

##########
# In order to see if DE and TATA are significant, I need to correct for ASE sequencing depth
# Add number of reads for ASE SNP

# First I need to get the allelic expression for all SNPs tested
LCL_counts <- data.frame(by_SNP.Ind_LCL %>% summarise(ref=sum(ref.reads1), alt=sum(alt.reads1)))
LCL_counts$total <- LCL_counts$ref + LCL_counts$alt
LCL_counts$rankTotal <- rank(LCL_counts$total)

#IPSCs
IPSC_counts <- data.frame(by_SNP.Ind_IPSC %>% summarise(ref=sum(ref.reads1), alt=sum(alt.reads1)))
IPSC_counts$total <- IPSC_counts$ref + IPSC_counts$alt
IPSC_counts$rankTotal <- rank(IPSC_counts$total)

#CMs
CM_counts <- data.frame(by_SNP.Ind_CM %>% summarise(ref=sum(ref.reads1), alt=sum(alt.reads1)))
CM_counts$total <- CM_counts$ref + CM_counts$alt
CM_counts$rankTotal <- rank(CM_counts$total)

# Now put them all together
LCL_counts$CellType <- "LCL"
IPSC_counts$CellType <- "IPSC"
CM_counts$CellType <- "CM"

all_counts <- rbind(LCL_counts, IPSC_counts, CM_counts)

vp_allCellTypes_gene_melt <- merge(vp_allCellTypes_gene_melt, all_counts, by=c("SNP_Individual", "CellType"))

# Test DE with SNP rank  for residual variance (also sig when I use total ASE reads, not just rank)
summary(lm(value ~ as.integer(DE) + rankTotal, data=vp_allCellTypes_gene_melt[ vp_allCellTypes_gene_melt$variable == "Residual",]))

Residuals:
     Min       1Q   Median       3Q      Max
-0.20100 -0.03489 -0.00354  0.02002  0.47447

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
(Intercept)     1.855e-01  1.809e-03  102.55   <2e-16 ***
as.integer(DE)  9.943e-03  7.509e-04   13.24   <2e-16 ***
rankTotal      -3.140e-05  2.357e-07 -133.22   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06771 on 21419 degrees of freedom
Multiple R-squared:  0.4621,    Adjusted R-squared:  0.4621
F-statistic:  9201 on 2 and 21419 DF,  p-value: < 2.2e-16

# Test DE with SNP rank for treatment variance (also sig when I use total ASE reads, not just rank)
summary(lm(value ~ as.integer(DE) + rankTotal, data=vp_allCellTypes_gene_melt[ vp_allCellTypes_gene_melt$variable == "Tr",]))

Residuals:
     Min       1Q   Median       3Q      Max
-0.10486 -0.03328 -0.00346  0.01715  1.02978

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
(Intercept)     1.025e-01  1.713e-03  59.835   <2e-16 ***
as.integer(DE)  1.054e-03  7.113e-04   1.481    0.139   # This was sig with X
rankTotal      -1.808e-05  2.233e-07 -80.990   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06413 on 21419 degrees of freedom
Multiple R-squared:  0.2366,    Adjusted R-squared:  0.2366
F-statistic:  3320 on 2 and 21419 DF,  p-value: < 2.2e-16

# What if I get rid of "More DE" category?
vp_allCellTypes_gene_melt_test <- vp_allCellTypes_gene_melt
vp_allCellTypes_gene_melt_test$DE <- as.character(vp_allCellTypes_gene_melt_test$DE)
vp_allCellTypes_gene_melt_test[vp_allCellTypes_gene_melt_test$DE == "More DE", "DE"]) <- "DE"

vp_allCellTypes_gene_melt_test$DE <- factor(vp_allCellTypes_gene_melt_test$DE, levels = c("Not_DE", "DE"))

# Test TATA with SNP rank  for residual variance (also sig when I use total ASE reads, not just rank)
summary(lm(value ~ TATA + rankTotal, data=vp_allCellTypes_gene_melt[ vp_allCellTypes_gene_melt$variable == "Residual",]))

Residuals:
     Min       1Q   Median       3Q      Max
-0.19656 -0.03505 -0.00356  0.01976  0.47482

Coefficients:
              Estimate Std. Error  t value Pr(>|t|)
(Intercept)  2.052e-01  9.959e-04  206.044   <2e-16 ***
TATATATA     2.403e-03  1.024e-03    2.347    0.019 *
rankTotal   -3.173e-05  2.364e-07 -134.249   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06797 on 21419 degrees of freedom
Multiple R-squared:  0.4579,    Adjusted R-squared:  0.4578
F-statistic:  9044 on 2 and 21419 DF,  p-value: < 2.2e-16

# Test TATA with SNP rank for treatment variance
summary(lm(value ~ TATA + rankTotal, data=vp_allCellTypes_gene_melt[ vp_allCellTypes_gene_melt$variable == "Tr",]))

Residuals:
     Min       1Q   Median       3Q      Max
-0.10570 -0.03328 -0.00364  0.01712  1.03010

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.043e-01  9.396e-04  111.02   <2e-16 ***
TATATATA     1.478e-03  9.660e-04    1.53    0.126
rankTotal   -1.814e-05  2.230e-07  -81.34   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06413 on 21419 degrees of freedom
Multiple R-squared:  0.2366,    Adjusted R-squared:  0.2366
F-statistic:  3320 on 2 and 21419 DF,  p-value: < 2.2e-16

# Do CpG and GC content
summary(lm(value ~ rankTotal + as.numeric(GC.), data=vp_allCellTypes_gene_melt_resid %>% filter(variable == "Tr")))

Residuals:
     Min       1Q   Median       3Q      Max
-0.10579 -0.03034 -0.00272  0.01609  1.02896

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)      9.904e-02  2.577e-03  38.425   <2e-16 ***
rankTotal       -1.795e-05  2.390e-07 -75.090   <2e-16 ***
as.numeric(GC.)  8.444e-03  3.743e-03   2.256   0.0241 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06109 on 17162 degrees of freedom
  (4257 observations deleted due to missingness)
Multiple R-squared:  0.2476,    Adjusted R-squared:  0.2475
F-statistic:  2823 on 2 and 17162 DF,  p-value: < 2.2e-16

summary(lm(value ~ rankTotal + as.numeric(GC.), data=vp_allCellTypes_gene_melt_resid %>% filter(variable == "Residual")))

Residuals:
     Min       1Q   Median       3Q      Max
-0.19118 -0.03171 -0.00262  0.01828  0.47988

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)      2.156e-01  2.695e-03   79.99  < 2e-16 ***
rankTotal       -3.080e-05  2.500e-07 -123.23  < 2e-16 ***
as.numeric(GC.) -2.340e-02  3.914e-03   -5.98 2.28e-09 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06388 on 17162 degrees of freedom
  (4257 observations deleted due to missingness)
Multiple R-squared:  0.4728,    Adjusted R-squared:  0.4727
F-statistic:  7696 on 2 and 17162 DF,  p-value: < 2.2e-16

summary(lm(value ~ rankTotal + as.numeric(CpG.), data=vp_allCellTypes_gene_melt_resid %>% filter(variable == "Tr")))

Residuals:
     Min       1Q   Median       3Q      Max
-0.10579 -0.03034 -0.00270  0.01614  1.02916

Coefficients:
                   Estimate Std. Error t value Pr(>|t|)
(Intercept)       1.029e-01  1.217e-03  84.611   <2e-16 ***
rankTotal        -1.795e-05  2.396e-07 -74.922   <2e-16 ***
as.numeric(CpG.)  1.784e-02  8.729e-03   2.044   0.0409 *
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0611 on 17166 degrees of freedom
  (4253 observations deleted due to missingness)
Multiple R-squared:  0.2472,    Adjusted R-squared:  0.2471
F-statistic:  2818 on 2 and 17166 DF,  p-value: < 2.2e-16

summary(lm(value ~ rankTotal + as.numeric(CpG.), data=vp_allCellTypes_gene_melt_resid %>% filter(variable == "Residual")))

Residuals:
     Min       1Q   Median       3Q      Max
-0.18940 -0.03186 -0.00273  0.01829  0.47918

Coefficients:
                   Estimate Std. Error t value Pr(>|t|)
(Intercept)       2.061e-01  1.272e-03  162.07  < 2e-16 ***
rankTotal        -3.074e-05  2.504e-07 -122.77  < 2e-16 ***
as.numeric(CpG.) -6.815e-02  9.123e-03   -7.47 8.38e-14 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06386 on 17166 degrees of freedom
  (4253 observations deleted due to missingness)
Multiple R-squared:  0.4732,    Adjusted R-squared:  0.4732
F-statistic:  7711 on 2 and 17166 DF,  p-value: < 2.2e-16


## Compare genes variance based on loss of function intolerance
# Load Gnomad data in:
gnomad <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/annotation/gnomad.v2.1.1.lof_metrics.by_gene.txt", header=T, stringsAsFactors=F, sep="\t")

# I think the columns I'm interested in are "oe_lof", "oe_lof_lower", and "oe_lof_upper"
vp_allCellTypes_gene_melt_lof <- merge(vp_allCellTypes_gene_melt, gnomad[,c("gene", "oe_lof", "oe_lof_lower", "oe_lof_upper")], by.x = "g.id", by.y="gene") # 4891/6310 genes had lof score

summary(lm(value ~ oe_lof + rankTotal, data=vp_allCellTypes_gene_melt_lof[ vp_allCellTypes_gene_melt_lof$variable == "Residual",]))

Residuals:
     Min       1Q   Median       3Q      Max
-0.18667 -0.03169 -0.00270  0.01824  0.47437

Coefficients:
              Estimate Std. Error  t value Pr(>|t|)
(Intercept)  1.953e-01  1.351e-03  144.587  < 2e-16 ***
oe_lof       8.759e-03  1.455e-03    6.019  1.8e-09 ***
rankTotal   -3.064e-05  2.548e-07 -120.250  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.0632 on 16430 degrees of freedom
  (67 observations deleted due to missingness)
Multiple R-squared:  0.4752,    Adjusted R-squared:  0.4752
F-statistic:  7440 on 2 and 16430 DF,  p-value: < 2.2e-16

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.LoF.allVar.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_lof, aes(x=as.numeric(oe_lof), y=log(value + 0.0001), color=CellType)) + geom_point() + theme_bw() + labs(title="ASE variance vs LoF OE") + xlab("LoF OE") + ylab("log(Variance)") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + facet_grid(variable ~ .) + scale_color_manual(values=c("#00BA38", "#619CFF", "#C77CFF")) #+ ylim(-12,3)
dev.off()

# Just residual variance
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.LoF.ResidVar.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_lof %>% filter(variable == "Residual"), aes(x=as.numeric(oe_lof), y=value, color=CellType)) + geom_point() + theme_bw() + labs(title="ASE variance vs LoF OE") + xlab("LoF OE") + ylab("Residual Variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + scale_color_manual(values=c("#00BA38", "#619CFF", "#C77CFF")) #+ ylim(-12,3)
dev.off()

# Those scatter plots didn't look good, so take residual variance and break into deciles.
vp_allCellTypes_gene_melt_lof %>% filter((variable == "Residual" | variable == "Tr") & !is.na(oe_lof)) -> vp_allCellTypes_gene_melt_lof_resid
vp_allCellTypes_gene_melt_lof_resid$decile <- cut2(vp_allCellTypes_gene_melt_lof_resid$oe_lof, g=10)

# I think I want to regress out the rankTotal so I can then make the plot
vp_allCellTypes_gene_melt_lof_resid$regress <- resid(lm(value ~ rankTotal, data = vp_allCellTypes_gene_melt_lof_resid))


pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.LoF.Resid.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_lof_resid, aes(x=decile, y=regress, fill=decile), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by LoF Decile") + xlab("LoF Decile") + ylab("ASE Residual Variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) # + ylim(0,0.6)
dev.off()

# Try rank(Residual variance)
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.LoF.RankResid.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_lof_resid, aes(x=decile, y=rank(regress), fill=decile), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by LoF Decile") + xlab("LoF Decile") + ylab("ASE Residual Variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) # + ylim(0,0.6)
dev.off()


# Try scatter
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.LoF.ResidVar.regress.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_lof_resid, aes(x=oe_lof, y=regress, color=CellType)) + geom_point() + theme_bw() + labs(title="ASE variance vs LoF OE") + xlab("LoF OE") + ylab("Residual Variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + scale_color_manual(values=c("#00BA38", "#619CFF", "#C77CFF")) #+ ylim(-12,3)
dev.off()

# Plot rank(Reisudal variance)
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.LoF.RankResidVar.regress.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_lof_resid, aes(x=oe_lof, y=rank(regress), color=CellType)) + geom_point() + theme_bw() + labs(title="ASE variance vs LoF OE") + xlab("LoF OE") + ylab("Residual Variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) + scale_color_manual(values=c("#00BA38", "#619CFF", "#C77CFF")) #+ ylim(-12,3)
dev.off()



### Do the same thing for TATA and DE genes, regressing out rankTotal before plotting
# I think I want to regress out the rankTotal so I can then make the plot
vp_allCellTypes_gene_melt %>% filter(variable == "Residual" | variable == "Tr") -> vp_allCellTypes_gene_melt_resid
vp_allCellTypes_gene_melt_resid$regress <- resid(lm(value ~ rankTotal, data = vp_allCellTypes_gene_melt_resid))

# TATA cells together, only plotting residual variance after removing rankTotal
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.TATA.ResidualVar.regress.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_resid, aes(x=TATA, y=regress, fill=TATA), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by TATA box") + xlab("Group") + ylab("ASE variance after regressing rankTotal") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) #+ ylim(0,0.7)
dev.off()

# Remove outliers
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.TATA.ResidualVar.regress.noOut.box.noX.pdf", width = 3.8)
ggplot(vp_allCellTypes_gene_melt_resid, aes(x=TATA, y=regress, fill=TATA), color="black") + geom_boxplot(outlier.shape=NA, notch=TRUE) + theme_bw() + labs(title="ASE Variance by TATA box") + xlab("Group") + ylab("ASE variance after regressing rankTotal") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + ylim(-0.15,0.15)
dev.off()

# For DE, plot just Residual and treatment variance since those 2 are significant
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.DE.Resid_Tr.regress.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_resid %>% filter(variable == "Tr" | variable == "Residual"), aes(x=DE, y=regress, fill=DE), color="black") + geom_boxplot(notch=TRUE) + theme_bw() + labs(title="ASE Variance by DE status") + xlab("Group") + ylab("ASE variance after regressing rankTotal") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + facet_grid(. ~ variable) #+ ylim(0,0.6)
dev.off()

# Remove outliers
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.DE.Resid_Tr.regress.noOut.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_resid %>% filter(variable == "Tr" | variable == "Residual"), aes(x=DE, y=regress, fill=DE), color="black") + geom_boxplot(outlier.shape = NA, notch=TRUE) + theme_bw() + labs(title="ASE Variance by DE status") + xlab("Group") + ylab("ASE variance after regressing rankTotal") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + facet_grid(. ~ variable) + ylim(0,0.2)
dev.off()



######################################
######################################
######################################
######
# I want to compare residual variance for the same SNP across cell types. Consider only SNPs with ASE by ANOVA in both.

vp_LCL <- merge(vp_LCL, LCL_counts, by=c("SNP_Individual", "CellType"))
vp_IPSC <- merge(vp_IPSC, IPSC_counts, by=c("SNP_Individual", "CellType"))
LCL_IPSC_var <- merge(vp_LCL, vp_IPSC, by=c("SNP_Individual", "rsID")) # LCL is .x
rcorr(LCL_IPSC_var$Residual.x, LCL_IPSC_var$Residual.y, type="spearman")
     x    y
x 1.00 0.59  # was 0.58
y 0.59 1.00

n= 1160 # Was 1265


P
  x  y
x     0
y  0

# Roger wants a 95% CI:
tanh(atanh(0.59) + 1.96/(sqrt(1160-3)))
[1] 0.62629
 tanh(atanh(0.59) - 1.96/(sqrt(1160-3)))
[1] 0.5511586

# Check if it's significant after accounting for expression
summary(lm(formula = Residual.x ~ Residual.y + rankTotal.x + rankTotal.y, data = LCL_IPSC_var))

Residuals:
     Min       1Q   Median       3Q      Max
-0.25627 -0.03179 -0.00221  0.01872  0.37024

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  2.162e-01  9.270e-03  23.319  < 2e-16 ***
Residual.y   1.664e-01  3.121e-02   5.332 1.17e-07 ***
rankTotal.x -3.819e-05  1.400e-06 -27.267  < 2e-16 ***
rankTotal.y  7.030e-06  1.814e-06   3.876 0.000112 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06461 on 1156 degrees of freedom
Multiple R-squared:  0.5432,    Adjusted R-squared:  0.542
F-statistic: 458.2 on 3 and 1156 DF,  p-value: < 2.2e-16



# LCLs vs CMs
vp_CM <- merge(vp_CM, CM_counts, by=c("SNP_Individual", "CellType"))

LCL_CM_var <- merge(vp_LCL, vp_CM, by=c("SNP_Individual", "rsID")) # LCL is .x
rcorr(LCL_CM_var$Residual.x, LCL_CM_var$Residual.y, type="spearman")
     x    y
x 1.00 0.55  # Was 0.53
y 0.55 1.00

n= 963 # was 1069


P
  x  y
x     0
y  0

# Roger wants a 95% CI:
tanh(atanh(0.55) + 1.96/(sqrt(963-3)))
[1] 0.5925846
 tanh(atanh(0.55) - 1.96/(sqrt(963-3)))
[1] 0.5043496

# Check if it's significant after accounting for expression
summary(lm(formula = Residual.x ~ Residual.y + rankTotal.x + rankTotal.y, data = LCL_CM_var))

Residuals:
     Min       1Q   Median       3Q      Max
-0.19769 -0.02914 -0.00316  0.01892  0.33214

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.765e-01  9.469e-03  18.638  < 2e-16 ***
Residual.y   3.495e-01  3.439e-02  10.161  < 2e-16 ***
rankTotal.x -3.767e-05  1.342e-06 -28.065  < 2e-16 ***
rankTotal.y  1.311e-05  1.897e-06   6.907 8.98e-12 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06422 on 959 degrees of freedom
Multiple R-squared:  0.5777,    Adjusted R-squared:  0.5764
F-statistic: 437.3 on 3 and 959 DF,  p-value: < 2.2e-16

# IPSCs vs CMs
IPSC_CM_var <- merge(vp_IPSC, vp_CM, by=c("SNP_Individual", "rsID")) # IPSC is .x
rcorr(IPSC_CM_var$Residual.x, IPSC_CM_var$Residual.y, type="spearman")
     x    y
x 1.00 0.58  # was 0.57
y 0.58 1.00

n= 1692 # Was 1801


P
  x  y
x     0
y  0

# Roger wants a 95% CI:
tanh(atanh(0.58) + 1.96/(sqrt(1692-3)))
[1] 0.6107735
 tanh(atanh(0.58) - 1.96/(sqrt(1692-3)))
[1] 0.5474769

# Check if it's significant after accounting for expression
summary(lm(formula = Residual.x ~ Residual.y + rankTotal.x + rankTotal.y, data = IPSC_CM_var))

Residuals:
     Min       1Q   Median       3Q      Max
-0.19235 -0.02876 -0.00108  0.01524  0.42127

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.880e-01  7.169e-03  26.230  < 2e-16 ***
Residual.y   2.195e-01  2.568e-02   8.549  < 2e-16 ***
rankTotal.x -3.648e-05  1.076e-06 -33.912  < 2e-16 ***
rankTotal.y  7.633e-06  1.428e-06   5.344 1.04e-07 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05893 on 1688 degrees of freedom
Multiple R-squared:  0.5681,    Adjusted R-squared:  0.5674
F-statistic: 740.2 on 3 and 1688 DF,  p-value: < 2.2e-16


# Now put all ASE only on the same plot
LCL_IPSC_var$Group <- rep("LCL_IPSC")
LCL_CM_var$Group <- rep("LCL_CM")
IPSC_CM_var$Group <- rep("IPSC_CM")

All_comp_var <- rbind(LCL_IPSC_var, LCL_CM_var, IPSC_CM_var)

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/residVar_betwen_cellTypes.scatter.noX.pdf")
p <- ggplot(All_comp_var, aes(log(Residual.x), log(Residual.y), fill=Group))
print(p + geom_point(size=3, pch=21) + labs(title = "Residual Variance Across Cell Types", x = "log(residual variance)", y = "log(residual variance)") + theme_bw() + theme(text = element_text(size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + geom_abline(intercept=0, slope=1))
dev.off()

########################
########################
# Calculate overlap between cASE genes and DEGs and DSGs
# Get DE genes
# Keep only the transcrpt with the greatest lfc per gene per treatment per cell type
all_genes <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/All/all_DEgene_results.txt", header=F, stringsAsFactors=F)
colnames(all_genes) <- c("enst", "DE_padj", "DE_pval", "DE_logFC", "ensg", "g.id", "type", "file.path")

# Remove all of the header lines spread throughout the file
all_genes <- all_genes[ all_genes$DE_logFC != "logFC", ]
all_genes$DE_padj <- as.numeric(all_genes$DE_padj)
all_genes$DE_pval <- as.numeric(all_genes$DE_pval)
all_genes$DE_logFC <- as.numeric(all_genes$DE_logFC)

# Keep only the transcrpt with the greatest lfc per gene per treatment per cell type
all_genes <- all_genes %>% arrange(ensg, file.path, desc(abs(DE_logFC)))

all_genes_one <- all_genes[ !duplicated(all_genes[,c("ensg", "file.path")]), ]

# Put all_genes_one in a format which can be used to make plot
all_genes_one$file <- sapply(strsplit(all_genes_one$file.path, "/"), "[", 6) # Useful function that's similar to cut in unix
all_genes_one$CellType <- sapply(strsplit(all_genes_one$file, "_"), "[", 1) # Useful function that's similar to cut in unix
all_genes_one$Treatment.ID <- sapply(strsplit(all_genes_one$file, "_"), "[", 6) # Useful function that's similar to cut in unix

all_genes_one$Treatment.ID <- gsub(".txt", "", all_genes_one$Treatment.ID)
all_genes_one$CellType <- gsub("D", "", all_genes_one$CellType)


# Add gene info to SNPs (start with CMs)
mylm_CM_5df_Tr_anovaASE$rsID <- sapply(strsplit(mylm_CM_5df_Tr_anovaASE$SNP_Individual, "_"), "[", 1)

mylm_CM_5df_Tr_anovaASE_gene <- merge(mylm_CM_5df_Tr_anovaASE, gene_info, by="rsID")

# Add T* treatment ID's to ASE, since that's what I use for genes
treatKey <- data.frame(Treatment.ID = c("T42C1", "T9C1", "T27C1", "T13C1", "T15C1", "T12C1", "T33C1", "T14C1", "T19C1", "T30C1", "T6C1", "T20C1"), term = c("TrAcetaminophin", "TrAldosterone", "TrCadmium", "TrCaffeine", "TrCopper", "TrDexamethasone", "TrInsulin", "TrNicotine", "TrSelenium", "TrTriclosan", "TrVitaminA", "TrZinc"))

mylm_CM_5df_Tr_anovaASE_gene <- merge(mylm_CM_5df_Tr_anovaASE_gene, treatKey, by="term")

CM_ASE_GE <- merge(mylm_CM_5df_Tr_anovaASE_gene, all_genes_one[, c("DE_padj", "DE_pval", "DE_logFC", "ensg", "CellType", "Treatment.ID")], by=c("CellType", "Treatment.ID", "ensg"))

# See if cASE is enriched in DE genes
# cASE(+) DE(+)
cASE_P <- sum(CM_ASE_GE$padj < 0.1)
cASE_N <- sum(CM_ASE_GE$padj >= 0.1)

cASE_P_DE_P <- sum(abs(CM_ASE_GE$DE_logFC) >= 0.25 & CM_ASE_GE$DE_padj < 0.1 & CM_ASE_GE$padj < 0.1)
cASE_P_DE_N <- cASE_P - cASE_P_DE_P

cASE_N_DE_P <- sum(abs(CM_ASE_GE$DE_logFC) >= 0.25 & CM_ASE_GE$DE_padj < 0.1 & CM_ASE_GE$padj >= 0.1)
cASE_N_DE_N <- cASE_N - cASE_N_DE_P

# Do fisher's test to see if cASE is enriched for DE
CM_enrich <-
matrix(c(cASE_P_DE_P, cASE_N_DE_P, cASE_P_DE_N, cASE_N_DE_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DE = c("DE+", "DE-")))
fish_test_CM <- fisher.test(CM_enrich, alternative = "two.sided")

# LCL DEGs
mylm_LCL_5df_Tr_anovaASE$rsID <- sapply(strsplit(mylm_LCL_5df_Tr_anovaASE$SNP_Individual, "_"), "[", 1)

mylm_LCL_5df_Tr_anovaASE_gene <- merge(mylm_LCL_5df_Tr_anovaASE, gene_info, by="rsID")

mylm_LCL_5df_Tr_anovaASE_gene <- merge(mylm_LCL_5df_Tr_anovaASE_gene, treatKey, by="term")

LCL_ASE_GE <- merge(mylm_LCL_5df_Tr_anovaASE_gene, all_genes_one[, c("DE_padj", "DE_pval", "DE_logFC", "ensg", "CellType", "Treatment.ID")], by=c("CellType", "Treatment.ID", "ensg"))

# See if cASE is enriched in DE genes
# cASE(+) DE(+)
cASE_P <- sum(LCL_ASE_GE$padj < 0.1)
cASE_N <- sum(LCL_ASE_GE$padj >= 0.1)

cASE_P_DE_P <- sum(abs(LCL_ASE_GE$DE_logFC) >= 0.25 & LCL_ASE_GE$DE_padj < 0.1 & LCL_ASE_GE$padj < 0.1)
cASE_P_DE_N <- cASE_P - cASE_P_DE_P

cASE_N_DE_P <- sum(abs(LCL_ASE_GE$DE_logFC) >= 0.25 & LCL_ASE_GE$DE_padj < 0.1 & LCL_ASE_GE$padj >= 0.1)
cASE_N_DE_N <- cASE_N - cASE_N_DE_P

# Do fisher's test to see if cASE is enriched for DE
LCL_enrich <-
matrix(c(cASE_P_DE_P, cASE_N_DE_P, cASE_P_DE_N, cASE_N_DE_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DE = c("DE+", "DE-")))
fish_test_LCL <- fisher.test(LCL_enrich, alternative = "two.sided")

# IPSC DEGs
mylm_IPSC_5df_Tr_anovaASE$rsID <- sapply(strsplit(mylm_IPSC_5df_Tr_anovaASE$SNP_Individual, "_"), "[", 1)

mylm_IPSC_5df_Tr_anovaASE_gene <- merge(mylm_IPSC_5df_Tr_anovaASE, gene_info, by="rsID")

mylm_IPSC_5df_Tr_anovaASE_gene <- merge(mylm_IPSC_5df_Tr_anovaASE_gene, treatKey, by="term")

IPSC_ASE_GE <- merge(mylm_IPSC_5df_Tr_anovaASE_gene, all_genes_one[, c("DE_padj", "DE_pval", "DE_logFC", "ensg", "CellType", "Treatment.ID")], by=c("CellType", "Treatment.ID", "ensg"))

# See if cASE is enriched in DE genes
# cASE(+) DE(+)
cASE_P <- sum(IPSC_ASE_GE$padj < 0.1)
cASE_N <- sum(IPSC_ASE_GE$padj >= 0.1)

cASE_P_DE_P <- sum(abs(IPSC_ASE_GE$DE_logFC) >= 0.25 & IPSC_ASE_GE$DE_padj < 0.1 & IPSC_ASE_GE$padj < 0.1)
cASE_P_DE_N <- cASE_P - cASE_P_DE_P

cASE_N_DE_P <- sum(abs(IPSC_ASE_GE$DE_logFC) >= 0.25 & IPSC_ASE_GE$DE_padj < 0.1 & IPSC_ASE_GE$padj >= 0.1)
cASE_N_DE_N <- cASE_N - cASE_N_DE_P

# Do fisher's test to see if cASE is enriched for DE
IPSC_enrich <-
matrix(c(cASE_P_DE_P, cASE_N_DE_P, cASE_P_DE_N, cASE_N_DE_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DE = c("DE+", "DE-")))
fish_test_IPSC <- fisher.test(IPSC_enrich, alternative = "two.sided")

### Now try differentially spliced genes. These files only have the differentially spliced genes
CM_DSG <- read.table("../annotation/CM_DSG.txt", stringsAsFactors=F, header=F)
colnames(CM_DSG) <- c("ensg", "CellType", "Treatment.ID")
CM_DSG$DSG <- 1

CM_ASE_DSG <- merge(mylm_CM_5df_Tr_anovaASE_gene, CM_DSG, by=c("CellType", "Treatment.ID", "ensg"), all.x=T)

CM_ASE_DSG$DSG[is.na(CM_ASE_DSG$DSG)] <- 0 # Non-DSG genes are 0

# See if cASE is enriched in DE genes
# cASE(+) DE(+)
cASE_P <- sum(CM_ASE_DSG$padj < 0.1)
cASE_N <- sum(CM_ASE_DSG$padj >= 0.1)

cASE_P_DSG_P <- sum(CM_ASE_DSG$DSG == 1 & CM_ASE_DSG$padj < 0.1)
cASE_P_DSG_N <- cASE_P - cASE_P_DSG_P

cASE_N_DSG_P <- sum(CM_ASE_DSG$DSG == 1 & CM_ASE_DSG$padj >= 0.1)
cASE_N_DSG_N <- cASE_N - cASE_N_DSG_P

# Do fisher's test to see if cASE is enriched for DSG
CM_enrich_DSG <-
matrix(c(cASE_P_DSG_P, cASE_N_DSG_P, cASE_P_DSG_N, cASE_N_DSG_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DSG = c("DSG+", "DSG-")))
fish_test_CM_DSG <- fisher.test(CM_enrich_DSG, alternative = "two.sided")

# Now do LCL
LCL_DSG <- read.table("../annotation/LCL_DSG.txt", stringsAsFactors=F, header=F)
colnames(LCL_DSG) <- c("ensg", "CellType", "Treatment.ID")
LCL_DSG$DSG <- 1

LCL_ASE_DSG <- merge(mylm_LCL_5df_Tr_anovaASE_gene, LCL_DSG, by=c("CellType", "Treatment.ID", "ensg"), all.x=T)

LCL_ASE_DSG$DSG[is.na(LCL_ASE_DSG$DSG)] <- 0 # Non-DSG genes are 0

# See if cASE is enriched in DE genes
# cASE(+) DE(+)
cASE_P <- sum(LCL_ASE_DSG$padj < 0.1)
cASE_N <- sum(LCL_ASE_DSG$padj >= 0.1)

cASE_P_DSG_P <- sum(LCL_ASE_DSG$DSG == 1 & LCL_ASE_DSG$padj < 0.1)
cASE_P_DSG_N <- cASE_P - cASE_P_DSG_P

cASE_N_DSG_P <- sum(LCL_ASE_DSG$DSG == 1 & LCL_ASE_DSG$padj >= 0.1)
cASE_N_DSG_N <- cASE_N - cASE_N_DSG_P

# Do fisher's test to see if cASE is enriched for DSG
LCL_enrich_DSG <-
matrix(c(cASE_P_DSG_P, cASE_N_DSG_P, cASE_P_DSG_N, cASE_N_DSG_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DSG = c("DSG+", "DSG-")))
fish_test_LCL_DSG <- fisher.test(LCL_enrich_DSG, alternative = "two.sided")

# Now do IPSC
IPSC_DSG <- read.table("../annotation/IPSC_DSG.txt", stringsAsFactors=F, header=F)
colnames(IPSC_DSG) <- c("ensg", "CellType", "Treatment.ID")
IPSC_DSG$DSG <- 1

IPSC_ASE_DSG <- merge(mylm_IPSC_5df_Tr_anovaASE_gene, IPSC_DSG, by=c("CellType", "Treatment.ID", "ensg"), all.x=T)

IPSC_ASE_DSG$DSG[is.na(IPSC_ASE_DSG$DSG)] <- 0 # Non-DSG genes are 0

# See if cASE is enriched in DE genes
# cASE(+) DE(+)
cASE_P <- sum(IPSC_ASE_DSG$padj < 0.1)
cASE_N <- sum(IPSC_ASE_DSG$padj >= 0.1)

cASE_P_DSG_P <- sum(IPSC_ASE_DSG$DSG == 1 & IPSC_ASE_DSG$padj < 0.1)
cASE_P_DSG_N <- cASE_P - cASE_P_DSG_P

cASE_N_DSG_P <- sum(IPSC_ASE_DSG$DSG == 1 & IPSC_ASE_DSG$padj >= 0.1)
cASE_N_DSG_N <- cASE_N - cASE_N_DSG_P

# Do fisher's test to see if cASE is enriched for DSG
IPSC_enrich_DSG <-
matrix(c(cASE_P_DSG_P, cASE_N_DSG_P, cASE_P_DSG_N, cASE_N_DSG_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DSG = c("DSG+", "DSG-")))
fish_test_IPSC_DSG <- fisher.test(IPSC_enrich_DSG, alternative = "two.sided")

# Put data together and make a forest plot
#for_plot <- data.frame(CellType = c("LCL", "IPSC", "CM", "LCL", "IPSC", "CM"), AnalysisType = c("Expression", "Expression", "Expression", "Splicing", "Splicing", "Splicing"), OR = c(fish_test_LCL$estimate, fish_test_IPSC$estimate, fish_test_CM$estimate, fish_test_LCL_DSG$estimate, fish_test_IPSC_DSG$estimate, fish_test_CM_DSG$estimate), Low = c(fish_test_LCL$conf.int[1], fish_test_IPSC$conf.int[1], fish_test_CM$conf.int[1], fish_test_LCL_DSG$conf.int[1], fish_test_IPSC_DSG$conf.int[1], fish_test_CM_DSG$conf.int[1]), Hi = c(fish_test_LCL$conf.int[2], fish_test_IPSC$conf.int[2], fish_test_CM$conf.int[2], fish_test_LCL_DSG$conf.int[2], fish_test_IPSC_DSG$conf.int[2], fish_test_CM_DSG$conf.int[2]), pval = c(fish_test_LCL$p.value, fish_test_IPSC$p.value, fish_test_CM$p.value, fish_test_LCL_DSG$p.value, fish_test_IPSC_DSG$p.value, fish_test_CM_DSG$p.value))

#for_plot$sampleID <- paste0(for_plot$CellType, "_", for_plot$AnalysisType)

# Plot
#fp <- ggplot(data=for_plot, aes(x=sampleID, y=OR, ymin=Low, ymax=Hi)) +
#geom_pointrange(aes(color=AnalysisType), size=0.8) +
#geom_hline(yintercept=1, lty=2) +
#coord_flip() +  # flip coordinates (puts labels on y axis)
#xlab("Cell Type") + ylab("Enrichment (95% CI)") + theme_bw() + geom_errorbar(aes(ymin=Low, ymax=Hi, col=AnalysisType), width=0.5) + theme(text = element_text(size=14, color="black"), legend.position = "none")

#pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/cASE_DSG_DEG_enrichment.forest.noX.pdf", height=4, width=4)
#fp
#dev.off()



##################
### For gene expression, split up vs downregulated genes

# See if cASE is enriched in upregulated genes
# cASE(+) DE(+)
cASE_P <- sum(CM_ASE_GE$padj < 0.1)
cASE_N <- sum(CM_ASE_GE$padj >= 0.1)

cASE_P_DE_P <- sum(CM_ASE_GE$DE_logFC >= 0.25 & CM_ASE_GE$DE_padj < 0.1 & CM_ASE_GE$padj < 0.1)
cASE_P_DE_N <- cASE_P - cASE_P_DE_P

cASE_N_DE_P <- sum(CM_ASE_GE$DE_logFC >= 0.25 & CM_ASE_GE$DE_padj < 0.1 & CM_ASE_GE$padj >= 0.1)
cASE_N_DE_N <- cASE_N - cASE_N_DE_P

# Do fisher's test to see if cASE is enriched for DE
CM_up_enrich <-
matrix(c(cASE_P_DE_P, cASE_N_DE_P, cASE_P_DE_N, cASE_N_DE_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DE = c("DE+", "DE-")))
fish_test_CM_up <- fisher.test(CM_up_enrich, alternative = "two.sided")

## CM Down
# See if cASE is enriched in upregulated genes
# cASE(+) DE(+)
cASE_P <- sum(CM_ASE_GE$padj < 0.1)
cASE_N <- sum(CM_ASE_GE$padj >= 0.1)

cASE_P_DE_P <- sum(CM_ASE_GE$DE_logFC <= -0.25 & CM_ASE_GE$DE_padj < 0.1 & CM_ASE_GE$padj < 0.1)
cASE_P_DE_N <- cASE_P - cASE_P_DE_P

cASE_N_DE_P <- sum(CM_ASE_GE$DE_logFC <= -0.25 & CM_ASE_GE$DE_padj < 0.1 & CM_ASE_GE$padj >= 0.1)
cASE_N_DE_N <- cASE_N - cASE_N_DE_P

# Do fisher's test to see if cASE is enriched for DE
CM_down_enrich <-
matrix(c(cASE_P_DE_P, cASE_N_DE_P, cASE_P_DE_N, cASE_N_DE_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DE = c("DE+", "DE-")))
fish_test_CM_down <- fisher.test(CM_down_enrich, alternative = "two.sided")


# Do IPSCs
cASE_P <- sum(IPSC_ASE_GE$padj < 0.1)
cASE_N <- sum(IPSC_ASE_GE$padj >= 0.1)

cASE_P_DE_P <- sum(IPSC_ASE_GE$DE_logFC >= 0.25 & IPSC_ASE_GE$DE_padj < 0.1 & IPSC_ASE_GE$padj < 0.1)
cASE_P_DE_N <- cASE_P - cASE_P_DE_P

cASE_N_DE_P <- sum(IPSC_ASE_GE$DE_logFC >= 0.25 & IPSC_ASE_GE$DE_padj < 0.1 & IPSC_ASE_GE$padj >= 0.1)
cASE_N_DE_N <- cASE_N - cASE_N_DE_P

# Do fisher's test to see if cASE is enriched for DE
IPSC_up_enrich <-
matrix(c(cASE_P_DE_P, cASE_N_DE_P, cASE_P_DE_N, cASE_N_DE_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DE = c("DE+", "DE-")))
fish_test_IPSC_up <- fisher.test(IPSC_up_enrich, alternative = "two.sided")

## IPSC Down
# See if cASE is enriched in upregulated genes
# cASE(+) DE(+)
cASE_P <- sum(IPSC_ASE_GE$padj < 0.1)
cASE_N <- sum(IPSC_ASE_GE$padj >= 0.1)

cASE_P_DE_P <- sum(IPSC_ASE_GE$DE_logFC <= -0.25 & IPSC_ASE_GE$DE_padj < 0.1 & IPSC_ASE_GE$padj < 0.1)
cASE_P_DE_N <- cASE_P - cASE_P_DE_P

cASE_N_DE_P <- sum(IPSC_ASE_GE$DE_logFC <= -0.25 & IPSC_ASE_GE$DE_padj < 0.1 & IPSC_ASE_GE$padj >= 0.1)
cASE_N_DE_N <- cASE_N - cASE_N_DE_P

# Do fisher's test to see if cASE is enriched for DE
IPSC_down_enrich <-
matrix(c(cASE_P_DE_P, cASE_N_DE_P, cASE_P_DE_N, cASE_N_DE_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DE = c("DE+", "DE-")))
fish_test_IPSC_down <- fisher.test(IPSC_down_enrich, alternative = "two.sided")

# Do LCLs
cASE_P <- sum(LCL_ASE_GE$padj < 0.1)
cASE_N <- sum(LCL_ASE_GE$padj >= 0.1)

cASE_P_DE_P <- sum(LCL_ASE_GE$DE_logFC >= 0.25 & LCL_ASE_GE$DE_padj < 0.1 & LCL_ASE_GE$padj < 0.1)
cASE_P_DE_N <- cASE_P - cASE_P_DE_P

cASE_N_DE_P <- sum(LCL_ASE_GE$DE_logFC >= 0.25 & LCL_ASE_GE$DE_padj < 0.1 & LCL_ASE_GE$padj >= 0.1)
cASE_N_DE_N <- cASE_N - cASE_N_DE_P

# Do fisher's test to see if cASE is enriched for DE
LCL_up_enrich <-
matrix(c(cASE_P_DE_P, cASE_N_DE_P, cASE_P_DE_N, cASE_N_DE_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DE = c("DE+", "DE-")))
fish_test_LCL_up <- fisher.test(LCL_up_enrich, alternative = "two.sided")

## LCL Down
# See if cASE is enriched in upregulated genes
# cASE(+) DE(+)
cASE_P <- sum(LCL_ASE_GE$padj < 0.1)
cASE_N <- sum(LCL_ASE_GE$padj >= 0.1)

cASE_P_DE_P <- sum(LCL_ASE_GE$DE_logFC <= -0.25 & LCL_ASE_GE$DE_padj < 0.1 & LCL_ASE_GE$padj < 0.1)
cASE_P_DE_N <- cASE_P - cASE_P_DE_P

cASE_N_DE_P <- sum(LCL_ASE_GE$DE_logFC <= -0.25 & LCL_ASE_GE$DE_padj < 0.1 & LCL_ASE_GE$padj >= 0.1)
cASE_N_DE_N <- cASE_N - cASE_N_DE_P

# Do fisher's test to see if cASE is enriched for DE
LCL_down_enrich <-
matrix(c(cASE_P_DE_P, cASE_N_DE_P, cASE_P_DE_N, cASE_N_DE_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       DE = c("DE+", "DE-")))
fish_test_LCL_down <- fisher.test(LCL_down_enrich, alternative = "two.sided")

# Put data together and make a forest plot for up vs down, GE only
for_plot_upDown <- data.frame(CellType = c("LCL", "IPSC", "CM", "LCL", "IPSC", "CM"), AnalysisType = c("Up", "Up", "Up", "Down", "Down", "Down"), OR = c(fish_test_LCL_up$estimate, fish_test_IPSC_up$estimate, fish_test_CM_up$estimate, fish_test_LCL_down$estimate, fish_test_IPSC_down$estimate, fish_test_CM_down$estimate), Low = c(fish_test_LCL_up$conf.int[1], fish_test_IPSC_up$conf.int[1], fish_test_CM_up$conf.int[1], fish_test_LCL_down$conf.int[1], fish_test_IPSC_down$conf.int[1], fish_test_CM_down$conf.int[1]), Hi = c(fish_test_LCL_up$conf.int[2], fish_test_IPSC_up$conf.int[2], fish_test_CM_up$conf.int[2], fish_test_LCL_down$conf.int[2], fish_test_IPSC_down$conf.int[2], fish_test_CM_down$conf.int[2]), pval = c(fish_test_LCL_up$p.value, fish_test_IPSC_up$p.value, fish_test_CM_up$p.value, fish_test_LCL_down$p.value, fish_test_IPSC_down$p.value, fish_test_CM_down$p.value))

for_plot_upDown$sampleID <- paste0(for_plot_upDown$CellType, "_", for_plot_upDown$AnalysisType)

for_plot_upDown <- for_plot_upDown %>% arrange(AnalysisType, CellType)
for_plot_upDown$sampleID <- factor(for_plot_upDown$sampleID, levels=for_plot_upDown$sampleID)

# Plot
fp <- ggplot(data=for_plot_upDown, aes(x=sampleID, y=OR, ymin=Low, ymax=Hi)) +
geom_pointrange(aes(color=AnalysisType), size=0.8) +
geom_hline(yintercept=1, lty=2) +
coord_flip() +  # flip coordinates (puts labels on y axis)
xlab("Cell Type") + ylab("Enrichment (95% CI)") + theme_bw() + geom_errorbar(aes(ymin=Low, ymax=Hi, col=AnalysisType), width=0.5) + theme(text = element_text(size=14, color="black"), legend.position = "none")

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/cASE_DEG_upDown_enrichment.forest.noX.pdf", height=4, width=4)
fp
dev.off()


# Remake plot for splicing only
# Put data together and make a forest plot
for_plot_splice <- data.frame(CellType = c("LCL", "IPSC", "CM"), AnalysisType = c("Splicing", "Splicing", "Splicing"), OR = c(fish_test_LCL_DSG$estimate, fish_test_IPSC_DSG$estimate, fish_test_CM_DSG$estimate), Low = c(fish_test_LCL_DSG$conf.int[1], fish_test_IPSC_DSG$conf.int[1], fish_test_CM_DSG$conf.int[1]), Hi = c(fish_test_LCL_DSG$conf.int[2], fish_test_IPSC_DSG$conf.int[2], fish_test_CM_DSG$conf.int[2]), pval = c(fish_test_LCL_DSG$p.value, fish_test_IPSC_DSG$p.value, fish_test_CM_DSG$p.value))

for_plot_splice$sampleID <- paste0(for_plot_splice$CellType, "_", for_plot_splice$AnalysisType)

# Plot
fp <- ggplot(data=for_plot_splice, aes(x=sampleID, y=OR, ymin=Low, ymax=Hi)) +
geom_pointrange(aes(), size=0.8) +
geom_hline(yintercept=1, lty=2) +
coord_flip() +  # flip coordinates (puts labels on y axis)
xlab("Cell Type") + ylab("Enrichment (95% CI)") + theme_bw() + geom_errorbar(aes(ymin=Low, ymax=Hi), width=0.5) + theme(text = element_text(size=14, color="black"), legend.position = "none")

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/cASE_DSG_enrichment.forest.noX.pdf", height=2, width=4)
fp
dev.off()

#########################
#########################
# Now do OMIM intersection
# Add OMIM genes
OMIM <- read.delim("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/annotation/OMIM_morbidmap.121919.txt", sep="\t", header=F, comment.char="#")
colnames(OMIM) <- c("Phenotype", "Gene_Symbols", "MIM_Number", "Cyto_Location")

OMIM <- separate_rows(OMIM, Gene_Symbols, sep = ",")

OMIM_CM_Tr <- merge(mylm_CM_5df_Tr_anovaASE_gene, OMIM, by.x="g.id", by.y="Gene_Symbols")

OMIM_CM_Tr[grep("card", OMIM_CM_Tr$Phenotype),] %>% filter(padj < 0.1)
OMIM_CM_Tr[grep("oxidat", OMIM_CM_Tr$Phenotype),] %>% filter(padj < 0.1)
OMIM_CM_Tr[grep("rrhyth", OMIM_CM_Tr$Phenotype),] %>% filter(padj < 0.1)
OMIM_CM_Tr[grep("Gluco", OMIM_CM_Tr$Phenotype),] %>% filter(padj < 0.1)

forest_plot_CM("rs59914360_GM19209")
forest_plot_CM("rs35783144_GM18858")
forest_plot_CM("rs17103147_GM19209")


#  Combined oxidative phosphorylation deficiency 8,  Arrhythmogenic right ventricular dysplasia 11, {Glucocorticoid therapy, response to}, 

OMIM_LCL_Tr <- merge(mylm_LCL_5df_Tr_anovaASE_gene, OMIM, by.x="g.id", by.y="Gene_Symbols")
OMIM_LCL_Tr %>% filter(padj < 0.1) %>% select(Phenotype) %>% unique

OMIM_IPSC_Tr <- merge(mylm_IPSC_5df_Tr_anovaASE_gene, OMIM, by.x="g.id", by.y="Gene_Symbols")
OMIM_IPSC_Tr %>% filter(padj < 0.1) %>% select(Phenotype) %>% unique

######## Overlap with GxE database
gxe_db <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/FastQTL/GxE_mapping/normalized_metas/external_validation/datasets/GxE_full_db_build1/GxE_db_full.txt", header=T, sep="\t", stringsAsFactors=F)

# Just do Knowles CM paper
knowles2018 <- gxe_db[gxe_db$study == "Knowles2018",]

CM_sig_cASE <- mylm_CM_5df_Tr_anovaASE_gene[mylm_CM_5df_Tr_anovaASE_gene$padj < 0.1, ]
sum(CM_sig_cASE$ensg %in% knowles2018$ensg)
[1] 16 # was 14

unique(CM_sig_cASE[ which(CM_sig_cASE$ensg %in% knowles2018$ensg), "g.id"])
 [1] "FRAS1"    "PDGFC"    "MPHOSPH6" "HIBCH"    "LGI2"     "DGKB"
 [7] "PCM1"     "BBS2"     "GLIS3"    "COL22A1"

CM_nomSig_cASE <- mylm_CM_5df_Tr_anovaASE_gene[mylm_CM_5df_Tr_anovaASE_gene$p.value < 0.05, ]

sum(knowles2018$ensg %in% mylm_CM_5df_Tr_anovaASE_gene$ensg)
[1] 105 # was 105
sum(!knowles2018$ensg %in% mylm_CM_5df_Tr_anovaASE_gene$ensg)
[1] 342 # was 342
sum(knowles2018$ensg %in% CM_nomSig_cASE)
[1] 0
sum(knowles2018$ensg %in% CM_nomSig_cASE$ensg)
[1] 79 # was 79
sum(knowles2018$ensg %in% CM_sig_cASE$ensg)
[1] 12 # was 10 with X
79/105
[1] 0.75 # Was 0.704918
sum(unique(knowles2018$ensg) %in% mylm_CM_5df_Tr_anovaASE_gene$ensg)
[1] 105 # was 105
length(unique(knowles2018$ensg))
[1] 447

# Which 12 genes replicate at FDR significance?
unique(CM_sig_cASE[which(CM_sig_cASE$ensg %in% knowles2018$ensg), c("g.id", "term")])

IPSC_sig_cASE <- mylm_IPSC_5df_Tr_anovaASE_gene[mylm_IPSC_5df_Tr_anovaASE_gene$padj < 0.1, ]
sum(IPSC_sig_cASE$ensg %in% knowles2018$ensg)
[1] 22 # was 22

IPSC_nomSig_cASE <- mylm_IPSC_5df_Tr_anovaASE_gene[mylm_IPSC_5df_Tr_anovaASE_gene$p.value < 0.05, ]

sum(knowles2018$ensg %in% mylm_IPSC_5df_Tr_anovaASE_gene$ensg)
[1] 108 # was 109
sum(!knowles2018$ensg %in% mylm_IPSC_5df_Tr_anovaASE_gene$ensg)
[1] 339 # was 338
sum(knowles2018$ensg %in% IPSC_nomSig_cASE$ensg)
[1] 83 # was 83
sum(knowles2018$ensg %in% IPSC_sig_cASE$ensg)
[1] 16 # was 16
83/108
[1] 0.77 # was 0.76

### Alan previously did the sQTLs from Knowles with all cell types, but I think we just want the CMs.
sqtl <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/sQTL_dox/dox_sQTL_genes_unique.bed") %>%
  setNames(c("chr", "start", "end", "id", "type", "V1")) #740
resqtl <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/sQTL_dox/dox_re_sQTL_genes_unique.bed") %>%
  setNames(c("chr", "start", "end", "id", "type", "V1")) #62

# (re)sQTL genes with nominally significant cASE
mylm_CM_5df_Tr_anovaASE_gene %>% filter(g.id %in% resqtl$V1, p.value < 0.05) %>% pull(g.id) %>% unique() %>% length()
[1] 11
mylm_CM_5df_Tr_anovaASE_gene %>% filter(g.id %in% sqtl$V1, p.value < 0.05) %>% pull(g.id) %>% unique() %>% length()
[1] 88

# (re)sQTL genes with nominally significant ASE
myanova_reduced_CM_5df_finite$rsID <- sapply(strsplit(myanova_reduced_CM_5df_finite$SNP_Individual, "_"), "[", 1)
myanova_reduced_CM_5df_finite_gene <- merge(myanova_reduced_CM_5df_finite, gene_info, by="rsID")
myanova_reduced_CM_5df_finite_gene %>% filter(g.id %in% resqtl$V1, p.value < 0.05) %>% pull(g.id) %>% unique() %>% length()
[1] 16
myanova_reduced_CM_5df_finite_gene %>% filter(g.id %in% sqtl$V1, p.value < 0.05) %>% pull(g.id) %>% unique() %>% length()
[1] 182

  
save(list=c("ASE", "cv", "myanova_reduced_allCell", "mylm_CM_5df", "mylm_IPSC_5df", "mylm_LCL_5df", "vp_allCellTypes_gene"), file="Cell_Type_Specific_noReadCov_DF_ANOVA_add1_noX/anova_lm_vp.Rd")






## Compare genes variance in differentially spliced genes
# Load in differentially spliced genes (DSGs)
all_DSG <- bind_rows(LCL_DSG, IPSC_DSG, CM_DSG) %>% select(-Treatment.ID) %>% unique()

vp_allCellTypes_gene_melt_lof_DSG <- merge(vp_allCellTypes_gene_melt_lof, all_DSG, by=c("ensg", "CellType"), all.x=T)

# Replace NA with 0
vp_allCellTypes_gene_melt_lof_DSG[which(is.na(vp_allCellTypes_gene_melt_lof_DSG$DSG)), "DSG"] <- 0

# DSG cells together
summary(lm(value ~ DSG + rankTotal, data=vp_allCellTypes_gene_melt_lof_DSG[ vp_allCellTypes_gene_melt_lof_DSG$variable == "Residual",]))
Residuals:
     Min       1Q   Median       3Q      Max
-0.18907 -0.03176 -0.00257  0.01823  0.47934

Coefficients:
              Estimate Std. Error  t value Pr(>|t|)
(Intercept)  2.002e-01  1.096e-03  182.596   <2e-16 ***
DSG         -2.062e-04  1.042e-03   -0.198    0.843
rankTotal   -3.081e-05  2.568e-07 -119.973   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06329 on 16497 degrees of freedom
Multiple R-squared:  0.4741,    Adjusted R-squared:  0.4741
F-statistic:  7437 on 2 and 16497 DF,  p-value: < 2.2e-16

summary(lm(value ~ DSG + rankTotal, data=vp_allCellTypes_gene_melt_lof_DSG[ vp_allCellTypes_gene_melt_lof_DSG$variable == "Tr",]))

Residuals:
     Min       1Q   Median       3Q      Max
-0.10609 -0.03005 -0.00267  0.01603  1.02996

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.044e-01  1.049e-03  99.525   <2e-16 ***
DSG          2.131e-03  9.973e-04   2.137   0.0326 *
rankTotal   -1.805e-05  2.457e-07 -73.460   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06057 on 16497 degrees of freedom
Multiple R-squared:  0.2507,    Adjusted R-squared:  0.2506
F-statistic:  2760 on 2 and 16497 DF,  p-value: < 2.2e-16

# Plot the treatment variance by whether the gene is DSG (after regressing out expression)
vp_allCellTypes_gene_melt_lof_DSG %>% filter(variable == "Tr") %>% mutate(regress = resid(lm(value ~ rankTotal, data = .))) -> vp_allCellTypes_gene_melt_lof_Tr_regress

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes.DSG.Regress.Tr.box.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_lof_Tr_regress %>% filter(variable == "Tr"), aes(x=factor(DSG), y=regress, fill=factor(DSG)), color="black") + geom_boxplot(notch=TRUE, outlier.shape=NA) + theme_bw() + labs(title="Treatment Variance by whether gene is DSG") + xlab("Group") + ylab("ASE variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + ylim(-0.12,0.12)
dev.off()


### 
# Try Mendelian cardiomyopathy genes (selected from Mayo Clinic Lab's Comprehensive Cardiomyopathy Multi-Gene Panel; https://www.mayocliniclabs.com/test-catalog/Overview/63164)
Mendelian_cardio <- c("ABCC9", "ACTC1", "ACTN2", "ANKRD1", "BRAF", "CAV3", "CBL", "CRYAB", "CSRP3", "DES", "DSC2", "DSG2", "DSP", "DTNA", "GLA", "HRAS", "JUP", "KRAS", "LAMA4", "LAMP2", "LDB3", "LMNA", "MAP2K1", "MAP2K2", "MYBPC3", "MYH6", "MYH7", "MYL2", "MYL3", "MYLK2", "MYOZ2", "MYPN", "NEXN", "NRAS", "PKP2", "PLN", "PRKAG2", "PTPN11", "RAF1", "RBM20", "RYR2", "SCN5A", "SGCD", "SHOC2", "SOS1", "TAZ", "TCAP", "TMEM43", "TNNC1", "TNNI3", "TNNT2", "TPM1", "TTN", "TTR", "VCL")

# I just want to do this for CM's I think. Create column for Mendelian cardiomyopathy disease yes/no

vp_allCellTypes_gene_melt_lof_resid %>% filter(CellType == "CM") %>% mutate(Mend_Card = "No") -> vp_CM_gene_melt_MendCard

vp_CM_gene_melt_MendCard[which(vp_CM_gene_melt_MendCard$g.id %in% Mendelian_cardio), "Mend_Card"] <- "Yes"

summary(lm(value ~ Mend_Card + rankTotal, data=vp_CM_gene_melt_MendCard %>% filter(variable == "Residual")))
Residuals:
     Min       1Q   Median       3Q      Max
-0.16716 -0.02835 -0.00260  0.01744  0.42684

Coefficients:
               Estimate Std. Error t value Pr(>|t|)
(Intercept)   1.831e-01  1.754e-03 104.427  < 2e-16 ***
Mend_CardYes  2.666e-02  9.028e-03   2.953  0.00316 **
rankTotal    -2.922e-05  4.168e-07 -70.092  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05826 on 5407 degrees of freedom
Multiple R-squared:  0.4761,    Adjusted R-squared:  0.4759
F-statistic:  2457 on 2 and 5407 DF,  p-value: < 2.2e-16

summary(lm(value ~ Mend_Card + rankTotal, data=vp_CM_gene_melt_MendCard %>% filter(variable == "Tr")))
Call:
lm(formula = value ~ Mend_Card + rankTotal, data = vp_CM_gene_melt_MendCard %>%
    filter(variable == "Tr"))

Residuals:
     Min       1Q   Median       3Q      Max
-0.10337 -0.02771 -0.00203  0.01681  0.59777

Coefficients:
               Estimate Std. Error t value Pr(>|t|)
(Intercept)   1.037e-01  1.738e-03  59.657   <2e-16 ***
Mend_CardYes -1.580e-03  8.946e-03  -0.177     0.86
rankTotal    -1.847e-05  4.130e-07 -44.711   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05772 on 5407 degrees of freedom
Multiple R-squared:  0.2701,    Adjusted R-squared:  0.2699
F-statistic:  1001 on 2 and 5407 DF,  p-value: < 2.2e-16

# Plot the residual variance by whether the gene is in a Mendelian cardiomyopathy gene (after regressing out expression)
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/CM.MendCard.Regress.Resid.box.noX.pdf")
ggplot(vp_CM_gene_melt_MendCard %>% filter(variable == "Residual"), aes(x=Mend_Card, y=regress, fill=Mend_Card), color="black") + geom_boxplot(notch=TRUE, outlier.shape=NA) + theme_bw() + labs(title="ASE Variance by whether gene is in Mendelian Cardiomyopathy") + xlab("Group") + ylab("ASE variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + ylim(-0.12,0.1)
dev.off()

### Try TWAS vs non-TWAS genes
TWAS <- read.table("/nfs/rprdata/ALOFT/AL1-6_ln/FastQTL/GxE_mapping/normalized_metas/TWAS_overlap/media-4(2).txt", header=T, stringsAsFactors=)
TWAS$ensg <- gsub("\\..*", "", TWAS$Gene)

# All of the TWAS genes in the file are significant
vp_allCellTypes_gene_melt_lof %>% mutate(TWAS = "No") -> vp_allCellTypes_gene_melt_TWAS

vp_allCellTypes_gene_melt_TWAS[which(vp_allCellTypes_gene_melt_TWAS$ensg %in% unique(TWAS$ensg)), "TWAS"] <- "Yes"

# Most genes are in TWAS
table(vp_allCellTypes_gene_melt_TWAS$TWAS)

   No   Yes
12048  53952     # was 14732 55580

# Test for relationship
summary(lm(value ~ TWAS + rankTotal, data=vp_allCellTypes_gene_melt_TWAS %>% filter(variable == "Residual")))
Residuals:
     Min       1Q   Median       3Q      Max
-0.18947 -0.03161 -0.00272  0.01838  0.48070

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  2.054e-01  1.482e-03  138.63  < 2e-16 ***
TWASYes     -6.641e-03  1.275e-03   -5.21 1.92e-07 ***
rankTotal   -3.078e-05  2.526e-07 -121.89  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06324 on 16497 degrees of freedom
Multiple R-squared:  0.475,     Adjusted R-squared:  0.4749
F-statistic:  7462 on 2 and 16497 DF,  p-value: < 2.2e-16

# Now try Tr variance
summary(lm(value ~ TWAS + rankTotal, data=vp_allCellTypes_gene_melt_TWAS %>% filter(variable == "Tr")))

Residuals:
     Min       1Q   Median       3Q      Max
-0.10465 -0.03016 -0.00265  0.01599  1.02940

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.044e-01  1.420e-03  73.540   <2e-16 ***
TWASYes      5.465e-04  1.221e-03   0.448    0.654   # No longer significant
rankTotal   -1.796e-05  2.419e-07 -74.244   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.06057 on 16497 degrees of freedom
Multiple R-squared:  0.2505,    Adjusted R-squared:  0.2504
F-statistic:  2757 on 2 and 16497 DF,  p-value: < 2.2e-16


# Plot the residual variance by whether the gene is in a Mendelian cardiomyopathy gene (after regressing out expression)
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/CM.MendCard.Regress.Resid.box.noX.pdf")
ggplot(vp_CM_gene_melt_MendCard %>% filter(variable == "Residual"), aes(x=Mend_Card, y=regress, fill=Mend_Card), color="black") + geom_boxplot(notch=TRUE, outlier.shape=NA) + theme_bw() + labs(title="ASE Variance by whether gene is in Mendelian Cardiomyopathy") + xlab("Group") + ylab("ASE variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + ylim(-0.12,0.1)
dev.off()

# Add dN/dS
dNdS <- read.table("/nfs/rprdata/Anthony/data/29mammals/Overall_dN_dS.bed", header=F, stringsAsFactors=F)
colnames(dNdS)[c(4,11)] <- c("g.id", "dNdS")

vp_allCellTypes_gene_melt_dNdS <- merge(vp_allCellTypes_gene_melt, dNdS[,c("g.id", "dNdS")], by = "g.id") # 2848/6310 genes had dN/dS score


# Plot dN/dS distribution for all reported genes in 29 mammals file and just those which have ANOVA ASE.
dNdS$Group <- "All"
vp_allCellTypes_gene_melt_dNdS %>% select(g.id, dNdS) %>% unique() %>% mutate(Group = "ANOVA ASE") -> dNdS_ANOVA

dNdS_forPlot <- rbind(dNdS[,c("g.id", "dNdS", "Group")], dNdS_ANOVA)

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/dNdS_dist.hist.noX.pdf")
ggplot(dNdS_forPlot, aes(x=dNdS, fill=Group)) + geom_histogram(position="dodge") + theme_bw() + labs(title="dN/dS distribution") + xlab("dN/dS") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()



summary(lm(value ~ dNdS + rankTotal, data=vp_allCellTypes_gene_melt_dNdS %>% filter(variable == "Residual")))
Residuals:
     Min       1Q   Median       3Q      Max
-0.17749 -0.02904 -0.00207  0.01702  0.42365

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.867e-01  1.645e-03 113.506  < 2e-16 ***
dNdS         1.614e-02  4.530e-03   3.562  0.00037 ***
rankTotal   -2.901e-05  3.245e-07 -89.384  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05927 on 8947 degrees of freedom
Multiple R-squared:  0.4733,    Adjusted R-squared:  0.4731
F-statistic:  4019 on 2 and 8947 DF,  p-value: < 2.2e-16

# Now try Tr variance
summary(lm(value ~ dNdS + rankTotal, data=vp_allCellTypes_gene_melt_dNdS %>% filter(variable == "Tr")))

Residuals:
     Min       1Q   Median       3Q      Max
-0.10661 -0.02859 -0.00222  0.01578  0.83722

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  1.076e-01  1.650e-03  65.185   <2e-16 ***
dNdS        -3.154e-03  4.545e-03  -0.694    0.488
rankTotal   -1.821e-05  3.256e-07 -55.919   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.05946 on 8947 degrees of freedom
Multiple R-squared:  0.2591,    Adjusted R-squared:  0.2589
F-statistic:  1565 on 2 and 8947 DF,  p-value: < 2.2e-16

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/all.dNdS.Resid.2dens.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_dNdS %>% filter(variable == "Residual"), aes(x=dNdS, y=value)) +  geom_bin2d(bins=50) + theme_bw() + labs(title="ASE Variance by dN/dS score") + xlab("dN/dS") + ylab("ASE variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) #+ ylim(-0.12,0.1)
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/all.dNdS.Tr.2dens.noX.pdf")
ggplot(vp_allCellTypes_gene_melt_dNdS %>% filter(variable == "Tr"), aes(x=dNdS, y=value)) +  geom_bin2d(bins=50) + theme_bw() + labs(title="ASE Treatment Variance by dN/dS score") + xlab("dN/dS") + ylab("ASE variance") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) #+ ylim(-0.12,0.1)
dev.off()


######################
######################
# Plot the ASE from linear model for genes belonging to groups which we used to analyze ASE variance (DEG, DSG, GC, TATA, etc)

# I want to average the treatment ASE measures and add them to the intercept value for each SNP_Individual
mylm_LCL_5df_Tr %>% group_by(SNP_Individual) %>% summarise(Tr_avg = mean(estimate)) -> LCL_avg_TrASE
mylm_LCL_5df %>% filter(term == "(Intercept)") %>% inner_join(LCL_avg_TrASE, by="SNP_Individual") %>% mutate(ASE_val = estimate + Tr_avg, CellType = "LCL") -> mylm_LCL_5df_ASE_val

mylm_IPSC_5df_Tr %>% group_by(SNP_Individual) %>% summarise(Tr_avg = mean(estimate)) -> IPSC_avg_TrASE
mylm_IPSC_5df %>% filter(term == "(Intercept)") %>% inner_join(IPSC_avg_TrASE, by="SNP_Individual") %>% mutate(ASE_val = estimate + Tr_avg, CellType = "IPSC") -> mylm_IPSC_5df_ASE_val

mylm_CM_5df_Tr %>% group_by(SNP_Individual) %>% summarise(Tr_avg = mean(estimate)) -> CM_avg_TrASE
mylm_CM_5df %>% filter(term == "(Intercept)") %>% inner_join(CM_avg_TrASE, by="SNP_Individual") %>% mutate(ASE_val = estimate + Tr_avg, CellType = "CM") -> mylm_CM_5df_ASE_val

# Add this info to table with ASE variance because it already has annotations for some categories.
rbind(mylm_LCL_5df_ASE_val, mylm_IPSC_5df_ASE_val, mylm_CM_5df_ASE_val) %>% select(SNP_Individual, ASE_val, CellType) %>% right_join(vp_allCellTypes_gene, by=c("SNP_Individual", "CellType")) -> vp_allCellTypes_gene_ASEval

vp_allCellTypes_gene_ASEval_counts <- merge(vp_allCellTypes_gene_ASEval, all_counts, by=c("SNP_Individual", "CellType"))


# Plot ecdf for different categories
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.DE.ecdf.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval, aes(x=ASE_val, color = DE)) + stat_ecdf() + theme_bw() + labs(title="DEG vs non DEG ASE: ecdf") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

# KS test for not_DE vs DE
ks.test(vp_allCellTypes_gene_ASEval %>% filter(DE == "Not_DE") %>% pull(ASE_val), vp_allCellTypes_gene_ASEval %>% filter(DE == "DE") %>% pull(ASE_val)) # D = 0.041, pval = 1.4 x 10^(-5) ; same as with X

# KS test for DE vs More DE
ks.test(vp_allCellTypes_gene_ASEval %>% filter(DE == "DE") %>% pull(ASE_val), vp_allCellTypes_gene_ASEval %>% filter(DE == "More DE") %>% pull(ASE_val)) # D = 0.047, pval = 1.28 x 10^(-5) # was 0.043

# Barplot of absolute value of ASE variance for DEG
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.DEG.box.noX.pdf", width=4.5)
ggplot(vp_allCellTypes_gene_ASEval, aes(x = DE, y=abs(ASE_val), fill = DE)) + geom_boxplot(notch=T) + theme_bw() + labs(title="Absolute ASE by DEG status") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

summary(lm(abs(ASE_val) ~ as.integer(DE) + rankTotal, data=vp_allCellTypes_gene_ASEval_counts))
Residuals:
     Min       1Q   Median       3Q      Max
-0.99490 -0.26572 -0.07284  0.20009  2.04620

Coefficients:
                 Estimate Std. Error t value Pr(>|t|)
(Intercept)     9.116e-01  1.135e-02  80.302  < 2e-16 ***
as.integer(DE)  3.565e-02  4.713e-03   7.563 4.09e-14 ***
rankTotal      -9.551e-05  1.480e-06 -64.549  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.425 on 21419 degrees of freedom
Multiple R-squared:  0.1688,    Adjusted R-squared:  0.1687
F-statistic:  2175 on 2 and 21419 DF,  p-value: < 2.2e-16



pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.TATA.ecdf.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval, aes(x=ASE_val, color = TATA)) + stat_ecdf() + theme_bw() + labs(title="TATA ASE: ecdf") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

# KS test
ks.test(vp_allCellTypes_gene_ASEval %>% filter(TATA == "TATA") %>% pull(ASE_val), vp_allCellTypes_gene_ASEval %>% filter(TATA == "No_TATA") %>% pull(ASE_val)) # D = 0.06, pval = 1.9 x 10^(-14)

summary(lm(abs(ASE_val) ~ TATA + rankTotal, data=vp_allCellTypes_gene_ASEval_counts))


# Barplot of absolute value of ASE variance for TATA
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.TATA.box.noX.pdf", width=3.8)
ggplot(vp_allCellTypes_gene_ASEval, aes(x = TATA, y=abs(ASE_val), fill = TATA)) + geom_boxplot(notch=T) + theme_bw() + labs(title="Absolute ASE by TATA status") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()


# For CpG and GC use a scatter plot
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.GC.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval, aes(x=as.numeric(GC.), y = ASE_val, fill = CellType)) + geom_point(pch=21) + theme_bw() + labs(title="ASE value based on GC") + xlab("Promoter GC content %") + ylab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.CpG.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval, aes(x=as.numeric(CpG.), y = ASE_val, fill = CellType)) + geom_point(pch=21) + theme_bw() + labs(title="ASE value based on CpG") + xlab("Promoter CpG proportion") + ylab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

# Also plot absolute value for scatter plot
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells_abs.GC.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval, aes(x=as.numeric(GC.), y = abs(ASE_val), fill = CellType)) + geom_point(pch=21) + theme_bw() + labs(title="Absolute ASE value based on GC") + xlab("Promoter GC content %") + ylab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells_abs.CpG.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval, aes(x=as.numeric(CpG.), y = abs(ASE_val), fill = CellType)) + geom_point(pch=21) + theme_bw() + labs(title="Absolute ASE value based on CpG") + xlab("Promoter CpG proportion") + ylab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

# See if it looks better as 2d density plot
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells_abs.GC.2dens.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval, aes(x=as.numeric(GC.), y = abs(ASE_val))) + geom_bin2d(bins=50) + theme_bw() + labs(title="Absolute ASE value based on GC") + xlab("Promoter GC content %") + ylab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + scale_fill_continuous(type = "viridis")
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells_abs.CpG.2dens.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval, aes(x=as.numeric(CpG.), y = abs(ASE_val))) + geom_bin2d(bins=50) + theme_bw() + labs(title="Absolute ASE value based on CpG") + xlab("Promoter CpG proportion") + ylab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + scale_fill_continuous(type = "viridis")
dev.off()

# Calculate correlations
rcorr(as.numeric(vp_allCellTypes_gene_ASEval$CpG.), abs(vp_allCellTypes_gene_ASEval$ASE_val), type="spearman") # cor = -0.09, p = 0
rcorr(as.numeric(vp_allCellTypes_gene_ASEval$GC.), abs(vp_allCellTypes_gene_ASEval$ASE_val), type="spearman") # cor = -0.07, p = 0

summary(lm(abs(ASE_val) ~ as.numeric(CpG.) + rankTotal, data=vp_allCellTypes_gene_ASEval_counts))

Residuals:
     Min       1Q   Median       3Q      Max
-0.95065 -0.24702 -0.06424  0.18223  1.87194

Coefficients:
                   Estimate Std. Error t value Pr(>|t|)
(Intercept)       9.747e-01  7.969e-03 122.312  < 2e-16 ***
as.numeric(CpG.) -4.240e-01  5.717e-02  -7.417 1.26e-13 ***
rankTotal        -9.151e-05  1.569e-06 -58.316  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4002 on 17166 degrees of freedom
  (4253 observations deleted due to missingness)
Multiple R-squared:  0.1719,    Adjusted R-squared:  0.1718
F-statistic:  1782 on 2 and 17166 DF,  p-value: < 2.2e-16

summary(lm(abs(ASE_val) ~ as.numeric(GC.) + rankTotal, data=vp_allCellTypes_gene_ASEval_counts))
Residuals:
     Min       1Q   Median       3Q      Max
-0.95856 -0.24799 -0.06417  0.18218  1.87418

Coefficients:
                  Estimate Std. Error t value Pr(>|t|)
(Intercept)      1.017e+00  1.690e-02  60.195  < 2e-16 ***
as.numeric(GC.) -1.190e-01  2.454e-02  -4.851 1.24e-06 ***
rankTotal       -9.205e-05  1.567e-06 -58.730  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4006 on 17162 degrees of freedom
  (4257 observations deleted due to missingness)
Multiple R-squared:  0.1705,    Adjusted R-squared:  0.1704
F-statistic:  1764 on 2 and 17162 DF,  p-value: < 2.2e-16


###
# Now do DSG and lof
vp_allCellTypes_gene_melt_lof_DSG %>% select(SNP_Individual, CellType, oe_lof, DSG) %>% distinct() -> vp_allCellTypes_gene_melt_lof_DSG_small

rbind(mylm_LCL_5df_ASE_val, mylm_IPSC_5df_ASE_val, mylm_CM_5df_ASE_val) %>% select(SNP_Individual, ASE_val, CellType) %>% right_join(vp_allCellTypes_gene_melt_lof_DSG_small, by=c("SNP_Individual", "CellType")) -> vp_allCellTypes_gene_ASEval_DSG_lof

vp_allCellTypes_gene_ASEval_DSG_lof_counts <- merge(vp_allCellTypes_gene_ASEval_DSG_lof, all_counts, by=c("SNP_Individual", "CellType"))


# DSG plot
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.DSG.ecdf.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval_DSG_lof, aes(x=ASE_val, color = factor(DSG))) + stat_ecdf() + theme_bw() + labs(title="DSG ASE: ecdf") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

# KS test for DSG vs Not DSG
ks.test(vp_allCellTypes_gene_ASEval_DSG_lof %>% filter(DSG == 0) %>% pull(ASE_val), vp_allCellTypes_gene_ASEval_DSG_lof %>% filter(DSG == 1) %>% pull(ASE_val)) # D = 0.041, pval = 4.8 x 10^(-6) # Was 0.054

summary(lm(abs(ASE_val) ~ DSG + rankTotal, data=vp_allCellTypes_gene_ASEval_DSG_lof_counts))
Residuals:
     Min       1Q   Median       3Q      Max
-0.91541 -0.24642 -0.06446  0.18064  1.88365

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  9.376e-01  6.933e-03 135.241   <2e-16 ***
DSG          5.711e-03  6.593e-03   0.866    0.386
rankTotal   -9.309e-05  1.624e-06 -57.316   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4001 on 16477 degrees of freedom
Multiple R-squared:   0.17,     Adjusted R-squared:  0.1699
F-statistic:  1688 on 2 and 16477 DF,  p-value: < 2.2e-16


# Absoulte value boxplot for DSG
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.DSG.box.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval_DSG_lof, aes(x = factor(DSG), y=abs(ASE_val), fill = factor(DSG))) + geom_boxplot(notch=T) + theme_bw() + labs(title="Absolute ASE by DSG status") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()


# For lof use a scatter plot
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.lof.scatter.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval_DSG_lof, aes(x=as.numeric(oe_lof), y = ASE_val, fill = CellType)) + geom_point(pch=21) + theme_bw() + labs(title="ASE value based on lof") + xlab("Observed / expected lof") + ylab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells_abs.lof.2dens.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval_DSG_lof, aes(x=as.numeric(oe_lof), y = abs(ASE_val))) + geom_bin2d(bins=50) + theme_bw() + labs(title="Absolute ASE value based on LoF") + xlab("LoF Observed / Expected") + ylab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + scale_fill_continuous(type = "viridis")
dev.off()

# Calculate correlations
rcorr(as.numeric(vp_allCellTypes_gene_ASEval_DSG_lof$oe_lof), abs(vp_allCellTypes_gene_ASEval_DSG_lof$ASE_val), type="spearman") # cor = 0.1, p = 0 ; was 0.07

summary(lm(abs(ASE_val) ~ oe_lof + rankTotal, data=vp_allCellTypes_gene_ASEval_DSG_lof_counts))
Residuals:
     Min       1Q   Median       3Q      Max
-0.93374 -0.24667 -0.06489  0.18171  1.91265

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  9.042e-01  8.545e-03 105.821  < 2e-16 ***
oe_lof       6.179e-02  9.207e-03   6.712 1.98e-11 ***
rankTotal   -9.149e-05  1.612e-06 -56.755  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.3996 on 16410 degrees of freedom
  (67 observations deleted due to missingness)
Multiple R-squared:  0.1719,    Adjusted R-squared:  0.1718
F-statistic:  1703 on 2 and 16410 DF,  p-value: < 2.2e-16

###
# Mendelian cardiomyopathy genes 
vp_CM_gene_melt_MendCard %>% select(SNP_Individual, CellType, Mend_Card) %>% distinct() -> vp_allCellTypes_gene_melt_lof_MendCard_small

mylm_CM_5df_ASE_val %>% select(SNP_Individual, ASE_val, CellType) %>% right_join(vp_allCellTypes_gene_melt_lof_MendCard_small, by=c("SNP_Individual", "CellType")) -> vp_allCellTypes_gene_ASEval_MendCard

vp_allCellTypes_gene_ASEval_MendCard_counts <- merge(vp_allCellTypes_gene_ASEval_MendCard, all_counts, by=c("SNP_Individual", "CellType"))


pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.MendCard.ecdf.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval_MendCard, aes(x=ASE_val, color = Mend_Card)) + stat_ecdf() + theme_bw() + labs(title="Mendelian Cardiomyopathy ASE: ecdf") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

# KS test for Mend_Card vs Not Mend_Card
ks.test(vp_allCellTypes_gene_ASEval_MendCard %>% filter(Mend_Card == "No") %>% pull(ASE_val), vp_allCellTypes_gene_ASEval_MendCard %>% filter(Mend_Card == "Yes") %>% pull(ASE_val)) # D = 0.11, pval = 0.66

summary(lm(abs(ASE_val) ~ Mend_Card + rankTotal, data=vp_allCellTypes_gene_ASEval_MendCard_counts))
Residuals:
     Min       1Q   Median       3Q      Max
-0.88761 -0.23422 -0.06497  0.17326  1.83960

Coefficients:
               Estimate Std. Error t value Pr(>|t|)
(Intercept)   9.094e-01  1.207e-02  75.322   <2e-16 ***
Mend_CardYes  1.003e-01  6.095e-02   1.645      0.1
rankTotal    -9.071e-05  2.873e-06 -31.577   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.3932 on 5198 degrees of freedom
Multiple R-squared:  0.161,     Adjusted R-squared:  0.1607
F-statistic: 498.8 on 2 and 5198 DF,  p-value: < 2.2e-16


# Absoulte value boxplot for Mend_Card
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.MendCard.box.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval_MendCard, aes(x = Mend_Card, y=abs(ASE_val), fill = Mend_Card)) + geom_boxplot(notch=T) + theme_bw() + labs(title="Absolute ASE by Mendelian Cardiomyopathy status") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()


###
# TWAS genes
vp_allCellTypes_gene_melt_TWAS %>% select(SNP_Individual, CellType, TWAS) %>% distinct() -> vp_allCellTypes_gene_melt_lof_TWAS_small

rbind(mylm_LCL_5df_ASE_val, mylm_IPSC_5df_ASE_val, mylm_CM_5df_ASE_val) %>% select(SNP_Individual, ASE_val, CellType) %>% right_join(vp_allCellTypes_gene_melt_lof_TWAS_small, by=c("SNP_Individual", "CellType")) -> vp_allCellTypes_gene_ASEval_TWAS

vp_allCellTypes_gene_ASEval_TWAS_counts <- merge(vp_allCellTypes_gene_ASEval_TWAS, all_counts, by=c("SNP_Individual", "CellType"))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.TWAS.ecdf.noX.pdf")
ggplot(vp_allCellTypes_gene_ASEval_TWAS, aes(x=ASE_val, color = TWAS)) + stat_ecdf() + theme_bw() + labs(title=" TWAS ASE: ecdf") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()

# KS test for TWAS vs Not TWAS
ks.test(vp_allCellTypes_gene_ASEval_TWAS %>% filter(TWAS == "No") %>% pull(ASE_val), vp_allCellTypes_gene_ASEval_TWAS %>% filter(TWAS == "Yes") %>% pull(ASE_val)) # D = 0.050, pval = 1.4 x 10^(-5) ; was 0.058

summary(lm(abs(ASE_val) ~ TWAS + rankTotal, data=vp_allCellTypes_gene_ASEval_TWAS_counts))
Residuals:
     Min       1Q   Median       3Q      Max
-0.93157 -0.24685 -0.06494  0.17997  1.88992

Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)  9.609e-01  9.548e-03 100.639  < 2e-16 ***
TWASYes     -2.829e-02  8.204e-03  -3.449 0.000565 ***
rankTotal   -9.238e-05  1.627e-06 -56.784  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4003 on 15937 degrees of freedom
Multiple R-squared:  0.1693,    Adjusted R-squared:  0.1691
F-statistic:  1623 on 2 and 15937 DF,  p-value: < 2.2e-16


# Absoulte value boxplot for TWAS
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells.TWAS.box.noX.pdf", width=3.5)
ggplot(vp_allCellTypes_gene_ASEval_TWAS, aes(x = TWAS, y=abs(ASE_val), fill = TWAS)) + geom_boxplot(notch=T) + theme_bw() + labs(title="Absolute ASE by TWAS status") + xlab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15))
dev.off()


### How many genes with CM treatment cASE are in a CVD TWAS?
mylm_CM_5df_Tr_anovaASE_gene %>% filter(padj < 0.1) %>% pull(ensg) %>% unique() -> CM_Tr_cASE_ensg

CVD_TWAS_ensg <- TWAS %>% filter(Trait %in% c("CARDIoGRAM_C4D_CAD", "GLGC_Mc_HDL", "GLGC_Mc_LDL", "GLGC_Mc_TG", "UKB_20002_1065_self_reported_hypertension", "UKB_20002_1473_self_reported_high_cholesterol", "HRGene_HeartRate", "ICBP_DiastolicPressure", "ICBP_SystolicPressure", "MAGNETIC_HDL.C", "MAGNETIC_IDL.TG", "MAGNETIC_LDL.C", "UKB_20002_1094_self_reported_deep_venous_thrombosis_dvt", "UKB_6150_1_Vascular_or_heart_problems_diagnosed_by_doctor_Heart_attack")) %>% pull(ensg) %>% unique()


sum(CM_Tr_cASE_ensg %in% CVD_TWAS_ensg)
[1] 52

## Do dN/dS mean ASE
rbind(mylm_LCL_5df_ASE_val, mylm_IPSC_5df_ASE_val, mylm_CM_5df_ASE_val) %>% select(SNP_Individual, ASE_val, CellType) %>% right_join(vp_allCellTypes_gene_melt_dNdS, by=c("SNP_Individual", "CellType")) -> vp_allCellTypes_gene_dNdS_ASEval 

vp_allCellTypes_gene_dNdS_ASEval_counts <- merge(vp_allCellTypes_gene_dNdS_ASEval, all_counts, by=c("SNP_Individual", "CellType"))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/ASEval_allCells_abs.dNdS.2dens.noX.pdf")
ggplot(vp_allCellTypes_gene_dNdS_ASEval, aes(x=as.numeric(dNdS), y = abs(ASE_val))) + geom_bin2d(bins=50) + theme_bw() + labs(title="Absolute ASE value based on dN/dS") + xlab("dN/dS") + ylab("ASE value") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + scale_fill_continuous(type = "viridis")
dev.off()

rcorr(as.numeric(vp_allCellTypes_gene_dNdS_ASEval$dNdS), abs(vp_allCellTypes_gene_dNdS_ASEval$ASE_val), type="spearman") # cor = 0.08, p = 0




##############################
##############################
# Write tables and files

# Make table of FDR significant cASE which is present in Knowles 2018
knowles_sig_write <- unique(CM_sig_cASE[which(CM_sig_cASE$ensg %in% knowles2018$ensg), ])
knowles_sig_write <- knowles_sig_write %>% select(SNP_Individual, term, Treatment.ID, estimate, std.error, statistic, p.value, padj, ensg, g.id)

write.table(knowles_sig_write, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/Knowles_sig_cASE_noX.tab", quote=F, row.names=F, col.names=T, sep="\t")




########TWAS 

# Plot TWAS Residual and Tr variance after correction
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.TWAS.Resid.regress.box.noX.pdf", width = 3.5)
ggplot(vp_allCellTypes_gene_melt_TWAS %>% filter(variable == "Residual"), aes(x=TWAS, y=regress, fill=TWAS), color="black") + geom_boxplot(notch=TRUE, outlier.shape=NA) + theme_bw() + labs(title="ASE Variance by TWAS status") + xlab("Group") + ylab("ASE variance after regressing rankTotal") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + ylim(-0.1, 0.1)
dev.off()

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/varPart/AllCellTypes_combined.TWAS.Tr.regress.box.noX.pdf", width=3.5)
ggplot(vp_allCellTypes_gene_melt_TWAS %>% filter(variable == "Tr"), aes(x=TWAS, y=regress, fill=TWAS), color="black") + geom_boxplot(notch=TRUE, outlier.shape=NA) + theme_bw() + labs(title="ASE Variance by TWAS status") + xlab("Group") + ylab("ASE variance after regressing rankTotal") + theme(text = element_text(size=12), axis.text.x = element_text(size=13), axis.text.y = element_text(size=15)) + ylim(-0.1, 0.1)
dev.off()









# Roger wants me to do Mann-Whitney (Wilcox) test to see if there's any difference between DE/non-DE genes and TWAS/non-TWAS genes and Tr variance
wilcox.test(regress ~ DE, vp_allCellTypes_gene_melt_lof_resid %>% filter(variable == "Tr", DE %in% c("Not_DE", "DE")))

        Wilcoxon rank sum test with continuity correction

data:  regress by DE
W = 19332102, p-value = 1.175e-05
alternative hypothesis: true location shift is not equal to 0

wilcox.test(regress ~ DE, vp_allCellTypes_gene_melt_lof_resid %>% filter(variable == "Tr", DE %in% c("DE", "More DE")))

        Wilcoxon rank sum test with continuity correction

data:  regress by DE
W = 13703378, p-value = 0.1152
alternative hypothesis: true location shift is not equal to 0

# Do the same for TWAS:
vp_allCellTypes_gene_melt_TWAS$regress <- resid(lm(value ~ rankTotal, data = vp_allCellTypes_gene_melt_TWAS))

wilcox.test(regress ~ TWAS, vp_allCellTypes_gene_melt_TWAS %>% filter(variable == "Tr"))



###########################
###########################
# Write tables
myanova_reduced_allCell_forWrite <- myanova_reduced_allCell %>% select(SNP_Individual, Res.Df, RSS, Sum.Sq, F.value, p.value, padj, cell)
write.table(myanova_reduced_allCell_forWrite, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/ANOVA_cellSep.tab", row.names=F, col.names=T, quote=F, sep="\t")

mylm_LCL_5df$cell <- "LCL"
mylm_IPSC_5df$cell <- "IPSC"
mylm_CM_5df$cell <- "CM"
mylm_all_5df <- rbind(mylm_LCL_5df, mylm_IPSC_5df, mylm_CM_5df)
write.table(mylm_all_5df, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/FixedEffect_cellSep.tab", row.names=F, col.names=T, quote=F, sep="\t")

vp_allCellTypes_gene_forWrite <- vp_allCellTypes_gene %>% select(SNP_Individual, ensg, g.id, Tr, PlateVar, Control.ID, Residual, CellType)
write.table(vp_allCellTypes_gene_forWrite, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/MixedEffect_cellSep.tab", row.names=F, col.names=T, quote=F, sep="\t")


### Roger wants a table with the number of tests for each variable in the linear model and number below 10% FDR
# For control
mylm_CM_5df_Control <- mylm_CM_5df[grep("Control.IDCO2", mylm_CM_5df$term),]
mylm_CM_5df_Control <- mylm_CM_5df_Control[order(mylm_CM_5df_Control$p.value),]
mylm_CM_5df_Control_anovaASE <- mylm_CM_5df_Control[which(mylm_CM_5df_Control$SNP_Individual %in% CM_anova_ASE),]
mylm_CM_5df_Control_anovaASE$padj <- p.adjust(mylm_CM_5df_Control_anovaASE$p.value, method="BH")
# LCLs
mylm_LCL_5df_Control <- mylm_LCL_5df[grep("Control.IDCO2", mylm_LCL_5df$term),]
mylm_LCL_5df_Control <- mylm_LCL_5df_Control[order(mylm_LCL_5df_Control$p.value),]
mylm_LCL_5df_Control_anovaASE <- mylm_LCL_5df_Control[which(mylm_LCL_5df_Control$SNP_Individual %in% LCL_anova_ASE),]
mylm_LCL_5df_Control_anovaASE$padj <- p.adjust(mylm_LCL_5df_Control_anovaASE$p.value, method="BH")
# IPSCs
mylm_IPSC_5df_Control <- mylm_IPSC_5df[grep("Control.IDCO2", mylm_IPSC_5df$term),]
mylm_IPSC_5df_Control <- mylm_IPSC_5df_Control[order(mylm_IPSC_5df_Control$p.value),]
mylm_IPSC_5df_Control_anovaASE <- mylm_IPSC_5df_Control[which(mylm_IPSC_5df_Control$SNP_Individual %in% IPSC_anova_ASE),]
mylm_IPSC_5df_Control_anovaASE$padj <- p.adjust(mylm_IPSC_5df_Control_anovaASE$p.value, method="BH")

mylm_CM_5df_Control_anovaASE %>% filter(padj < 0.1) %>% pull(SNP_Individual) %>% unique() %>% length() # 24
mylm_CM_5df_Control_anovaASE %>% pull(SNP_Individual) %>% unique() %>% length() # 6811
mylm_IPSC_5df_Control_anovaASE %>% filter(padj < 0.1) %>% pull(SNP_Individual) %>% unique() %>% length() # 55
mylm_IPSC_5df_Control_anovaASE %>% pull(SNP_Individual) %>% unique() %>% length() # 6794
mylm_LCL_5df_Control_anovaASE %>% filter(padj < 0.1) %>% pull(SNP_Individual) %>% unique() %>% length() # 2
mylm_LCL_5df_Control_anovaASE %>% pull(SNP_Individual) %>% unique() %>% length() # 7127

mylm_CM_5df_PlateVar_anovaASE %>% filter(padj < 0.1) %>% pull(SNP_Individual) %>% unique() %>% length() # 68
mylm_CM_5df_PlateVar_anovaASE %>% pull(SNP_Individual) %>% unique() %>% length() # 6811
mylm_IPSC_5df_PlateVar_anovaASE %>% filter(padj < 0.1) %>% pull(SNP_Individual) %>% unique() %>% length() # 98
mylm_IPSC_5df_PlateVar_anovaASE %>% pull(SNP_Individual) %>% unique() %>% length() # 6794
mylm_LCL_5df_PlateVar_anovaASE %>% filter(padj < 0.1) %>% pull(SNP_Individual) %>% unique() %>% length() # 0
mylm_LCL_5df_PlateVar_anovaASE %>% pull(SNP_Individual) %>% unique() %>% length() # 7127

# Treatments
mylm_CM_5df_Tr_anovaASE %>% group_by(term) %>% summarise(n = n())
mylm_CM_5df_Tr_anovaASE %>% group_by(term) %>% filter(padj < 0.1) %>% summarise(n = n())
mylm_IPSC_5df_Tr_anovaASE %>% group_by(term) %>% summarise(n = n())
mylm_IPSC_5df_Tr_anovaASE %>% group_by(term) %>% filter(padj < 0.1) %>% summarise(n = n())
mylm_LCL_5df_Tr_anovaASE %>% group_by(term) %>% summarise(n = n())
mylm_LCL_5df_Tr_anovaASE %>% group_by(term) %>% filter(padj < 0.1) %>% summarise(n = n())


