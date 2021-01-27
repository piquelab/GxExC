# Do GTEx overlap with cASE genes. Also, is there enrichment of cASE genes in loss of function genes? This has been updated to use the data which excludes sex chromosomes.

# First I need to get all genes which are an eGene in any tissue:




# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/annotation

less /wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/*egenes* | awk '$29 < 0.05' | cut -f1 | sort | uniq > all_tissues.sig_eGenes.txt # 34,548 significant eGenes in any tissue

less /wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/*egenes* | cut -f1 | sort | uniq > all_tissues.tested_eGenes.txt # 38,833 genes tested for being eGenes


#############################
#############################
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
library(upsetR)
library(lme4)

load("interaction_noReadCov_DF_ANOVA_add1_noX/anova_lm_vp.Rd")

# From model with all cell types:
# Genes with ANOVA ASE: myanova_all_5df_finite_gene %>% filter(padj < 0.1) %>% pull(ensg) %>% unique()

# Genes tested by GTEx:
GTEx_tested <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/annotation/all_tissues.tested_eGenes.txt", header=F, stringsAsFactors=F)
# Remove the decimal: 
GTEx_tested <- gsub("\\..*", "", GTEx_tested$V1)

# GTEx sig
GTEx_eGenes <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/annotation/all_tissues.sig_eGenes.txt", header=F, stringsAsFactors=F)
# Remove the decimal: 
GTEx_eGenes <- gsub("\\..*", "", GTEx_eGenes$V1)

# ANOVA ASE tested
gene_anno <- read.table("../../Quasar_output/ASE_SNPs.genes.uniq.txt", header=F, stringsAsFactors=F)
colnames(gene_anno) <- c("chr", "pos0", "pos1", "ref", "alt", "rsID", "ensg", "gene_type", "g.id")

myanova_all_5df_finite$rsID <- sapply(strsplit(myanova_all_5df_finite$SNP_Individual, "_"), "[", 1)
myanova_all_5df_finite_gene <- merge(myanova_all_5df_finite, gene_anno, by="rsID")

myanova_all_5df_finite_gene %>% pull(ensg) %>% unique() -> ANOVA_ASE_tested_ensg

# ANOVA ASE sig
myanova_all_5df_finite_gene %>% filter(padj < 0.1) %>% pull(ensg) %>% unique() -> ANOVA_ASE_sig_ensg
length(ANOVA_ASE_sig_ensg)
[1] 5646

# Only keep genes tested in both:
GTEx_tested <- GTEx_tested[(GTEx_tested %in% ANOVA_ASE_tested_ensg)]
ANOVA_ASE_tested_ensg <- ANOVA_ASE_tested_ensg[(ANOVA_ASE_tested_ensg %in% GTEx_tested)]

ANOVA_ASE_sig_ensg <- ANOVA_ASE_sig_ensg[(ANOVA_ASE_sig_ensg %in% GTEx_tested)]
GTEx_eGenes <- GTEx_eGenes[(GTEx_eGenes %in% ANOVA_ASE_tested_ensg)]

###
# PROBLEM:
# Of the genes which were tested for ASE and for being an eGene (9,594), 99% (9,516) are eGenes. So I can't do the enrichment by considering if a gene is an eGene in any tissue. Try enrichment for ASE/eGene for CM tissue only.

GTEx_LV <- read.table("/wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/Heart_Left_Ventricle.v8.egenes.txt.gz", header=T, stringsAsFactors=F, sep="\t")

GTEx_LV_tested <- unique(GTEx_LV$gene_id)
GTEx_LV_tested <- unique(gsub("\\..*", "", GTEx_LV_tested)) # 21,353 genes

GTEx_LV %>% filter(qval < 0.05) %>% pull(gene_id) -> GTEx_LV_eGenes 
GTEx_LV_eGenes <- unique(gsub("\\..*", "", GTEx_LV_eGenes)) # 9,642 genes

# Only keep genes tested in both:
GTEx_LV_tested <- GTEx_LV_tested[(GTEx_LV_tested %in% ANOVA_ASE_tested_ensg)]
ANOVA_ASE_tested_LV_ensg <- ANOVA_ASE_tested_ensg[(ANOVA_ASE_tested_ensg %in% GTEx_LV_tested)]

ANOVA_ASE_LV_sig_ensg <- ANOVA_ASE_sig_ensg[(ANOVA_ASE_sig_ensg %in% GTEx_LV_tested)]
GTEx_LV_eGenes <- GTEx_LV_eGenes[(GTEx_LV_eGenes %in% ANOVA_ASE_tested_LV_ensg)]

# Do Fisher's exact test
ASE_P_GTEx_P <- sum(GTEx_LV_eGenes %in% ANOVA_ASE_LV_sig_ensg)
ASE_P_GTEx_N <- sum(! ANOVA_ASE_LV_sig_ensg %in% GTEx_LV_eGenes)

ASE_N_GTEx_P <- sum(! GTEx_LV_eGenes %in% ANOVA_ASE_LV_sig_ensg)
ASE_N_GTEx_N <- length(GTEx_LV_tested) - ASE_P_GTEx_P - ASE_P_GTEx_N - ASE_N_GTEx_P

# Do fisher's test to see if cASE is enriched for DE
ASE_GTEx_LV_enrich <-
matrix(c(ASE_P_GTEx_P, ASE_N_GTEx_P, ASE_P_GTEx_N, ASE_N_GTEx_N),
       nrow = 2,
       dimnames = list(ASE = c("ASE+", "ASE-"),
                       GTEx_LV = c("GTEx+", "GTEx-")))
fish_test_ASE_GTEx_LV <- fisher.test(ASE_GTEx_LV_enrich, alternative = "two.sided")

#########################
#########################
# Now see if treatment cASE is more enriched.

ASE <- fread("../../Quasar_output/all_allOutput_noHeaders.sorted.bed", header=F)
colnames(ASE) <- c("chr", "pos0", "pos", "ref", "alt", "rsID", "af", "cell.line", "treatment", "ref.reads", "alt.reads", "beta", "beta.se", "pval", "qval")
ASE <- data.frame(ASE)

# Add 1 to every ref and alt read and calculate log ratio
ASE$ref.reads1 <- ASE$ref.reads + 1
ASE$alt.reads1 <- ASE$alt.reads + 1
ASE$beta1 <- log(ASE$ref.reads1 / ASE$alt.reads1)

cv <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/all_covar.txt", header=T, stringsAsFactors=F)

ASE <- merge(ASE, cv, by.x = "treatment", by.y = "Filename")

# For DCM1R2, change water control to ethanol
ASE <- ASE[!(ASE$Plate.ID == "DCM1R2" & ASE$Treatment.Name == "Water"),]
ASE[(ASE$Plate.ID == "DCM1R2" & ASE$Control.ID == "CO1"), "Control.ID"] <- "CO2"
ASE$SNP_Individual <- paste0(ASE$rsID, "_", ASE$Individual)

# Get cASE SNPs
# Only keep SNPs which have at least one water and one ethanol control for each cell type (except DCM1R2, which needs just one ethanol)

ASE_etoh <- ASE[ASE$Treatment.Name == "Ethanol",]
ASE_etoh_cast <- dcast(ASE_etoh, SNP_Individual ~ CellType) # Creates df of # of etoh per cell type
good_ASE_etoh <- ASE_etoh_cast[ASE_etoh_cast$CM > 0 & ASE_etoh_cast$IPSC > 0 & ASE_etoh_cast$LCL > 0, "SNP_Individual"]
length(good_ASE_etoh)
[1] 69133

ASE_h2o <- ASE[ASE$Treatment.Name == "Water",]
ASE_h2o_cast <- dcast(ASE_h2o, SNP_Individual ~ CellType) # Creates df of # of h2o per cell type
good_ASE_h2o <- ASE_h2o_cast[ASE_h2o_cast$CM > 0 & ASE_h2o_cast$IPSC > 0 & ASE_h2o_cast$LCL > 0, "SNP_Individual"]
length(good_ASE_h2o)
[1] 67163 # was 69882

all_anova_ASE <- unique(myanova_all_5df_finite[myanova_all_5df_finite$padj < 0.1, "SNP_Individual"])

mylm_all_5df %>% filter(SNP_Individual %in% all_anova_ASE, grepl("^Tr", term)) %>% mutate(group = "Treatment Effect", padj = p.adjust(p.value, method="BH")) %>% filter(padj < 0.1 & SNP_Individual %in% good_ASE_h2o & SNP_Individual %in% good_ASE_etoh) -> Tr_cASE_etoh_h2o

dim(Tr_cASE_etoh_h2o)
[1] 1409   10  # 1102 unique SNP_Ind

Tr_cASE_etoh_h2o$rsID <- sapply(strsplit(Tr_cASE_etoh_h2o$SNP_Individual, "_"), "[", 1)
Tr_cASE_etoh_h2o_gene <- merge(Tr_cASE_etoh_h2o, gene_anno, by="rsID")

### I'm going to try the enrichment 2 ways. First, I'll include all genes tested for ASE. Then, I'll subset to just the genes tested for cASE (which means they had to be ASE +).

# Only keep genes tested in both:
GTEx_LV_tested <- GTEx_LV_tested[(GTEx_LV_tested %in% ANOVA_ASE_tested_ensg)]
ANOVA_ASE_tested_LV_ensg <- ANOVA_ASE_tested_ensg[(ANOVA_ASE_tested_ensg %in% GTEx_LV_tested)]

cASE_sig_ensg <- Tr_cASE_etoh_h2o_gene %>% filter(ensg %in% GTEx_LV_tested) %>% pull(ensg) %>% unique() # 881 genes; was 900 with X
GTEx_LV_eGenes <- GTEx_LV_eGenes[(GTEx_LV_eGenes %in% ANOVA_ASE_tested_LV_ensg)] 

# Do Fisher's exact test
cASE_P_GTEx_P <- sum(GTEx_LV_eGenes %in% cASE_sig_ensg)
cASE_P_GTEx_N <- sum(! cASE_sig_ensg %in% GTEx_LV_eGenes)

cASE_N_GTEx_P <- sum(! GTEx_LV_eGenes %in% cASE_sig_ensg)
cASE_N_GTEx_N <- length(GTEx_LV_tested) - cASE_P_GTEx_P - cASE_P_GTEx_N - cASE_N_GTEx_P

# Do fisher's test to see if cASE is enriched for DE
cASE_GTEx_LV_enrich <-
matrix(c(cASE_P_GTEx_P, cASE_N_GTEx_P, cASE_P_GTEx_N, cASE_N_GTEx_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       GTEx_LV = c("GTEx+", "GTEx-")))
fish_test_cASE_GTEx_LV <- fisher.test(cASE_GTEx_LV_enrich, alternative = "two.sided")

# Now subset to just those genes which were positive for ASE and tested for cASE

GTEx_LV_tested_ASE <- GTEx_LV_tested[(GTEx_LV_tested %in% ANOVA_ASE_LV_sig_ensg)]
# ANOVA_ASE_tested_LV_ensg <- ANOVA_ASE_tested_ensg[(ANOVA_ASE_tested_ensg %in% GTEx_LV_tested)] # This is now just ANOVA_ASE_LV_sig_ensg

# cASE_sig_ensg <- Tr_cASE_etoh_h2o_gene %>% filter(ensg %in% GTEx_LV_tested) %>% pull(ensg) %>% unique() # 900 genes. This stays the same as before
GTEx_LV_eGenes_ASE <- GTEx_LV_eGenes[(GTEx_LV_eGenes %in% ANOVA_ASE_LV_sig_ensg)] 

# Do Fisher's exact test. This overwrites previous
cASE_P_GTEx_P <- sum(GTEx_LV_eGenes_ASE %in% cASE_sig_ensg)
cASE_P_GTEx_N <- sum(! cASE_sig_ensg %in% GTEx_LV_eGenes_ASE)

cASE_N_GTEx_P <- sum(! GTEx_LV_eGenes_ASE %in% cASE_sig_ensg)
cASE_N_GTEx_N <- length(GTEx_LV_tested_ASE) - cASE_P_GTEx_P - cASE_P_GTEx_N - cASE_N_GTEx_P

# Do fisher's test to see if ccASE is enriched for DE
cASE_GTEx_LV_enrich_ASEonly <-
matrix(c(cASE_P_GTEx_P, cASE_N_GTEx_P, cASE_P_GTEx_N, cASE_N_GTEx_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       GTEx_LV = c("GTEx+", "GTEx-")))
fish_test_cASE_GTEx_LV_ASEonly <- fisher.test(cASE_GTEx_LV_enrich_ASEonly, alternative = "two.sided")

#############
# Are treatment cASE genes enriched in loss of function-tolerant genes?
gnomad <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/annotation/gnomad.v2.1.1.lof_metrics.by_gene.txt", header=T, stringsAsFactors=F, sep="\t")

# gnomAD genes tested for ASE, with column for if it's ASE or cASE:
myanova_all_5df_finite_gene %>% pull(g.id) %>% unique() -> ANOVA_ASE_tested_gid
# ANOVA ASE sig
myanova_all_5df_finite_gene %>% filter(padj < 0.1) %>% pull(g.id) %>% unique() -> ANOVA_ASE_sig_gid


gnomad_testedASE <- gnomad %>% filter(gene %in% ANOVA_ASE_tested_gid) %>% mutate(ASE = if_else(gene %in% ANOVA_ASE_sig_gid, 1, 0), cASE = if_else(gene %in% Tr_cASE_etoh_h2o_gene$g.id, 1, 0))

# Plot oe_lof_upper for ASE/no ASE and cASE/no cASE
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/LOEUF.ASE.box.pdf")
ggplot(gnomad_testedASE, aes(x=factor(ASE), y=as.numeric(oe_lof_upper), fill=factor(ASE))) + geom_boxplot(notch=T) + theme_bw() + labs(title="LOEUF values by ASE status (All cell types") + xlab("ASE") + ylab("LOEUF") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) #+ ylim(-12,3)
dev.off()

# Check if it's significant (it's not):
summary(lm(oe_lof_upper ~  factor(ASE), data=gnomad_testedASE))
Coefficients:
             Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.776995   0.008131  95.564   <2e-16 ***
factor(ASE)1  0.02000    0.01075   1.861   0.0627 .


# Now do cASE
pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/all_cell/DF_cutoff/ANOVA_add1/LOEUF.cASE.box.pdf")
ggplot(gnomad_testedASE, aes(x=factor(cASE), y=as.numeric(oe_lof_upper), fill=factor(cASE))) + geom_boxplot(notch=T) + theme_bw() + labs(title="LOEUF values by cASE status (All cell types") + xlab("cASE") + ylab("LOEUF") + theme(text = element_text(size=12), axis.text.x = element_text(size=15), axis.text.y = element_text(size=15)) #+ ylim(-12,3)
dev.off()

# Check if it's significant (it's not):
summary(lm(oe_lof_upper ~  factor(cASE), data=gnomad_testedASE))
Coefficients:
              Estimate Std. Error t value Pr(>|t|)
(Intercept)    0.79060    0.00562 140.670   <2e-16 ***
factor(cASE)1 -0.02064    0.01748  -1.181    0.238


#########################
#########################
# Now do CM ASE/cASE only:

load("Cell_Type_Specific_noReadCov_DF_ANOVA_add1_noX/anova_lm_vp.Rd") # c("ASE", "cv", "myanova_reduced_allCell", "mylm_CM_5df", "mylm_IPSC_5df", "mylm_LCL_5df", "vp_allCellTypes_gene")

# Get CM ASE
myanova_reduced_allCell %>% filter(cell == "CM") -> myanova_reduced_CM

# Only keep SNPs with 5 degrees of freedom from linear model.
myanova_reduced_CM_5df <- myanova_reduced_CM[ myanova_reduced_CM$SNP_Individual %in% unique(mylm_CM_5df$SNP_Individual), ]
myanova_reduced_CM_5df_finite <- myanova_reduced_CM_5df[ is.finite(myanova_reduced_CM_5df$p.value),]
myanova_reduced_CM_5df_finite <- myanova_reduced_CM_5df_finite[order(myanova_reduced_CM_5df_finite$p.value),]
myanova_reduced_CM_5df_finite$padj <- p.adjust(myanova_reduced_CM_5df_finite$p.value, method="BH")
sum(myanova_reduced_CM_5df_finite$padj < 0.1) # How much ASE with FDR < 0.1; 6811

myanova_reduced_CM_5df_finite$rsID <- sapply(strsplit(myanova_reduced_CM_5df_finite$SNP_Individual, "_"), "[", 1)
myanova_reduced_CM_5df_finite_gene <- merge(myanova_reduced_CM_5df_finite, gene_anno, by="rsID")


# I need to reload the GTEx data because I overwrote it when I was subsetting to genes present in both.
GTEx_LV <- read.table("/wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/Heart_Left_Ventricle.v8.egenes.txt.gz", header=T, stringsAsFactors=F, sep="\t")

GTEx_LV_tested <- unique(GTEx_LV$gene_id)
GTEx_LV_tested <- unique(gsub("\\..*", "", GTEx_LV_tested)) # 21,353 genes

GTEx_LV %>% filter(qval < 0.05) %>% pull(gene_id) -> GTEx_LV_eGenes 
GTEx_LV_eGenes <- unique(gsub("\\..*", "", GTEx_LV_eGenes)) # 9,642 genes

# Now get left atrium
GTEx_AA <- read.table("/wsu/home/groups/piquelab/GTEx_official_data/V8/single_tissue_results/GTEx_Analysis_v8_eQTL/Heart_Atrial_Appendage.v8.egenes.txt.gz", header=T, stringsAsFactors=F, sep="\t")

GTEx_AA_tested <- unique(GTEx_AA$gene_id)
GTEx_AA_tested <- unique(gsub("\\..*", "", GTEx_AA_tested)) # 21,353 genes

GTEx_AA %>% filter(qval < 0.05) %>% pull(gene_id) -> GTEx_AA_eGenes 
GTEx_AA_eGenes <- unique(gsub("\\..*", "", GTEx_AA_eGenes)) # 9,642 genes

# ANOVA ASE tested
myanova_reduced_CM_5df_finite_gene %>% pull(ensg) %>% unique() -> CM_ASE_tested_ensg

# ANOVA ASE sig
myanova_reduced_CM_5df_finite_gene %>% filter(padj < 0.1) %>% pull(ensg) %>% unique() -> CM_ASE_sig_ensg
length(CM_ASE_sig_ensg) # 3033; was 3098 with X

# How many ASE genes in LV?
sum(CM_ASE_sig_ensg %in% GTEx_LV_eGenes) # 1519; was 1550

# How many ASE genes in AA?
sum(CM_ASE_sig_ensg %in% GTEx_AA_eGenes) # 1619; was 1651

# For enrichment, only keep genes tested in both. Start with LV:
GTEx_LV_tested_ASE <- GTEx_LV_tested[(GTEx_LV_tested %in% CM_ASE_tested_ensg)]
CM_ASE_tested_ensg_LV <- CM_ASE_tested_ensg[(CM_ASE_tested_ensg %in% GTEx_LV_tested)]

CM_ASE_sig_ensg_LV <- CM_ASE_sig_ensg[(CM_ASE_sig_ensg %in% GTEx_LV_tested_ASE)]
GTEx_LV_eGenes_ASE <- GTEx_LV_eGenes[(GTEx_LV_eGenes %in% CM_ASE_tested_ensg_LV)]

# Do Fisher's exact test
CM_P_LV_P <- sum(GTEx_LV_eGenes_ASE %in% CM_ASE_sig_ensg_LV)
CM_P_LV_N <- sum(! CM_ASE_sig_ensg_LV %in% GTEx_LV_eGenes_ASE)

CM_N_LV_P <- sum(! GTEx_LV_eGenes_ASE %in% CM_ASE_sig_ensg_LV)
CM_N_LV_N <- length(GTEx_LV_tested_ASE) - CM_P_LV_P - CM_P_LV_N - CM_N_LV_P

# Do fisher's test to see if ASE is enriched for LV eQTLs
CM_ASE_GTEx_LV_enrich <-
matrix(c(CM_P_LV_P, CM_N_LV_P, CM_P_LV_N, CM_N_LV_N),
       nrow = 2,
       dimnames = list(ASE = c("ASE+", "ASE-"),
                       GTEx_LV = c("GTEx+", "GTEx-")))
fish_test_CM_ASE_GTEx_LV <- fisher.test(CM_ASE_GTEx_LV_enrich, alternative = "two.sided")

### Now do AA
GTEx_AA_tested_ASE <- GTEx_AA_tested[(GTEx_AA_tested %in% CM_ASE_tested_ensg)]
CM_ASE_tested_ensg_AA <- CM_ASE_tested_ensg[(CM_ASE_tested_ensg %in% GTEx_AA_tested)]

CM_ASE_sig_ensg_AA <- CM_ASE_sig_ensg[(CM_ASE_sig_ensg %in% GTEx_AA_tested_ASE)]
GTEx_AA_eGenes_ASE <- GTEx_AA_eGenes[(GTEx_AA_eGenes %in% CM_ASE_tested_ensg_AA)]

# Do Fisher's exact test
CM_P_AA_P <- sum(GTEx_AA_eGenes_ASE %in% CM_ASE_sig_ensg_AA)
CM_P_AA_N <- sum(! CM_ASE_sig_ensg_AA %in% GTEx_AA_eGenes_ASE)

CM_N_AA_P <- sum(! GTEx_AA_eGenes_ASE %in% CM_ASE_sig_ensg_AA)
CM_N_AA_N <- length(GTEx_AA_tested_ASE) - CM_P_AA_P - CM_P_AA_N - CM_N_AA_P

# Do fisher's test to see if ASE is enriched for AA eQTLs
CM_ASE_GTEx_AA_enrich <-
matrix(c(CM_P_AA_P, CM_N_AA_P, CM_P_AA_N, CM_N_AA_N),
       nrow = 2,
       dimnames = list(ASE = c("ASE+", "ASE-"),
                       GTEx_AA = c("GTEx+", "GTEx-")))
fish_test_CM_ASE_GTEx_AA <- fisher.test(CM_ASE_GTEx_AA_enrich, alternative = "two.sided")

############
############
### Now do CM cASE

CM_anova_ASE <- unique(myanova_reduced_CM_5df_finite[myanova_reduced_CM_5df_finite$padj < 0.1, "SNP_Individual"])
mylm_CM_5df_Tr <- mylm_CM_5df[grep("Tr", mylm_CM_5df$term),]
mylm_CM_5df_Tr <- mylm_CM_5df_Tr[order(mylm_CM_5df_Tr$p.value),]

# Subset cASE for SNPs which display ASE
mylm_CM_5df_Tr_anovaASE <- mylm_CM_5df_Tr[which(mylm_CM_5df_Tr$SNP_Individual %in% CM_anova_ASE),]
mylm_CM_5df_Tr_anovaASE$padj <- p.adjust(mylm_CM_5df_Tr_anovaASE$p.value, method="BH")

# See how many significant interactions and SNP_Individuals
sum(mylm_CM_5df_Tr_anovaASE$padj < 0.1) # 492; was 441
length(unique(mylm_CM_5df_Tr_anovaASE[ mylm_CM_5df_Tr_anovaASE$padj < 0.1, "SNP_Individual"])) # 352; was 320

# Add gene information
mylm_CM_5df_Tr_anovaASE$rsID <- sapply(strsplit(mylm_CM_5df_Tr_anovaASE$SNP_Individual, "_"), "[", 1)
mylm_CM_5df_Tr_anovaASE_gene <- merge(mylm_CM_5df_Tr_anovaASE, gene_anno, by="rsID")

# cASE tested
mylm_CM_5df_Tr_anovaASE_gene %>% pull(ensg) %>% unique() -> CM_cASE_tested_ensg

# cASE sig
mylm_CM_5df_Tr_anovaASE_gene %>% filter(padj < 0.1) %>% pull(ensg) %>% unique() -> CM_cASE_sig_ensg
length(CM_cASE_sig_ensg) # 338; was 308

# How many cASE genes in LV?
sum(CM_cASE_sig_ensg %in% GTEx_LV_eGenes) # 170

# How many ASE genes in AA?
sum(CM_cASE_sig_ensg %in% GTEx_AA_eGenes) # 178

# For enrichment, only keep genes tested in both. Start with LV:
GTEx_LV_tested_cASE <- GTEx_LV_tested[(GTEx_LV_tested %in% CM_cASE_tested_ensg)]
CM_cASE_tested_ensg_LV <- CM_cASE_tested_ensg[(CM_cASE_tested_ensg %in% GTEx_LV_tested)]

CM_cASE_sig_ensg_LV <- CM_cASE_sig_ensg[(CM_cASE_sig_ensg %in% GTEx_LV_tested_cASE)]
GTEx_LV_eGenes_cASE <- GTEx_LV_eGenes[(GTEx_LV_eGenes %in% CM_cASE_tested_ensg_LV)]

# Do Fisher's exact test
CM_P_LV_P <- sum(GTEx_LV_eGenes_cASE %in% CM_cASE_sig_ensg_LV)
CM_P_LV_N <- sum(! CM_cASE_sig_ensg_LV %in% GTEx_LV_eGenes_cASE)

CM_N_LV_P <- sum(! GTEx_LV_eGenes_cASE %in% CM_cASE_sig_ensg_LV)
CM_N_LV_N <- length(GTEx_LV_tested_cASE) - CM_P_LV_P - CM_P_LV_N - CM_N_LV_P

# Do fisher's test to see if ccASE is enriched for DE
CM_cASE_GTEx_LV_enrich <-
matrix(c(CM_P_LV_P, CM_N_LV_P, CM_P_LV_N, CM_N_LV_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       GTEx_LV = c("GTEx+", "GTEx-")))
fish_test_CM_cASE_GTEx_LV <- fisher.test(CM_cASE_GTEx_LV_enrich, alternative = "two.sided")
p-value = 0.901
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.7943748 1.3151661
sample estimates:
odds ratio
   1.02115


### Now do AA
GTEx_AA_tested_cASE <- GTEx_AA_tested[(GTEx_AA_tested %in% CM_cASE_tested_ensg)]
CM_cASE_tested_ensg_AA <- CM_cASE_tested_ensg[(CM_cASE_tested_ensg %in% GTEx_AA_tested)]

CM_cASE_sig_ensg_AA <- CM_cASE_sig_ensg[(CM_cASE_sig_ensg %in% GTEx_AA_tested_cASE)]
GTEx_AA_eGenes_cASE <- GTEx_AA_eGenes[(GTEx_AA_eGenes %in% CM_cASE_tested_ensg_AA)]

# Do Fisher's exact test
CM_P_AA_P <- sum(GTEx_AA_eGenes_cASE %in% CM_cASE_sig_ensg_AA)
CM_P_AA_N <- sum(! CM_cASE_sig_ensg_AA %in% GTEx_AA_eGenes_cASE)

CM_N_AA_P <- sum(! GTEx_AA_eGenes_cASE %in% CM_cASE_sig_ensg_AA)
CM_N_AA_N <- length(GTEx_AA_tested_cASE) - CM_P_AA_P - CM_P_AA_N - CM_N_AA_P

# Do fisher's test to see if ccASE is enriched for DE
CM_cASE_GTEx_AA_enrich <-
matrix(c(CM_P_AA_P, CM_N_AA_P, CM_P_AA_N, CM_N_AA_N),
       nrow = 2,
       dimnames = list(cASE = c("cASE+", "cASE-"),
                       GTEx_AA = c("GTEx+", "GTEx-")))
fish_test_CM_cASE_GTEx_AA <- fisher.test(CM_cASE_GTEx_AA_enrich, alternative = "two.sided")
p-value = 0.712
alternative hypothesis: true odds ratio is not equal to 1
95 percent confidence interval:
 0.7441638 1.2223232
sample estimates:
odds ratio
 0.9528136


############
### Plot CM only results
# Put data together and make a forest plot for up vs down, GE only
for_plot_CM <- data.frame(ASE_cASE = c("ASE", "ASE", "cASE", "cASE"), Tissue = c("LV", "AA", "LV", "AA"), OR = c(fish_test_CM_ASE_GTEx_LV$estimate, fish_test_CM_ASE_GTEx_AA$estimate, fish_test_CM_cASE_GTEx_LV$estimate, fish_test_CM_cASE_GTEx_AA$estimate), Low = c(fish_test_CM_ASE_GTEx_LV$conf.int[1], fish_test_CM_ASE_GTEx_AA$conf.int[1], fish_test_CM_cASE_GTEx_LV$conf.int[1], fish_test_CM_cASE_GTEx_AA$conf.int[1]), Hi = c(fish_test_CM_ASE_GTEx_LV$conf.int[2], fish_test_CM_ASE_GTEx_AA$conf.int[2], fish_test_CM_cASE_GTEx_LV$conf.int[2], fish_test_CM_cASE_GTEx_AA$conf.int[2]), pval = c(fish_test_CM_ASE_GTEx_LV$p.value, fish_test_CM_ASE_GTEx_AA$p.value, fish_test_CM_cASE_GTEx_LV$p.value, fish_test_CM_cASE_GTEx_AA$p.value))

for_plot_CM$sampleID <- paste0(for_plot_CM$ASE_cASE, "_", for_plot_CM$Tissue)

for_plot_CM <- for_plot_CM %>% arrange(desc(ASE_cASE), Tissue)
for_plot_CM$sampleID <- factor(for_plot_CM$sampleID, levels=for_plot_CM$sampleID)

# Plot
fp <- ggplot(data=for_plot_CM, aes(x=sampleID, y=OR, ymin=Low, ymax=Hi)) + geom_errorbar(aes(ymin=Low, ymax=Hi), width=0.3, size=1.2) +
geom_pointrange(aes(fill=ASE_cASE), shape=21, color="black", size=1) +
geom_hline(yintercept=1, lty=2) +
coord_flip() +  # flip coordinates (puts labels on y axis)
xlab("Analysis") + ylab("Enrichment (95% CI)") + theme_bw() + theme(text = element_text(size=14, color="black")) + scale_fill_manual(values = c("Red", "Blue"))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/cASE_Linear_Model/models/cell_type_specific/DF_cutoff/ANOVA_add1/ASE_GTEx_enrichment.forest.noX.pdf", height=2, width=6)
fp
dev.off()
