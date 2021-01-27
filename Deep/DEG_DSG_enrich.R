# Anthony Findley
# 12/10/2020

# Are DEGs enriched for DSGs?
# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/

library(tidyverse)
library(data.table)

# First get number of DS genes per treatment:

#### CM
samplelist <- list(VitaminA = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T6C1_vs_CO2/significant_intron_gene_names.txt",
                   Aldosterone = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T9C1_vs_CO2/significant_intron_gene_names.txt",
                   Dexamethasone = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T12C1_vs_CO2/significant_intron_gene_names.txt",
                   Caffeine = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T13C1_vs_CO1/significant_intron_gene_names.txt",
                   Nicotine = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T14C1_vs_CO1/significant_intron_gene_names.txt",
                   Copper = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T15C1_vs_CO1/significant_intron_gene_names.txt",
                   Selenium = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T19C1_vs_CO1/significant_intron_gene_names.txt",
                   Zinc = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T20C1_vs_CO1/significant_intron_gene_names.txt",
                   Cadmium = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T27C1_vs_CO1/significant_intron_gene_names.txt",
                   Triclosan = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T30C1_vs_CO2/significant_intron_gene_names.txt",
                   Insulin = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T33C1_vs_CO1/significant_intron_gene_names.txt",
                   Acetaminophen = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/CM/CM_T42C1_vs_CO1/significant_intron_gene_names.txt")

samples <- lapply(names(samplelist),function(preprocessing){
  aux <- read.table(samplelist[[preprocessing]], header = T)
})

names(samples) <- names(samplelist)
cm_splice <- map_df(samples, ~as.data.frame(.x), .id="id") %>% mutate(CellType = "CM")


#### IPSC
samplelist <- list(VitaminA = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T6C1_vs_CO2/significant_intron_gene_names.txt",
                   Aldosterone = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T9C1_vs_CO2/significant_intron_gene_names.txt",
                   Dexamethasone = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T12C1_vs_CO2/significant_intron_gene_names.txt",
                   Caffeine = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T13C1_vs_CO1/significant_intron_gene_names.txt",
                   Nicotine = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T14C1_vs_CO1/significant_intron_gene_names.txt",
                   Copper = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T15C1_vs_CO1/significant_intron_gene_names.txt",
                   Selenium = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T19C1_vs_CO1/significant_intron_gene_names.txt",
                   Zinc = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T20C1_vs_CO1/significant_intron_gene_names.txt",
                   Cadmium = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T27C1_vs_CO1/significant_intron_gene_names.txt",
                   Triclosan = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T30C1_vs_CO2/significant_intron_gene_names.txt",
                   Insulin = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T33C1_vs_CO1/significant_intron_gene_names.txt",
                   Acetaminophen = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/IPSC/IPSC_T42C1_vs_CO1/significant_intron_gene_names.txt")

samples <- lapply(names(samplelist),function(preprocessing){
  aux <- read.table(samplelist[[preprocessing]], header = T)
})

names(samples) <- names(samplelist)
ipsc_splice <- map_df(samples, ~as.data.frame(.x), .id="id") %>% mutate(CellType = "IPSC")


#### LCL
samplelist <- list(VitaminA = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T6C1_vs_CO2/significant_intron_gene_names.txt",
                   Aldosterone = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T9C1_vs_CO2/significant_intron_gene_names.txt",
                   Dexamethasone = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T12C1_vs_CO2/significant_intron_gene_names.txt",
                   Caffeine = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T13C1_vs_CO1/significant_intron_gene_names.txt",
                   Nicotine = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T14C1_vs_CO1/significant_intron_gene_names.txt",
                   Copper = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T15C1_vs_CO1/significant_intron_gene_names.txt",
                   Selenium = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T19C1_vs_CO1/significant_intron_gene_names.txt",
                   Zinc = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T20C1_vs_CO1/significant_intron_gene_names.txt",
                   Cadmium = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T27C1_vs_CO1/significant_intron_gene_names.txt",
                   Triclosan = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T30C1_vs_CO2/significant_intron_gene_names.txt",
                   Insulin = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T33C1_vs_CO1/significant_intron_gene_names.txt",
                   Acetaminophen = "/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/LCL/LCL_T42C1_vs_CO1/significant_intron_gene_names.txt")

samples <- lapply(names(samplelist),function(preprocessing){
  aux <- read.table(samplelist[[preprocessing]], header = T)
})

names(samples) <- names(samplelist)
lcl_splice <- map_df(samples, ~as.data.frame(.x), .id="id") %>% mutate(CellType = "LCL")

# Combine all cell types
all_splice <- rbind(cm_splice, ipsc_splice, lcl_splice)

### Now get DE genes

Deep_genes <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/All/DE_genes.txt", header=F, stringsAsFactors=F)
colnames(Deep_genes) <- c("enst", "DE_padj", "DE_pval", "DE_logFC", "ensg", "g.id", "type", "file.path")

# Put Deep_genes in a format which can be used to make plot
Deep_genes$file <- sapply(strsplit(Deep_genes$file.path, "/"), "[", 6) # Useful function that's similar to cut in unix
Deep_genes$CellType <- sapply(strsplit(Deep_genes$file, "_"), "[", 1) # Useful function that's similar to cut in unix
Deep_genes$Treatment.ID <- sapply(strsplit(Deep_genes$file, "_"), "[", 6) # Useful function that's similar to cut in unix

Deep_genes$Treatment.ID <- gsub(".txt", "", Deep_genes$Treatment.ID)
Deep_genes$CellType <- gsub("D", "", Deep_genes$CellType)

# Add Treatment.Name from covariate file
cv_Deep <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/DIPSC1R1_covar.txt", header=T, stringsAsFactors=F)
Deep_genes_merge <- merge(Deep_genes, unique(cv_Deep[,c("Treatment.ID", "Treatment.Name")]), by="Treatment.ID")

# Correct spelling of acetaminophen
Deep_genes_merge[which(Deep_genes_merge$Treatment.Name =="Acetaminophin"), "Treatment.Name"] <- "Acetaminophen"

# Get number of DE and DS genes per treatment and cell type
DEGs <- Deep_genes_merge %>% group_by(CellType, Treatment.Name) %>% count(name = "DEG")
DSGs <- all_splice %>% group_by(CellType, id) %>% count(name = "DSG")

# Get number of genes which are both DEG and DSG
Deep_genes_splice <- inner_join(Deep_genes_merge[,c("Treatment.Name", "CellType", "g.id")], all_splice, by = c("g.id" = "Gene", "Treatment.Name" = "id", "CellType" = "CellType"))
DEG_DSGs <- Deep_genes_splice %>% group_by(CellType, Treatment.Name) %>% count(name = "DEG_DSG")

# Get total number of genes tested in both splicing and GE
all_DE_results <- fread("/nfs/rprdata/Anthony/GxE_iPSC/Deep/All/all_DEgene_results.txt", header=T, stringsAsFactors=F, data.table=F)
colnames(all_DE_results)[8] <- "file.path"

all_DE_results_uniq <- all_DE_results %>% select(file.path, g.id) %>% distinct()

# Put all_DE_results_uniq in a format which can be used to make plot
all_DE_results_uniq$file <- sapply(strsplit(all_DE_results_uniq$file.path, "/"), "[", 6) # Useful function that's similar to cut in unix
all_DE_results_uniq$CellType <- sapply(strsplit(all_DE_results_uniq$file, "_"), "[", 1) # Useful function that's similar to cut in unix
all_DE_results_uniq$Treatment.ID <- sapply(strsplit(all_DE_results_uniq$file, "_"), "[", 6) # Useful function that's similar to cut in unix

all_DE_results_uniq$Treatment.ID <- gsub(".txt", "", all_DE_results_uniq$Treatment.ID)
all_DE_results_uniq$CellType <- gsub("D", "", all_DE_results_uniq$CellType)

# Add Treatment.Name from covariate file
all_DE_results_uniq_merge <- merge(all_DE_results_uniq[,c("CellType", "Treatment.ID", "g.id")], unique(cv_Deep[,c("Treatment.ID", "Treatment.Name")]), by="Treatment.ID")

### Get tested splicing genes
Alan_LCL <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/LCL_tested_genes_unique.txt", stringsAsFactors=F, header=F)
Alan_LCL$CellType <- "LCL"

Alan_CM <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/CM_tested_genes_unique.txt", stringsAsFactors=F, header=F)
Alan_CM$CellType <- "CM"

Alan_IPSC <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/VariancePartitioning/IPSC_tested_genes_unique.txt", stringsAsFactors=F, header=F)
Alan_IPSC$CellType <- "IPSC"

Alan_all <- rbind(Alan_LCL, Alan_IPSC, Alan_CM)

# Join genes tested for splicing and expression to see how many per treatment
GE_splice_test <- inner_join(all_DE_results_uniq_merge, Alan_all, by = c("g.id" = "V1", "CellType" = "CellType"))
GE_splice_test_count <- GE_splice_test %>% group_by(CellType, Treatment.Name) %>% count(name = "Tested_genes")

# Subset DE genes to just those tested for splicing
Deep_genes_merge %>% select(CellType, Treatment.Name, g.id) %>% distinct() %>% inner_join(Alan_all, by = c("g.id" = "V1", "CellType" = "CellType")) -> DE_genes_tested_for_splice
DEGs <- DE_genes_tested_for_splice %>% group_by(CellType, Treatment.Name) %>% count(name = "DEG")

# Subset DS genes to just those tested for expression
all_splice %>% inner_join(all_DE_results_uniq_merge, by = c("Gene" = "g.id", "id" = "Treatment.Name", "CellType" = "CellType")) -> all_splice_tested_for_GE
DSGs <- all_splice_tested_for_GE %>% group_by(CellType, id) %>% count(name = "DSG")


# Combine all the data
df_1 <- inner_join(DEGs, DSGs, by = c("CellType" = "CellType", "Treatment.Name" = "id"))
df_2 <- inner_join(df_1, DEG_DSGs)
final_df <- inner_join(df_2, GE_splice_test_count)

# Subtract genes both DEG and DSG (DEG_DSG) from DEG and DSG category to get number of genes only DEG/DSG, and subtract resulting from all tested genes to see how many genes are not in any category
final_df %>% mutate(DEG_only = DEG - DEG_DSG, DSG_only = DSG - DEG_DSG, Not_DEG_DSG = Tested_genes - DEG_DSG - DEG_only - DSG_only) -> final_df



fisher_test_row <- function(data){
	enrich <- matrix(c(data$DEG_DSG, data$DSG_only, data$DEG_only, data$Not_DEG_DSG),
		nrow = 2,
		dimnames = list(Ancestry = c("Nean", "nonNean"), Type = c("centiSNP", "Footprint")))
	fish_test <- fisher.test(enrich, alternative = "two.sided")
	names(fish_test$estimate) <- "Odds_Ratio"
	final_results <- c(fish_test$estimate, Lower = fish_test$conf.int[1], Upper = fish_test$conf.int[2], Pval = fish_test$p.value, DEG_DSG = data$DEG_DSG, DSG_only = data$DSG_only, DEG_only = data$DEG_only, Not_DEG_DSG = data$Not_DEG_DSG)
}

aux <- final_df %>% group_by(CellType, Treatment.Name) %>% nest()  
mylist <- aux$data
names(mylist) <- paste0(aux$CellType, "_", aux$Treatment.Name)
myfisher <- map_dfr(mylist,fisher_test_row)

myfisher$Sample <- names(mylist)
myfisher$padj <- p.adjust(myfisher$Pval, method="BH")

myfisher$CellType <- sapply(strsplit(myfisher$Sample, "_"), "[", 1) # Useful function that's similar to cut in unix
myfisher$Treatment <- sapply(strsplit(myfisher$Sample, "_"), "[", 2) # Useful function that's similar to cut in unix



fp <- ggplot(data=myfisher, aes(x=Treatment, y=log(Odds_Ratio), ymin=log(Lower), ymax=log(Upper))) +
geom_pointrange(size=0.8) +
geom_hline(yintercept=log(1), lty=2) +
coord_flip() +  # flip coordinates (puts labels on y axis)
xlab("Cell Type") + ylab("Enrichment (95% CI); log-transformed") + theme_bw() + geom_errorbar(aes(ymin=log(Lower), ymax=log(Upper)), width=0.5) + facet_grid(. ~ CellType) + theme(text = element_text(size=14, color="black"), legend.position = "none")

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/plots/DEG_DSG_enrich.fp.pdf", height=5, width=10)
fp
dev.off()
