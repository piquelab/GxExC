### Make ASE barplot

# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/fixed_effect

R

library(data.table)
library(tidyverse)
library(RColorBrewer)

ASE <- fread("../../Quasar_output/all_allOutput_noHeaders.sorted.bed", header=F, data.table=F)
colnames(ASE) <- c("chr", "pos0", "pos", "ref", "alt", "rsID", "af", "cell.line", "treatment", "ref.reads", "alt.reads", "beta", "beta.se", "pval", "qval")

# Remove sex chromsomes and add 1 to every ref and alt read and calculate log ratio
ASE %>% filter(chr != "X" & chr != "Y") -> ASE

cv <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/all_covar.txt", header=T, stringsAsFactors=F)

ASE <- merge(ASE, cv, by.x = "treatment", by.y = "Filename")

# For DCM1R2, change water control to ethanol
ASE <- ASE[!(ASE$Plate.ID == "DCM1R2" & ASE$Treatment.Name == "Water"),]

# Take only the significant ASE (q-value < 0.1)
ASE_sig <- ASE[which(ASE$qval < 0.1),]

# Aggregate data so it's in the format needed for ggplot to make a barplot
for_plot <- aggregate(cbind(count = rsID) ~ CellType + Treatment.Name, data=ASE_sig, FUN = function(x){NROW(x)})

for_plot$Treatment.Name <- gsub("Acetaminophin", "Acetaminophen", for_plot$Treatment.Name)

for_plot$CellType <- factor(for_plot$CellType, levels=c("LCL", "IPSC", "CM"))

# Treatments are currently in alphabetical order. I want to keep that order, but put ethanol and water at the end.
library(forcats)
for_plot$Treatment.Name <- factor(for_plot$Treatment.Name, levels=unique(for_plot$Treatment.Name))
for_plot$Treatment.Name <- forcats::fct_relevel(for_plot$Treatment.Name, "Ethanol", after = Inf)
for_plot$Treatment.Name <- forcats::fct_relevel(for_plot$Treatment.Name, "Water", after = Inf)

# Put for_plot in same order as it will be plotted so colors will line up.
for_plot <- for_plot[order(for_plot$CellType, for_plot$Treatment.Name),]

p <- ggplot(data=for_plot, aes(x=Treatment.Name, y=count, color=CellType, fill=CellType)) +
  geom_bar(stat="identity") + facet_grid(CellType ~ .) +
  theme_bw() + theme(axis.text.x = element_text(angle = 60, vjust=0.5)) + scale_color_manual(values=c("#C77CFF", "#619CFF", "#00BA38")) + scale_fill_manual(values=c("#C77CFF", "#619CFF", "#00BA38"))

pdf("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/QuASAR/ASE_Barplot_CellColors.pdf")
p
dev.off()
