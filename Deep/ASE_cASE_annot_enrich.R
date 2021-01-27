# Anthony Findley
# 11/9/2020

# This is Alan's script for doing the enrichment of ASE/cASE SNPs in genomic annotations. I've modified it to remove SNPs on the X chromosome.
# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/fixed_effect/interaction_noReadCov_DF_ANOVA_add1_noX/Alan

# Script to add functional annotations to SNPs using
# sequence ontology from dbSNP

library(tidyverse)
library(data.table)

conv <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/dbSNP/sequence_ontology_conversion_table.txt",
                   sep = "\t", header = T)
dbsnp <- fread("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/dbSNP/dbSNP_153_common.bed.gz",
                    sep = "\t")
annotated <- merge(conv, dbsnp, by.x = "ID", by.y = "V13") %>%
  select(V1, V2, V3, V4, Name)
write.table(annotated, "dbsnp_annotated.txt",sep = "\t", row.names = F, col.names = F, quote = F)
aux <- mutate(annotated, ID = paste0(V1, ":", V2, "-", V3)) %>% select(ID, V4, Name)



# All tested SNPs
all <- read.table("/wsu/home/groups/piquelab/Alan/splicing/leafcutter/CM/stranded/meta/ASE_SNPs_input.bed.txt",
                  sep = "\t") %>% mutate(ID = paste0("chr", V1, ":", V2, "-", V3)) %>% select(ID)
all_annotated <- inner_join(all, aux, by = "ID")
write.table(all_annotated, "input_SNPs_dbSNP_annotated.txt", sep = "\t", row.names = F, col.names = F, quote = F)


# ASE
ase <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/fixed_effect/interaction_noReadCov_DF_ANOVA_add1_noX/ASE_SNPs.txt") %>%
  setNames("V4")
ase_annotated <- inner_join(ase, all_annotated, by = "V4")
write.table(ase_annotated, "ASE_SNPs_dbSNP_annotated.txt", sep = "\t", row.names = F, col.names = F, quote = F)
ase_notmatch_annotated <- anti_join(all_annotated, ase_annotated, by = "V4")
write.table(ase_notmatch_annotated, "nomatch_ASE_SNPs_dbSNP_annotated.txt", sep = "\t", row.names = F, col.names = F, quote = F)

# cASE
case <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Quasar/cASE_Linear_Model/models/fixed_effect/interaction_noReadCov_DF_ANOVA_add1_noX/Tr_cASE_SNPs.txt") %>%
  setNames("V4")
case_annotated <- inner_join(case, all_annotated , by = "V4")
write.table(case_annotated, "cASE_SNPs_dbSNP_annotated.txt", sep = "\t", row.names = F, col.names = F, quote = F)
case_notmatch_annotated <- anti_join(all_annotated, case_annotated, by = "V4")
write.table(case_notmatch_annotated, "nomatch_cASE_SNPs_dbSNP_annotated.txt", sep = "\t", row.names = F, col.names = F, quote = F)


# Create summary files
system("cat input_SNPs_dbSNP_annotated.txt | cut -f3 | sort | uniq -c > input_SNPs_dbSNP_annotated_summary.txt")
system("cat ASE_SNPs_dbSNP_annotated.txt | cut -f3 | sort | uniq -c > ASE_SNPs_dbSNP_annotated_summary.txt")
system("cat nomatch_ASE_SNPs_dbSNP_annotated.txt | cut -f3 | sort | uniq -c > nomatch_ASE_SNPs_dbSNP_annotated_summary.txt")
system("cat cASE_SNPs_dbSNP_annotated.txt | cut -f3 | sort | uniq -c > cASE_SNPs_dbSNP_annotated_summary.txt")
system("cat nomatch_cASE_SNPs_dbSNP_annotated.txt | cut -f3 | sort | uniq -c > nomatch_cASE_SNPs_dbSNP_annotated_summary.txt")


all_list <- read.table("input_SNPs_dbSNP_annotated_summary.txt", quote = "\"", comment.char = "", stringsAsFactors = F)$V2



# Build color list
colors <- data.frame("Event" = as.data.frame(all_list),
                     "col" = as.data.frame(scales::hue_pal()(18))) %>%
  setNames(c("Event", "col"))






# ASE test
nomatch <- read.table("nomatch_ASE_SNPs_dbSNP_annotated_summary.txt", quote = "\"", comment.char = "", stringsAsFactors = F) %>%
  setNames(c("no_ase", "event"))
match <- read.table("ASE_SNPs_dbSNP_annotated_summary.txt", quote = "\"", comment.char = "", stringsAsFactors = F) %>%
  setNames(c("yes_ase", "event"))
together <- full_join(match, nomatch, by = "event") %>% select(event, yes_ase, no_ase)
together[is.na(together)] <- 0

df <- data.frame()
for (i in all_list){
  temp <- filter(together, event == i)
  df2 <- data.frame("event" = c(temp$yes_ase, temp$no_ase),
                    "not_event" = c(nrow(ase_annotated) - temp$yes_ase, nrow(ase_notmatch_annotated) - temp$no_ase))
  test <- fisher.test(df2, alternative = "two.sided")
  df2 <- data.frame("Event" = i,
                    "OR" = test$estimate[[1]],
                    "Bottom" = test$conf.int[1],
                    "Upper" = test$conf.int[2],
                    "p" = test$p.value)

  df <- rbind(df, df2)
}

df <- mutate(df, padj = p.adjust(p, method = "BH"))

ggplot(data = df, aes(x = Event, OR, y = log10(OR+0.001), ymin = log10(Bottom+0.001), ymax = log10(Upper+0.001))) +
  geom_text(aes(label = ifelse(as.numeric(padj)<0.1, "*", ""), y = -3.25), size = 14, nudge_x = -0.3) +
  geom_pointrange(aes(col = Event)) +
  xlab("") +
  ylab("log10 Odds Ratio") +
  geom_errorbar(aes(ymin = log10(Bottom+0.001), ymax = log10(Upper+0.001), col = Event), width = 0.5, cex = 1) +
  theme(strip.text.y = element_text(hjust = 0, vjust = 1, angle = 180, face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_bw() +
  coord_flip()

ggsave("forestplot_dbsnp_ASEvsINPUT.pdf",
       plot = last_plot(),
       device = "pdf",
       width = 10,
       height = 7,
       dpi = 300)




df2 <- filter(df, padj < 0.1) %>% merge(colors, by = "Event")
ggplot(data = df2, aes(x = Event, y = log10(OR+0.001), ymin = log10(Bottom+0.001), ymax = log10(Upper+0.001))) +
  geom_pointrange(aes(col = Event)) +
  xlab("") +
  ylab("log10 Odds Ratio") +
  geom_errorbar(aes(ymin = log10(Bottom+0.001), ymax = log10(Upper+0.001), col = Event), width = 0.5, cex = 1) +
  theme(strip.text.y = element_text(hjust = 0, vjust = 1, angle = 180, face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_bw() +
  coord_flip()

ggsave("forestplot_dbsnp_ASEvsINPUT_sig_only.pdf",
       plot = last_plot(),
       device = "pdf",
       width = 10,
       height = 7,
       dpi = 300)









# cASE test
nomatch <- read.table("nomatch_cASE_SNPs_dbSNP_annotated_summary.txt", quote = "\"", comment.char = "", stringsAsFactors = F) %>%
  setNames(c("no_cASE", "event"))
match <- read.table("cASE_SNPs_dbSNP_annotated_summary.txt", quote = "\"", comment.char = "", stringsAsFactors = F) %>%
  setNames(c("yes_cASE", "event"))
together <- full_join(match, nomatch, by = "event") %>% select(event, yes_cASE, no_cASE)
together[is.na(together)] <- 0

df <- data.frame()
for (i in all_list){
  temp <- filter(together, event == i)
  df2 <- data.frame("event" = c(temp$yes_cASE, temp$no_cASE),
                    "not_event" = c(nrow(case_annotated) - temp$yes_cASE, nrow(case_notmatch_annotated) - temp$no_cASE))
  test <- fisher.test(df2, alternative = "two.sided")
  df2 <- data.frame("Event" = i,
                    "OR" = test$estimate[[1]],
                    "Bottom" = test$conf.int[1],
                    "Upper" = test$conf.int[2],
                    "p" = test$p.value)

  df <- rbind(df, df2)
}

df <- mutate(df, padj = p.adjust(p, method = "BH"))

ggplot(data = df, aes(x = Event, OR, y = log10(OR+0.001), ymin = log10(Bottom+0.001), ymax = log10(Upper+0.001))) +
  geom_text(aes(label = ifelse(as.numeric(padj)<0.1, "*", ""), y = -3.25), size = 14, nudge_x = -0.3) +
  geom_pointrange(aes(col = Event)) +
  xlab("") +
  ylab("log10 Odds Ratio") +
  geom_errorbar(aes(ymin = log10(Bottom+0.001), ymax = log10(Upper+0.001), col = Event), width = 0.5, cex = 1) +
  theme(strip.text.y = element_text(hjust = 0, vjust = 1, angle = 180, face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_bw() +
  coord_flip()

ggsave("forestplot_dbsnp_cASEvsINPUT.pdf",
       plot = last_plot(),
       device = "pdf",
       width = 10,
       height = 7,
       dpi = 300)




df2 <- filter(df, padj < 0.1) %>% merge(colors, by = "Event")
ggplot(data = df2, aes(x = Event, y = log10(OR+0.001), ymin = log10(Bottom+0.001), ymax = log10(Upper+0.001))) +
  geom_pointrange(aes(col = Event)) +
  xlab("") +
  ylab("log10 Odds Ratio") +
  geom_errorbar(aes(ymin = log10(Bottom+0.001), ymax = log10(Upper+0.001), col = Event), width = 0.5, cex = 1) +
  theme(strip.text.y = element_text(hjust = 0, vjust = 1, angle = 180, face = "bold")) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  theme_bw() +
  coord_flip()

ggsave("forestplot_dbsnp_cASEvsINPUT_sig_only.pdf",
       plot = last_plot(),
       device = "pdf",
       width = 10,
       height = 7,
       dpi = 300)
