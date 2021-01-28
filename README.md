# GxExC
This repository contains scripts used to analyze data from exposing LCLs, IPSCs, and CMs to 28 treatments, as detailed in Findley et al (2021). It is divided into two folders: Shallow and Deep. The Shallow folder contains scripts necessary for aligning the shallow RNA-sequencing reads and running DESeq2 to identify differentially expressed genes. The Deep folder contains scripts to align deep RNA-sequencing reads, identify differentially expressed and differentially spliced genes, identify instances of allele-specific expression (ASE) and conditional ASE (cASE), and partition the variance in gene expression, splicing, and ASE. 

# Shallow

1. align.sh: Aligns fastq's to the human genome (build GRCh37) using HISAT2 and performs QC and deduplication
2. coverageBed_2.25.0.sh: Counts number of aligned reads per transcript
3. bed2GeneCounts.R: Generates gene expression count matrix, where rows are transcripts and columns are sequencing libraries, to be used in DESeq2
4. DESeq.R: Run DESeq2 on gene expression count matrix to identify differentially expressed transcripts

# Deep

1. 
