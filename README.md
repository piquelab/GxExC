# GxExC
This repository contains scripts used to analyze data from exposing LCLs, IPSCs, and CMs to 28 treatments, as detailed in Findley et al (2021). It is divided into two folders: Shallow and Deep. The Shallow folder contains scripts necessary for aligning the shallow RNA-sequencing reads and running DESeq2 to identify differentially expressed genes. The Deep folder contains scripts to align deep RNA-sequencing reads, identify differentially expressed and differentially spliced genes, identify instances of allele-specific expression (ASE) and conditional ASE (cASE), and partition the variance in gene expression, splicing, and ASE. 

## Shallow

1. align.sh: Aligns fastq's to the human genome (build GRCh37) using HISAT2 and performs QC and deduplication
2. coverageBed_2.25.0.sh: Counts number of aligned reads per transcript
3. bed2GeneCounts.R: Generates gene expression count matrix, where rows are transcripts and columns are sequencing libraries, to be used in DESeq2
4. DESeq.R: Run DESeq2 on gene expression count matrix to identify differentially expressed transcripts

## Deep

### Alignment
1. align.sh: Aligns fastq's to the human genome (build GRCh37) using HISAT2 and performs QC and deduplication
    - merge_bam.sh: Combines bam files from 2 rounds of deep sequencing for CM plates
2. coverageBed_2.25.0.sh: Counts number of aligned reads per transcript
3. bed2GeneCounts.R: Generates gene expression count matrix, where rows are transcripts and columns are sequencing libraries, to be used in DESeq2

### Differential gene expression and splicing
4. DESeq.R: Run DESeq2 on gene expression count matrix to identify differentially expressed transcripts in response to each treatment
    - GO_DEGs.R: Gene ontology analysis comparing differentially expressed genes in each condition to the background of all expressed genes
5. DESeq_cell_treat_interact.R: Run DESeq2 on gene expression count matrix to identify treatment x cell type interactions using a likelihood ratio test.
6. DEG_DSG_enrich.R: Calculate enrichment of differentially expressed genes in differentially spliced genes

### Variance partitioning of gene expression and splicing
7. VarPart_cellTogether_GE.R: Variance partitioning of gene expression on all cell types together
8. VarPart_cellSep_GE.R: Variance partitioning of gene expression on cell types separately
9. VarPart_Splice.R: Variance partitioning of splicing, including both all cell types together and each cell type separately

### ASE analysis
10. ai_processing.R : Creates pileup files, which describes the number of reads mapping to each allele at heterozygous sites
    - Pileup_makefile: Used for submitting each sequencing library to ai_processing.R
11. QuASAR_prep.R: Prepares makefiles to be analyzed by QuASAR
12. combine_controls.R: Combines technical replicates of the 2 controls within each plate
13. QuASAR_pipeline.R: Run QuASAR for each individual-plate combination
14. ASE_barplot.R: Create barplot of ASE per treatment in Figure 3A
15. ASE_cor.R: Calculate correlations in ASE across individual, treatment, etc. and make Figure 3B

### cASE analysis
16. cASE_bigAnalysis_allCells.R: Analysis of cASE and ASE variance on all cell types analyzed together
16. cASE_bigAnalysis_CellSep.R: Analysis of cASE and ASE variance on cell types analyzed separately
17. ASE_cASE_annot_enrich.R: Calculate the enrichment of ASE and cASE SNPs in genomic annotations (Figure 3G)
