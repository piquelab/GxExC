# Anthony Findley
# 3/30/20

# Test for cell type specific gene expression. Throw out water controls from CM1R2 because there seemed to be a problem with them. 
# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/All/DESeq

library(tidyverse)
library(DESeq2)
library(qvalue)
library(reshape)
require(parallel)
require(BiocParallel)
source('/nfs/rprscratch/Allison/FUNGEI1again/GxE_pipeline/misc/getArgs.R')

anno <- read.table("/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.tsv.gz", header=T, sep="\t", stringsAsFactors=T)

## Get command-line arguments.
defaultList = list(
  cores=28,
bedTranscriptome="/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.tsv.gz",
# gcContentFile="/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.2014.faCount.gz",
platePrefix="All_CellTypeSpecific_03.30.20"
  )
args <- getArgs(defaults=defaultList)
platePrefix      <- args$platePrefix
cores            <- as.numeric(args$cores)
bedTranscriptome <- args$bedTranscriptome
gcContentFile    <- args$gcContentFile

print(args)

timestamp()

ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}

## To run DESeq2 in parallel, using the
## BiocParallel library
register(MulticoreParam(cores))

# Read data in 
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC1R1/P1_R2/fpkm/GeneCounts/DIPSC1R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC2R1/P1_R2/fpkm/GeneCounts/DIPSC2R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC1R2/P1_R2/fpkm/GeneCounts/DIPSC1R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DIPSC2R2/P1_R2/fpkm/GeneCounts/DIPSC2R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM2R2/R2_R3/fpkm/GeneCounts/DCM2R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM1R2/R2_R3/fpkm/GeneCounts/DCM1R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM2R1/R2_R3/fpkm/GeneCounts/DCM2R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DCM1R1/R2_R3/fpkm/GeneCounts/DCM1R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL1R1/P1_R2/fpkm/GeneCounts/DLCL1R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL2R1/P1_R2/fpkm/GeneCounts/DLCL2R1.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL1R2/P1_R2/fpkm/GeneCounts/DLCL1R2.data.Rd")
load("/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL2R2/P1_R2/fpkm/GeneCounts/DLCL2R2.data.Rd")

data <- cbind(DIPSC1R1counts, DIPSC2R1counts, DIPSC2R2counts, DIPSC1R2counts, DCM1R1counts, DCM2R1counts, DCM1R2counts, DCM2R2counts, DLCL1R1counts, DLCL2R1counts, DLCL2R2counts, DLCL1R2counts)

#Add covariate information to dge$samples
cv <- read.table("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/all_covar.txt", stringsAsFactors=F, header=T)

# Remove water control from DCM1R2
cv <- cv[ -which(cv$Plate.ID == "DCM1R2" & cv$Treatment.ID == "CO1"), ]
cv[ which(cv$Plate.ID == "DCM1R2" & cv$Control.ID == "CO1"), "Control.ID"] <- "CO2"

data <- data[,which(colnames(data) %in% cv$Filename)]
cv <- cv[which(cv$Filename %in% colnames(data)),]

# cv <- cv[which(cv$Filename %in% colnames(counts(ddsFull))),] # Use only if you're loading in ddsFull object

n.barcodes <- dim(data)[2]

## master list of treatment names and IDs
treatmentKey <- read.table(paste('/wsu/home/groups/piquelab/gmb/GxE_full/expression/treatmentKey.txt', sep=''), sep='\t',as.is=TRUE, header=TRUE, comment.char="")

cv$Control.ID <- treatmentKey[cv$Treatment.ID,"Control.ID"]
cv$Treatment <- treatmentKey[cv$Treatment.ID,"Common_Name"]

##################################################################
## assign variables, load data, and load experiment information
topDirectory <- 'out_data_'
outDir <- paste(topDirectory, platePrefix, sep='')
system(paste("mkdir -p",outDir))
##
plotsDir <- paste(outDir, '/plots', sep='')
system(paste("mkdir -p", plotsDir))
##
statsDir <- paste(outDir, '/stats', sep='')
system(paste("mkdir -p", statsDir))
##
dataDir <- paste(outDir, '/data_objects', sep='')
system(paste("mkdir -p", dataDir))
##
degDir <- paste(outDir, '/data_DEG', sep='')
system(paste("mkdir -p", degDir))
##
qcDir <- paste(outDir, '/QC', sep='')
system(paste("mkdir -p", qcDir))

############ QC ############
myNormedData <- ParallelSapply(seq_along(1:n.barcodes), FUN=function(ii){
      qqnorm(rank(data[, ii], ties.method = "random"), plot = F)$x
  })
myCor <- cor(myNormedData, method='pearson')
medianCor <- apply(myCor,1,median)
medmedCor <- median(medianCor)

## median absolute deviation
robust.sd <- median(abs(medianCor-medmedCor))*1.4826
thresh2 <- medmedCor-robust.sd*5

thresh <- .2
idxBlackListed <- medianCor<thresh
cat("##Num blacklist:",sum(idxBlackListed),"\n")

hist.name <- paste(qcDir, '/', platePrefix, "_QC_corr_hist", ".pdf", sep="")
plotCor <- data.frame(correlations=c(myCor))
pdf(hist.name)
m <- ggplot(plotCor, aes(x=correlations))
m + geom_histogram() + ggtitle(paste("Library correlations  ", platePrefix, sep='')) + geom_vline(xintercept = thresh2, colour="red", linetype="longdash")
dev.off()

blackBarcodes <- which(idxBlackListed)
if(length(blackBarcodes > 0)){
  fullBarcodes <- data.frame(blackListedBC=paste('HT', blackBarcodes, sep=''))
  bcName <- paste(qcDir, '/', platePrefix, "_QC_blackList", ".txt", sep="")
  ## write.table(fullBarcodes, file=bcName, quote=FALSE, row.names=FALSE, col.names=FALSE)
}

heat.name <- paste(qcDir, '/', platePrefix, "_QC_corr_heat", ".pdf", sep="")
pdf(heat.name)
image(x=1:n.barcodes, y=1:n.barcodes, myCor, axes=FALSE,xlab="",ylab="")
axis(1,at=c(1,(1:12)*8),las=1,cex=0.5,lwd.ticks=2)
axis(1,at=1:n.barcodes,label=NA,las=1,cex=0.5)
axis(2,at=c(1,(1:12)*8),las=1,cex=0.5,lwd.ticks=2)
axis(2,at=1:n.barcodes,label=NA,las=1,cex=0.5)
title(main=paste("Quantile Normalized Correlation ", platePrefix, sep=''))
dev.off()
############ /QC ############

cv$BlackList=FALSE
cv$BlackList[idxBlackListed] <- TRUE

##################################################################
## anno is an object in P*.data.Rd
## Naming the rows with the transcript ID.
rownames(anno) <- anno$Transcript.stable.ID

## Annotating transcript length in bp and GC content.
# anno2 <- read.table(gcContentFile,as.is=T,sep="\t",header=T,comment="")
# rownames(anno2) <- gsub("hg19_ensGene_", "", anno2$X.seq)
# anno2$avg.cg <- (anno2$C+anno2$G)/anno2$len
#coding length
# anno$codLen <- anno2[anno$t.id,"len"]
#transcript length
# anno$txLen <- (anno$c.start-anno$start)+(anno$stop-anno$c.stop)+anno$codLen
# anno2$avg.cg <- (anno2$C+anno2$G)/anno2$len
#Avg. CG on the coding part.
# anno$avg.cg <- anno2[anno$t.id,"avg.cg"]
# rm(anno2)

##################################################################
## Manual conversion to R factor objects:
cv$Treatment.ID <- factor(cv$Treatment.ID)
TreatmentLevels <- levels(cv$Treatment.ID)
cat("#",TreatmentLevels,"\n")
cv$Control.ID <- factor(cv$Control.ID)
ControlLevels <- levels(cv$Control.ID)
cat("#",ControlLevels,"\n")
TreatmentOnlyLevels <- TreatmentLevels[!(TreatmentLevels %in% ControlLevels)]
cat("#",TreatmentOnlyLevels,"\n")
ControlOnlyLevels <- TreatmentLevels[(TreatmentLevels %in% ControlLevels)]
cat("#",ControlOnlyLevels,"\n")

cv$Plate.ID <- factor(cv$Plate.ID)
cv$Individual <- factor(cv$Individual)
Individual.idLevels <- levels(cv$Individual)
cat("#",Individual.idLevels,"\n")
PlateLevels <- levels(cv$Plate.ID)
cat("#",PlateLevels,"\n")
cv$CellType <- factor(cv$CellType)

##Won't run with individual + plate name because some of them are linear combinations of the others, so we'll make a new variable combining the information...it will really be more like individual*plate though...so I dunno. Otherwise I would have to break it up into many variables (kind of like how I did the LRT model for microNEW. Try this way first because it's easier

## Just print some tables to see that cv looks fine:
table(cv$Control.ID,cv$Treatment.ID)
table(cv$Treatment.ID,cv$Individual)

##################################################################
## Preparing data for DEseq:
## Combine processed data into a DESeqDataSet
## & remove genes with very low coverage/expression

allColSamples <- paste(cv$Individual, cv$Treatment.ID, sep=".")
cat("#",allColSamples,"\n")

cv <- (cv[which(cv$Filename %in% colnames(data)),]) #So cv and data match
cv <- cv[order(cv$Filename),]
data <- data[,order(colnames(data))]
colnames(data)==cv$Filename #Check that files are in order

cv$Ind.Plate <- factor(paste0(cv$Individual,".",cv$Plate.ID))

rownames(cv) <- cv$Filename # DESeq now wants rownames of cv to be colnames of data

ddsFull <- DESeqDataSetFromMatrix(
    countData = round(data),
    colData = cv,
    design = ~ CellType + Treatment.ID + CellType:Treatment.ID)
keep <- rowSums(counts(ddsFull)) > 0
ddsFull <- ddsFull[keep,]
# colnames(ddsFull) <- cv$Filename

## Fit the model on the whole plate
system.time(ddsFull <- DESeq(ddsFull, test="LRT", reduced = ~ CellType + Treatment.ID, parallel=TRUE))

## Error?: 1 rows did not converge in beta, labelled in mcols(object)$fullBetaConv. Use larger maxit argument with nbinomLRT


##That took so fucking long, I don't want to risk losing it
save(list="ddsFull", file=paste(dataDir, '/DESeq2_', platePrefix, '.RData', sep=''))

res <- results(ddsFull,parallel=TRUE)

res_noNA <- res[!is.na(res$padj),]
res_noNA <- merge(data.frame(res_noNA), anno[,c("Transcript.stable.ID", "Gene.type", "Gene.stable.ID", "Gene.name")], by.x=0, by.y="Transcript.stable.ID")

len <- length(unique(as.character(res_noNA[which(res_noNA$padj < 0.1 & res_noNA$Gene.type == "protein_coding" & (res_noNA$log2FoldChange > 0.25 | res_noNA$log2FoldChange < -0.25)),"Gene.name"])))

res_noNA <- res_noNA[order(res_noNA$pvalue),]
res_noNA$exp <- ppoints(length(res_noNA$pvalue))

res_noNA %>% filter(Gene.type == "protein_coding") %>% select(Row.names, padj, pvalue, log2FoldChange, Gene.stable.ID, Gene.name) -> res_forWrite
colnames(res_forWrite) <- c("t.id", "padj", "pval", "logFC", "ensg", "g.id")

write.table(res_forWrite, file="/nfs/rprdata/Anthony/GxE_iPSC/Deep/Paper/Supp_Tables/DESeq_CellSpecific.tab", col.names=T, row.names=F, quote=F, sep="\t")
