# Anthony Findley
# 4/19/19

# Folder: /nfs/rprdata/Anthony/GxE_iPSC/Deep/LCL_combined/P1_R2
# This runs DESeq for each treatment against it's appropriate control in LCLs.

require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(qvalue)
require(plyr)
library(reshape)
require(parallel)
require(BiocParallel)
source('/nfs/rprscratch/Allison/FUNGEI1again/GxE_pipeline/misc/getArgs.R')

anno <- read.table("/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.tsv.gz", header=T, sep="\t", stringsAsFactors=T)

## Get command-line arguments.
defaultList = list(
  cores=14,
bedTranscriptome="/nfs/rprdata/Anthony/data/hg38_annotations/hg38_transcriptome.anno.tsv.gz",
# gcContentFile="/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.2014.faCount.gz",
platePrefix="DLCL_P1R2_04.19.19"
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

## Gene counts: this is our data for anlaysis
readCounts1 <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL1R1/P1_R2/fpkm/GeneCounts/DLCL1R1.data.Rd"
readCounts2 <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL1R2/P1_R2/fpkm/GeneCounts/DLCL1R2.data.Rd"
readCounts3 <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL2R1/P1_R2/fpkm/GeneCounts/DLCL2R1.data.Rd"
readCounts4 <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/DLCL2R2/P1_R2/fpkm/GeneCounts/DLCL2R2.data.Rd"

##anno is just the annotation for all of the transcripts...should be the same for all
load(readCounts1)
load(readCounts2)
load(readCounts3)
load(readCounts4)

data <- cbind(DLCL1R1counts, DLCL1R2counts, DLCL2R1counts, DLCL2R2counts)

## covariates file with experimental information on the samples
cov.file1 <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL1R1_covar.txt"
cov.file2 <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL1R2_covar.txt"
cov.file3 <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL2R1_covar.txt"
cov.file4 <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL2R2_covar.txt"

cv1 <- read.table(cov.file1 , stringsAsFactors=F, header=T, comment="")
cv2 <- read.table(cov.file2 , stringsAsFactors=F, header=T, comment="")
cv3 <- read.table(cov.file3 , stringsAsFactors=F, header=T, comment="")
cv4 <- read.table(cov.file4 , stringsAsFactors=F, header=T, comment="")

cv <- rbind(cv1, cv2, cv3, cv4)

data <- data[,colnames(data) %in% cv$Filename]

data <- data[,-which(colSums(data) < 3000000)] # Gets rid of control DLCL1R2_HT92 with few reads

cv <- cv[cv$Filename %in% colnames(data),]
cv <- cv[order(cv$Filename),]
data <- data[,order(colnames(data))]
which(!colnames(data)==cv$Filename)


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
    design = ~ Ind.Plate + Treatment.ID)
keep <- rowSums(counts(ddsFull)) > 0
ddsFull <- ddsFull[keep,]
# colnames(ddsFull) <- cv$Filename

## Fit the model on the whole plate
system.time(ddsFull <- DESeq(ddsFull,parallel=TRUE))

##That took so fucking long, I don't want to risk losing it
save(list="ddsFull", file=paste(dataDir, '/DESeq2_', platePrefix, '.RData', sep=''))

#####
# Initialize a list so I can save the number of DE genes on the wwwShare folder

DEgene_list <- list()

## Contrast each treatment to its respective control
##res <- ParallelSapply(TreatmentOnlyLevels,function(t){
# res <- sapply(TreatmentOnlyLevels,function(t){

for (i in 1:length(TreatmentOnlyLevels)){
  t <- TreatmentOnlyLevels[i]

  ## Select appropiate control for this treatment
  c <- ControlLevels[cv$Control.ID[cv$Treatment.ID==t][1]]
  cat("-----------------------------------------------------------","\n")
  cat("Processing treatment:",t,treatmentKey[t,"Treatment_Name"],"\n")
  cat("Using control:",c,treatmentKey[c,"Treatment_Name"],"\n")

  res <- results(ddsFull, contrast=c("Treatment.ID",t,c),parallel=TRUE)
  head(res)

  ## Use a filter to help with power
  fThresh <- attr(res,"filterThreshold")
  use <- res$baseMean > fThresh

  ## Plot filter
  fname=paste(plotsDir, '/', platePrefix, "_", gsub(' ', '_', t),"_Filter", ".pdf", sep="")
  pdf(fname);
  aux = attr(res,"metadata");#"filterNumRej")
  plot(aux$filterNumRej,type="b",ylab="number of rejections")
  abline(v=aux$theta[which.max(aux$numRej)],lty=3)
  dev.off();

  fname=paste(plotsDir, '/', platePrefix, "_", gsub(' ', '_', t), ".pdf", sep="")
  pdf(fname, height=500, width=500)
  qqplot(-log10(ppoints(length(res$pvalue))),-log10(res$pvalue),pch=20,cex=0.5)
  title(main=paste0(t,", FDR10%=",sum(res$padj<0.1,na.rm=T),", ",treatmentKey[t,"Short_Name"],", ",sum(use),"/",length(use)))
  abline(0, 1, col='red')
  dev.off()

  ## Get the top differentially expressed genes
  fname=paste(statsDir, '/', platePrefix,"_","DEG_stats_", gsub(' ', '_', t), ".txt",sep="")
  sub.table <- data.frame(res@rownames, res$'padj', res$'pvalue', res$'log2FoldChange', stringsAsFactors=FALSE)
  names(sub.table) <- c('t.id', 'padj', 'pval', 'logFC')
  sub.table$ensg <- anno[sub.table$t.id, 'Gene.stable.ID']
  sub.table$g.id <- anno[sub.table$t.id, 'Gene.name']
  sub.table$gene_type <- anno[sub.table$t.id, 'Gene.type']
  sub.table <- sub.table[!is.na(sub.table$padj), ]
  write.table(sub.table, file=fname, quote=FALSE, row.names=FALSE)
  cat("BH diff. expressed trx.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
  cat("BH unique genes diffexpr.", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), length(unique(anno$Gene.name[keep][(res$padj)<tr])) )),"\n")

  ## ASF: I'm adding output to Allison's script
  ## First, add number of DE genes to initialized DEgene_df
  DE_genes <- as.character(unique(as.character(sub.table[which(abs(sub.table$logFC) > 0.25 & sub.table$padj < 0.1),"g.id"])))
  DEgene_append <- data.frame(t(matrix(c(t, treatmentKey[t,"Treatment_Name"], length(DE_genes)))), stringsAsFactors=F)
  colnames(DEgene_append) <- c("Treatment.ID", "Treatment.Name", "DE.Genes")
  DEgene_list[[i]] <- DEgene_append

  ## Now make QQ-plot
  sub.table <- sub.table[order(sub.table$pval),]
  sub.table$exp <- ppoints(length(sub.table$pval))

  system(paste0("mkdir -p /nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/QQ_plots/P1_R2/", platePrefix))

  pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/QQ_plots/P1_R2/", platePrefix, "/", t, ".qq.pdf"))
  p <- ggplot(sub.table, aes(-log10(exp), -log10(pval)))
  print(p + geom_point() + labs(title = paste0("QQ-plot of DESeq2 p-values for ", platePrefix, ": ", t, " vs ", c), x = "-log10(exp.pvalue)", y = "-log10(obs.pvalue)") +geom_abline(intercept=0, slope=1))
  dev.off()

  }

  ## Return table with p-values, fold diff, qvalue and so on.
  #res
#}); names(res) <- TreatmentOnlyLevels

##Overwrite previous file
##Now include treatment vs control contrast (ACROSS INDIVIDUALS)
save(list=c("res", "ddsFull"), file=paste(dataDir, '/DESeq2_', platePrefix, '.RData', sep=''))

########### Do control vs control

t="CO1"
c="CO2"
con <- results(ddsFull, contrast=c("Treatment.ID","CO1","CO2"),parallel=TRUE)

  fname=paste(statsDir, '/', platePrefix,"_","DEG_stats_", gsub(' ', '_', t), ".txt",sep="")
  sub.table <- data.frame(con@rownames, con$'padj', con$'pvalue', con$'log2FoldChange', stringsAsFactors=FALSE)
  names(sub.table) <- c('t.id', 'padj', 'pval', 'logFC')
  sub.table$ensg <- anno[sub.table$t.id, 'Gene.stable.ID']
  sub.table$g.id <- anno[sub.table$t.id, 'Gene.name']
  sub.table$gene_type <- anno[sub.table$t.id, 'Gene.type']
  sub.table <- sub.table[!is.na(sub.table$padj), ]
  write.table(sub.table, file=fname, quote=FALSE, row.names=FALSE)
  cat("BH diff. expressed trx.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(con$padj)<tr))),"\n")
  cat("BH unique genes diffexpr.", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), length(unique(anno$Gene.name[keep][(con$padj)<tr])) )),"\n")

  ## ASF: I'm adding output to Allison's script
  ## First, add number of DE genes to initialized DEgene_df
  DE_genes <- as.character(unique(as.character(sub.table[which(abs(sub.table$logFC) > 0.25 & sub.table$padj < 0.1),"g.id"])))
  DEgene_append <- data.frame(t(matrix(c(t, treatmentKey[t,"Treatment_Name"], length(DE_genes)))), stringsAsFactors=F)
  colnames(DEgene_append) <- c("Treatment.ID", "Treatment.Name", "DE.Genes")
  DEgene_list[[length(TreatmentOnlyLevels) + 1]] <- DEgene_append

  ## Now make QQ-plot
  sub.table <- sub.table[order(sub.table$pval),]
  sub.table$exp <- ppoints(length(sub.table$pval))

  system(paste0("mkdir -p /nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/QQ_plots/P1_R2/", platePrefix))

  pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/QQ_plots/P1_R2/", platePrefix, "/", t, ".qq.pdf"))
  p <- ggplot(sub.table, aes(-log10(exp), -log10(pval)))
  print(p + geom_point() + labs(title = paste0("QQ-plot of DESeq2 p-values for ", platePrefix, ": ", t, " vs ", c), x = "-log10(exp.pvalue)", y = "-log10(obs.pvalue)") +geom_abline(intercept=0, slope=1))
  dev.off()


####

DEgene_table <- do.call(rbind, DEgene_list)

system("mkdir -p /nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/DE_Genes/P1_R2/")
write.table(DEgene_table, file=paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/DESeq2/DE_Genes/P1_R2/", platePrefix, "_DEgenes.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE)
