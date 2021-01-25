# This is an example of how to run DESeq2 on the shallow sequencing from the LCLs.

require(ggplot2) ## Other packages need to overwrite certain 1.0.1.993 functions
library(DESeq2)
library(qvalue)
require(plyr)
library(reshape)
require(parallel)
require(BiocParallel)
source('/nfs/rprscratch/Allison/FUNGEI1again/GxE_pipeline/misc/getArgs.R')

load("~/rprdata/Allison/FUNGEI_hisat/DESeq2/GeneCounts/FUNGEI1.data.Rd") #Just to get "anno" file to be used later on
rm(data) #Get rid of Allison's "data" file

## Get command-line arguments.
defaultList = list(
  cores=28,
bedTranscriptome="/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.2014.bed.gz",
gcContentFile="/wsu/home/groups/piquelab/data/RefTranscriptome/ensGene.hg19.2014.faCount.gz",
platePrefix="LCL_all_DESeq2_1.22.1"
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
readCounts1 <- "/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/LCL1R1/fpkm/GeneCounts/LCL1R1.data.Rd"
readCounts2 <- "/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/LCL1R2/fpkm/GeneCounts/LCL1R2.data.Rd"
readCounts3 <- "/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/LCL2R1/fpkm/GeneCounts/LCL2R1.data.Rd"
readCounts4 <- "/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/LCL2R2/fpkm/GeneCounts/LCL2R2.data.Rd"


##anno is just the annotation for all of the transcripts...should be the same for all
load(readCounts1)
load(readCounts2)
load(readCounts3)
load(readCounts4)

data <- cbind(LCL1R1counts, LCL1R2counts, LCL2R1counts, LCL2R2counts )

## covariates file with experimental information on the samples
cov.file1 <- "/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/Docs/LCL1R1_covar.txt"
cov.file2 <- "/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/Docs/LCL1R2_covar.txt"
cov.file3 <- "/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/Docs/LCL2R1_covar.txt"
cov.file4 <- "/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/Docs/LCL2R2_covar.txt"
cv1 <- read.table(cov.file1 , as.is=T, header=T, comment="")
cv2 <- read.table(cov.file2 , as.is=T, header=T, comment="")
cv3 <- read.table(cov.file3 , as.is=T, header=T, comment="")
cv4 <- read.table(cov.file4 , as.is=T, header=T, comment="")

cv <- rbind(cv1,cv2,cv3,cv4)



data <- data[,colnames(data) %in% cv$Filename]

n.barcodes <- dim(data)[2]

## master list of treatment names and IDs
treatmentKey <- read.table(paste('~/piquelab/gmb/GxE_full/expression/treatmentKey.txt', sep=''), sep='\t',as.is=TRUE, header=TRUE, comment.char="")

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
rownames(anno) <- anno$t.id

## Annotating transcript length in bp and GC content.
anno2 <- read.table(gcContentFile,as.is=T,sep="\t",header=T,comment="")
rownames(anno2) <- gsub("hg19_ensGene_", "", anno2$X.seq)
anno2$avg.cg <- (anno2$C+anno2$G)/anno2$len
#coding length
anno$codLen <- anno2[anno$t.id,"len"]
#transcript length
anno$txLen <- (anno$c.start-anno$start)+(anno$stop-anno$c.stop)+anno$codLen
anno2$avg.cg <- (anno2$C+anno2$G)/anno2$len
#Avg. CG on the coding part.
anno$avg.cg <- anno2[anno$t.id,"avg.cg"]
rm(anno2)

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


ddsFull <- DESeqDataSetFromMatrix(
    countData = round(data),
    colData = cv,
    design = ~ Ind.Plate + Treatment.ID)
keep <- rowSums(counts(ddsFull)) > 0
ddsFull <- ddsFull[keep,]
colnames(ddsFull) <- cv$Filename


## Fit the model on the whole plate
system.time(ddsFull <- DESeq(ddsFull,parallel=TRUE))

##That took so fucking long, I don't want to risk losing it
save(list="ddsFull", file=paste(dataDir, '/DESeq2_', platePrefix, '.RData', sep=''))

## Contrast each treatment to its respective control
##res <- ParallelSapply(TreatmentOnlyLevels,function(t){
res <- sapply(TreatmentOnlyLevels,function(t){

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
  sub.table$ensg <- anno[sub.table$t.id, 'ensg']
  sub.table$g.id <- anno[sub.table$t.id, 'g.id']
  sub.table <- sub.table[!is.na(sub.table$padj), ]
  write.table(sub.table, file=fname, quote=FALSE, row.names=FALSE)
  cat("BH diff. expressed trx.  ", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), sum(na.omit(res$padj)<tr))),"\n")
  #cat("BH unique genes diffexpr.", sapply(c(0.01,0.05,0.1,0.2),function (tr) c(paste(",",tr*100,"%FDR->"), length(unique(anno$g.id[keep][(res$padj)<tr])) )),"\n")

  ## Return table with p-values, fold diff, qvalue and so on.
  res
}); names(res) <- TreatmentOnlyLevels

##Overwrite previous file
##Now include treatment vs control contrast (ACROSS INDIVIDUALS)
save(list=c("res", "ddsFull"), file=paste(dataDir, '/DESeq2_', platePrefix, '.RData', sep=''))

