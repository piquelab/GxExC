##################################################################
## Prepare data for input to QuASAR for joint genotyping and ASE
##
## Arguments: cell.type
## Return Values:
##################################################################

library('QuASAR')
library('ggplot2')
library('parallel')


folder <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/"
platePrefix <- "DCM1R1"
cores = as.integer(Sys.getenv("NCPUS",1)) # default


##################################################################
## specify covariate file and pileup directories
##################################################################
#covDir <- '/wsu/home/fx/fx78/fx7820/rprdata/Anthony/GxE_iPSC/Shallow/Docs'

cov.path <- paste0(folder, "Docs/", platePrefix, "_covar.txt")
cv <- read.table(cov.path, header=TRUE, as.is=TRUE)


## Finding the pileup files

pileup.folder <- paste0(folder, platePrefix, "/R2_R3/Pileups/")
fileList <- dir(pileup.folder,"pileup.clean.bed.gz$")
fNames <- gsub(".pileup.clean.bed.gz","",fileList)
fileList <- paste0(pileup.folder,fileList)
names(fileList) <- fNames

## Helper functions
ParallelSapply <- function(...,mc.cores=cores){
  ## mc.preschedule=F because it doesn't work otherwise
  simplify2array(mclapply(...,mc.cores=mc.cores, mc.preschedule=F))
}

## function UnionExtractFields
## First, get the unique list of SNP positions being interrogated for
## this individual. Then,
UnionExtractFields <- function (fileList, combine = FALSE)
  {
    tmpFile <- scan(pipe("mktemp -t"), character(0))
    system(
      paste("zcat ",
            paste(fileList, collapse = " "),
#            " | grep -v -w '^chr\\|^NA' | cut -f 1-4,6-7 | sortBed -i stdin | uniq | gzip > ",
            " | grep -v -w '^chr\\|^NA' | cut -f 1-4,6-7 | sortBed -i stdin | uniq > ",
            tmpFile))
    #anno <- read.table(gzfile(tmpFile), sep = "\t", as.is = T)
    anno <- read.table(tmpFile, sep = "\t", as.is = T)
    ##
    aux <- ParallelSapply(fileList, function(fn) {
    ##aux <- sapply(fileList, function(fn) {
      cat("Processing:", fn, "\n")
      command = paste("intersectBed -a ", tmpFile, " -b ",
        fn, " -wao | cut -f 1-3,14-16 ", sep = "")
      aa <- read.table(pipe(command), sep = "\t", as.is = T,
                       na.strings = ".")
      aa[is.na(aa)] <- 0
      stopifnot(identical(aa[, 1:3], anno[, 1:3]))
      aa[, -(1:3)]
    })
    colnames(anno) = c("chr", "pos0", "pos", "ref", "rsID", "af")
    Ref <- as.matrix(do.call(cbind, aux[1, ]))
    Alt <- as.matrix(do.call(cbind, aux[2, ]))
    Err <- as.matrix(do.call(cbind, aux[3, ]))
    return.list <- list(ref = Ref, alt = Alt, err = Err, anno = anno)
    if (combine == TRUE) {
      allRef <- apply(Ref, MARGIN = 1, sum)
      allAlt <- apply(Alt, MARGIN = 1, sum)
      allErr <- apply(Err, MARGIN = 1, sum)
      return.list$all <- as.matrix(cbind(allRef, allAlt, allErr))
    }
    return(return.list)
  }


##################################################################
## prepare the relevant samples
##################################################################

ase.dat <- UnionExtractFields(fileList, combine=TRUE)

colnames(ase.dat$ref) <- names(fileList)
colnames(ase.dat$alt) <- names(fileList)
colnames(ase.dat$err) <- names(fileList)

rownames(cv) <- cv$sample.ID

cv <- cv[names(fileList),]

save(list=c('ase.dat', 'cv'), file=paste0(folder, platePrefix, "/R2_R3/GenotypeQC/",platePrefix,"_quasarIn.Rd"))

####################################################################


## the end ##

(base) [fx7820@warrior GenotypeQC]$ cat QuASAR_prep.R ^C
(base) [fx7820@warrior GenotypeQC]$ ls
DCM1R1_quasarIn.Rd  Genotype.QC.R  QuASAR_prep.R
(base) [fx7820@warrior GenotypeQC]$ cat Genotype.QC.R
library(data.table)

library(parallel)


folder <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/"
platePrefix <- "DCM1R1"
cores = as.integer(Sys.getenv("NCPUS",1)) # default

cargs <- commandArgs(trail=TRUE);
if(length(cargs)>=1)
  platePrefix <- cargs[1]
if(length(cargs)>=2)
  cores <- cargs[2]


load(file=paste0(folder, platePrefix, "/R2_R3/GenotypeQC/",platePrefix,"_quasarIn.Rd"))

tot <- (ase.dat$ref + ase.dat$alt)
rho.mat <- (ase.dat$ref+1)  / (ase.dat$ref + ase.dat$alt+2)
rownames(rho.mat) <- ase.dat$anno$rsID
print( head( ase.dat$anno$rsID ) )

dat <- rho.mat[rowSums(tot>15)>round(ncol(tot)*0.5), ]
snps <- length(rownames(dat))


##Correlation matrix
dat.cor <- cor(dat)

cn <- colnames(dat.cor)

library(pheatmap)

##rownames( cv ) <-   cv$sample.ID

cov.path <- paste0(folder, "Docs/", platePrefix, "_covar.txt")
cv <- read.table(cov.path, header=TRUE, as.is=TRUE)
rownames(cv) <- cv$Filename


col_anno <- data.frame(row.names=colnames(dat.cor),dbgap=cv[cn,"Individual"])

pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/GenotypeQC/R2_R3/", platePrefix, "_pheatmap.pdf"))
pheatmap(dat.cor,annotation_col=col_anno)
dev.off()

pdf(paste0("/nfs/rprdata/Anthony/GxE_iPSC/Deep/plots/GenotypeQC/R2_R3/", platePrefix, "_pheatmap.pdf"))
pheatmap(dat.cor,annotation_col=col_anno)
dev.off()

pc <- princomp(dat.cor)

library(ggplot2)

df <- data.frame(x=pc$loadings[,1],y=pc$loadings[,2],CellLine=cv[cn,"Individual"])

pdf(paste0("/nfs/rprscratch/wwwShare/Anthony/GxE_iPSC/Deep/GenotypeQC/R2_R3/", platePrefix, "_PCscatter.pdf"))
ggplot(df,aes(x=x,y=y,color=CellLine)) + geom_point()
dev.off()

pdf(paste0("/nfs/rprdata/Anthony/GxE_iPSC/Deep/plots/GenotypeQC/R2_R3/", platePrefix, "_PCscatter.pdf"))
ggplot(df,aes(x=x,y=y,color=CellLine)) + geom_point()
dev.off()
