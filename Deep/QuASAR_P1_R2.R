library('QuASAR')
library('ggplot2')
library('parallel')


## function UnionExtractFields
## First, get the unique list of SNP positions being interrogated for
## this individual. Then,
UnionExtractFields <- function (fileList, combine = FALSE)
  {
    tmpFile <- scan(pipe("mktemp -t"), character(0))
      system(
        paste("zcat ",
              paste(fileList, collapse = " "),
              " | grep -v -w '^chr\\|^NA' | cut -f 1-7 | sort -k 1,1 -k2,2n | uniq > ",
              tmpFile ))

    anno <- read.table(tmpFile, sep = "\t", as.is = T)
    aux <- ParallelSapply(fileList, function(fn) {
      cat("Processing:", fn, "\n")
      command = paste("module load bedtools; intersectBed -a ", tmpFile, " -b ",
        fn, " -wao | cut -f 1-3,15-17 ", sep = "")
##################################Probably is 'cut 1-4,14-16'!!!! then I can activate again stopifnot(identical()) but I have to check that!
      cat(command,"\n")
      aa <- read.table(pipe(command), sep = "\t", as.is = T,
                       na.strings = ".")
      aa[is.na(aa)] <- 0
      stopifnot(identical(aa[, 1:3], anno[, 1:3]))
      aa[, -(1:3)]
    })

    colnames(anno) = c("chr", "pos0", "pos", "ref", "alt", "rsID", "af")
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

cores = 8 # default
cargs <- commandArgs(trail=TRUE);
##covariateFile <- '../../0-Treatment_sequencing/FUNGEI_RNAseq_20161212.csv'
#covariateFile <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/DLCL2R2_covar.txt"

# Plate="DLCL1R2"
# Individual="GM18855"

homeDir <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/"

if(length(cargs)>=1)
  #!#covariateFile <- cargs[1]
  Individual <- cargs[ 1 ]
if(length(cargs)>=2)
  #!#covariateFile <- cargs[1]
  Plate <- cargs[ 2 ]
if(length(cargs)>=3)
  cores <- cargs[3]


## Helper functions
ParallelSapply <- function(...,mc.cores=cores){
  simplify2array(mclapply(...,mc.cores=mc.cores))
}



covariateFile <- paste0("/nfs/rprdata/Anthony/GxE_iPSC/Deep/Docs/",Plate,"_covar.txt")


##################################################################
## specify covariate file and pileup directories
##################################################################
##cmd <- paste0('cat ', covDir, 'GxE_DP*_covariates.txt | grep -w ', cell.type)
#cv <- read.table(covariateFile, as.is=TRUE, sep='\t', header=TRUE, row.name=1)
cv <- read.table(covariateFile, stringsAsFactors=FALSE, sep='\t', header=TRUE)#!#

# ASF 8/31/18: I included the next two commands because the script fails if everything in the covariate file doesn't have a pileup. So this removes those that don't have a pileup file

existingBarcodes <- system(paste0("ls ./combined_pileups/", Plate, "/*.pileup.clean.bed.gz | cut -d/ -f4 | cut -d. -f1"), intern=T)
cv <- cv[which(cv$Filename %in% existingBarcodes), ]

##cv$Treatment <- gsub(' ',  '_', cv$Treatment)

## Ideally, this could be implemented several ways:
## a) one covariate file per cell type
## b) one master covariate file, cell type as input & filter the cv object
## This pipeline currently expects (a)

cell.type <- unique(cv$CellType)

## Loop over individuals
## TODO: ultimately, should rename CellLine -> Sample in the cov file
#!#ParallelSapply(unique(cv$CellLine), function(sample) {
#!#for( sample in unique( cv$CellLine ) ) {
#!#for( sample in unique( cv$dbGaP.id ) ){






ind <- cv$Individual==Individual
fileList <- paste0("./combined_pileups/", Plate, "/", cv$Filename[ind],".pileup.clean.bed.gz")
samples <- cv$Filename[ind]
cov.file <- cv[ind, ]


ase.dat <- UnionExtractFields(fileList, combine=TRUE)

colnames(ase.dat$ref)=samples
colnames(ase.dat$alt)=samples
colnames(ase.dat$err)=samples



outFolder = "./QuASAR_input/"
system(paste("mkdir -p ",outFolder))

save(list=c("ase.dat","cov.file"),file=paste0(outFolder,Plate,"-",Individual,"_quasarIn.Rd"))


## the end ##
