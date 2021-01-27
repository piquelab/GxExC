##################################################################
## Software suite for joint genotyping and ASE inference on multiple
## experiments and samples
##
## Arguments: plate, cell.line, covariate.file, pileup.dir, out.dir
## Return Values: Genotypes, model.convergence, inference, metaData
##################################################################

library('QuASAR')
library('ggplot2')
library('data.table')
library('parallel')
library('qqman')
require('qvalue')

cores = 8 #!# as.numeric(Sys.getenv("NCPUS","1")) # default enviroment or 1

quasarInFolder= "./QuASAR_input/"
homeDir <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/"

# Plate="DLCL2R2"
# Individual="GM19209"

cargs <- commandArgs(trail=TRUE)
if(length(cargs)>=1)
  #!#covariateFile <- cargs[1]
  Individual <- cargs[ 1 ]
if(length(cargs)>=2)
  #!#covariateFile <- cargs[1]
  Plate <- cargs[ 2 ]
if(length(cargs)>=3)
  cores <- cargs[3]

load(paste0(quasarInFolder,Plate,"-",Individual,"_quasarIn.Rd"))

## Helper functions
ParallelSapply <- function(...,mc.cores=cores){
    simplify2array(mclapply(...,mc.cores=mc.cores))
  }


output.folder <- paste0('QuASAR_results_no500Thresh/', Individual, '/')
system(paste0('mkdir -p ', output.folder, '/plots/QQ'))
system(paste0('mkdir -p ', output.folder, '/plots/QC/oraf'))
system(paste0('mkdir -p ', output.folder, '/plots/QC/coverage'))



ase.dat$anno$af <- as.numeric(ase.dat$anno$af)
ase.dat$anno$af[ase.dat$anno$af<0.01] <- 0.01
ase.dat$anno$af[ase.dat$anno$af>0.99] <- 0.99
ase.dat$anno$af[is.na(ase.dat$anno$af)] <- 0.5

#!#Necessary variables to run the combine of the controls
##barcodes <- paste0(cov.file$Plate.ID, "-HT", cov.file$Barcode.ID) #cov.file$Barcode.ID

#!#treatments <- cov.file$Treatment
treatments <- cov.file$Treatment.ID
n.treatments <- length(treatments)

#!#controls <- unique(cov.file$Control.ID)
controls <- c( 'CO1', 'CO2' )

min.cov <- 15
ase.dat.gt <- PrepForGenotyping(ase.dat, min.coverage=min.cov)
str(ase.dat.gt)



  ##################################################################
  ## QC
  ##################################################################
  N <- ncol(ase.dat.gt$ref)
  if(TRUE){
    for(ii in 1:N){
      ##
      oraf <- ase.dat.gt$ref[, ii]/(ase.dat.gt$ref[, ii]+ase.dat.gt$alt[, ii])
      ind <- ((ase.dat.gt$ref[, ii] + ase.dat.gt$alt[, ii]) > min.cov) & (oraf>0) & (oraf<1)
      oraf <- oraf[ind]
      aux <- data.frame(oraf=oraf)
      ##
      pdf.file <- paste0(output.folder, '/plots/QC/oraf/', colnames(ase.dat.gt$ref)[ii],'_QC_oraf.pdf', sep='')
      pdf(file=pdf.file)
      ##hist(oraf, breaks=100)
      qc_plot <- ggplot(aux, aes(x=oraf))
      print(qc_plot +
             geom_histogram(binwidth=0.01) +
               ggtitle(paste(Individual, '_', ii, sep=''))) +
               theme_bw()
      dev.off()
    }}


  ## ##################################################################
  ## ## Collapse the controls
  ## ##################################################################
  ## plates <- unique(gsub('-HT.*', '', barcodes))

  ## controlref <- do.call('cbind', lapply(plates, function(p) {
  ##  do.call('cbind', lapply(controls, function(this){rowSums(ase.dat.gt$ref[, barcodes[ intersect(grep(this, treatment.IDs), grep(p, barcodes))  ] ])}))
  ## }))
  ## colnames(controlref) <- sapply(plates, function(p) { sapply(controls, function(c) { paste0(p, '-', c) }) })

  ## controlalt <- do.call('cbind', lapply(plates, function(p) {
  ##  do.call('cbind', lapply(controls, function(this){rowSums(ase.dat.gt$alt[, barcodes[ intersect(grep(this, treatment.IDs), grep(p, barcodes))  ]  ])}))
  ## }))
  ## colnames(controlalt) <- sapply(plates, function(p) { sapply(controls, function(c) { paste0(p, '-', c) }) })

  ## ## If any are zero, i.e., in one plate & not the other, remove
  ## torm = which(apply(controlref, 2, sum)==0)
  ## if ( length(torm) > 0 ) {
  ##  controlref <- controlref[, -which(apply(controlref, 2, sum) == 0) ]
  ## }
  ## torm = which(apply(controlalt, 2, sum)==0)
  ## if ( length(torm) > 0 ) {
  ##  controlalt <- controlalt[, -which(apply(controlalt, 2, sum) == 0) ]
  ## }

  ## finalref <- cbind(controlref, ase.dat.gt$ref[, barcodes[ -grep("CO", treatment.IDs) ]])
  ## finalalt <- cbind(controlalt, ase.dat.gt$alt[, barcodes[ -grep("CO", treatment.IDs) ]])

  ##################################################################
  ## QuASAR Model Fitting
  ## ase.joint ~ object for joint genotyoping across all samples
  ##################################################################
  ## Ensure the input to fitAseNullMulti is identical to ase.dat.final below

  finalref <- ase.dat.gt$ref
  finalalt <- ase.dat.gt$alt
  ase.joint <- fitAseNullMulti(finalref,finalalt , log.gmat=log(ase.dat.gt$gmat))


  ## Mandatory objects for inference
  ase.joint.eps <- ase.joint$eps
  n.eps <- length(ase.joint.eps)

  #!#sample.names <- colnames(ase.dat.gt$ref)
  sample.names <- colnames( finalref )

##cv2 <- as.data.frame(cv)
##rownames(cv2) <- cv2$Filename
#!#treatments.final <- cv2[cv2$Filename,]$Treatment.ID
##treatments.final <- cv2[cv2$Filename,]$`Plate.ID-Treatment.ID`

##  sample.names <- treatments_collapsed
##  treatments.final <- treatments_collapsed

  ## probably not correct, ase.dat.gt$gmat should be ase.joint$gmat
 ## ase.dat.final <- list(ref=finalref, alt=finalalt, gmat=ase.dat.gt$gmat, annotations=ase.dat.gt$annotations)


ase.dat.final <- list(ref=finalref, alt=finalalt, gmat=ase.joint$gmat, annotations=ase.dat.gt$annotations)

  ##################################################################
  ## Output model data; genotypes, etc.
  ## ase.joint ~ object for joint genotyoping across all samples
  ##################################################################
  out.gts <- data.frame(rsID=ase.dat.final$annotations$rsID, g0=ase.joint$gt[, 'g0'], g1=ase.joint$gt[, 'g1'], g2=ase.joint$gt[, 'g2'])

  ###EXPLANATION OF THE REASON WHY IT HAS BEEN ADDED THE FOLLOWING LINES#######
  #QuASAR seems to be too conservative for the homozygotes SNPs, so it tends to
  #call homozygotes SNPs with allelic imbalance lower than 5 reads (that could
  #be consider heterozygous). Enableing these lines the homozygotes with allelic
  #allelic imbalance lower than 5 reads are deleted from the genotyping file.
  #!#reads_out.gts <- data.frame( out.gts, ref = rowSums( ase.dat.final$ref ), alt = rowSums( ase.dat.final$alt ) )
  #!#out.gts_AS <- reads_out.gts[ reads_out.gts$g1 > 0.9 | ( ( reads_out.gts$g0 > 0.9 | reads_out.gts$g2 > 0.9 ) & abs( reads_out.gts$ref - reads_out.gts$alt ) > 5 ), ]
  #!#out.gts <- out.gts_AS[ , c( 'rsID', 'g0', 'g1', 'g2' ) ]
  ###END OF THE CHANGE#########################################################

  dat.name <- paste(output.folder, "/",Individual, "_", Plate, '_genotypes.txt.gz', sep='')
  write.table(out.gts, file=gzfile(dat.name), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

  #This genotype file contain the allele count, it could be useful, if you want to run RASQUAL, to generate a VCF file starting from this one do the imputation and then run RASQUAL directly (you have to enable the previous 10 lines to make it work)
  #!#dat.name <- paste(output.folder, "/",cell.line,'_genotypesAS.txt.gz', sep='')
  #!#write.table(out.gts_AS, file=gzfile(dat.name), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")

## Here we may need ot merge with 1000G genotypes and decide what we want to use of what we consider real heterozygotes

totref <- rowSums(ase.dat.final$ref)
totalt <- rowSums(ase.dat.final$alt)

hetindic <- ase.joint$gt[,"g1"]>0.99

# I'm commenting out the below threshold because we think it limits the number of SNPs we test.
truehet <- hetindic # & (totref>=5 & totalt>=5 & ( totref + totalt ) <= 500 )

table(truehet,hetindic)


  ##################################################################
  ## inference
  ##################################################################
  inference.data <- lapply(seq_along(1:n.eps), function(ii){
    ##################################################################
    ## sample ~ current sample to assess ASE
    ## this.sample ~ sample name
    ## coverage ~ coverage for this sample
    ## coverage.floor ~ minimum sample wide coverage
    ## coverage.ind ~ indicator for sufficient coverage of this sample
    ## ref ~ reference count for this sample with sufficient covergae
    ## alt ~ alternate count for this sample with sufficient covergae
    ## phi ~ genotype priors for this sample
    ## eps ~ jointly inferred error rate for this sample
    ## het ~ jointly inferred heterozygote probabilities for sites with
    ##       sufficient coverage
    ## het.ind ~ indicator for heterozygotes with p > .99
    ## annotations ~ annotations filtered by coverage and het probability
    ## q.thresh ~ q-value threshold for declaring signifigance
    ## DEBUGGING
    ## ii <- 1

    sample <- ii
    this.sample <- sample.names[sample]
    coverage <- (ase.dat.final$ref[, sample] + ase.dat.final$alt[, sample])
    coverage.floor <- 5
    coverage.ind <- (coverage>coverage.floor)
    ref <- ase.dat.final$ref[coverage.ind, sample]
    alt <- ase.dat.final$alt[coverage.ind, sample]
    tot <- ref + alt
    phi <- ase.dat.final$gmat[coverage.ind]
    eps <- ase.joint.eps[sample]

    #het <- ase.joint$gt[coverage.ind, 2]
##    het <- ase.joint$gt[ coverage.ind, 'imputed' ]

##    het.ind <- (het > 0.99)
    het.ind <- truehet[coverage.ind]
    numb.hets <- sum(het.ind)
    annotations <- ase.dat.final$annotations[coverage.ind, ][het.ind, ]
    q.thresh <- 0.2
##    this.treatment <- treatmentIDs_final[sample]
##    this.treatment <- treatments.final[ii]

    cat("#Hets:", this.sample, numb.hets, "\n")

    res <- fitQuasarMpra(ref[het.ind],alt[het.ind],rep(0.5,sum(het.ind)),eps=eps)

    ## output complete information for every heterozygote
    complete.dat <- annotations
##    complete.dat$cell.type <- cell.type
    complete.dat$cell.line <- Individual
    complete.dat$treatment <- this.sample ##this.treatment
    complete.dat$ref.reads <- ref[het.ind]
    complete.dat$alt.reads <- alt[het.ind]
    complete.dat$beta <- res$betas.beta.binom
    complete.dat$beta.se <- res$betas_se
    complete.dat$pval <- res$pval3
    complete.dat$qval <- res$padj_quasar

    filename.all <- paste0(output.folder, '/', Individual, '_', this.sample, '_allOutput.txt')
    write.table(complete.dat, file=filename.all, row.names=FALSE, col.name=TRUE, quote=FALSE)

    ## return data frame
##    rsID <- annotations
##    betas <- betas.beta.binom
##    dat <- data.frame(annotations$rsID, annotations$chr, annotations$pos0, rho=rho3, betas, betas.se, qv=qvals.qv3, pval=pval3, refCount=ref[het.ind], altCount=alt[het.ind])
    dat <- data.frame(annotations$rsID, annotations$chr, annotations$pos0,
                      rho=plogis(res$betas.beta.binom), betas=res$betas.beta.binom, betas.se=res$betas_se,
                      qv=res$padj_quasar,pval=res$pval3, refCount=ref[het.ind], altCount=alt[het.ind])

    meta.dat <- data.frame(hets=numb.hets, pi0=NA,
                           qv.05=sum(dat$qv<0.05),
                           qv.1=sum(dat$qv<0.1),
                           qv.2=sum(dat$qv<0.2),
                           qv.5=sum(dat$qv<0.5),
                           mean.rho=mean(dat$rho),
                           median.rho=median(dat$rho))
    ##print(meta.dat)
    temp <- list(dat=dat, meta.dat=meta.dat)
##    temp <- list(dat=dat)

  }) ## Returns a list of data & metaData
############################### END of inference loop lapply



  #!#names(inference.data) <- treatmentIDs_final
  names(inference.data) <- sample.names

  dat.name <- paste(output.folder, "/",Individual, "_", Plate, "_cov", min.cov, '_inference.RData', sep='')
  save(inference.data, file=dat.name)
  str(inference.data)
  ##load(dat.name)

  ##########################################
  ## 0.) plots for average rho aross the individual
  ##########################################
  all_rho_hat <- Reduce(c, sapply(seq_along(1:length(inference.data)), FUN=function(ii){inference.data[[ii]]$dat$rho}))
  mean_rho_hat <- round(mean(all_rho_hat), 4)
  median_rho_hat <- round(median(all_rho_hat), 4)
  se_rho_hat <- sd(all_rho_hat)

  all_ref <- Reduce(c, sapply(seq_along(1:length(inference.data)), FUN=function(ii){inference.data[[ii]]$dat$refCount}))
  all_alt <- Reduce(c, sapply(seq_along(1:length(inference.data)), FUN=function(ii){inference.data[[ii]]$dat$altCount}))
  all_coverage <- '+'(all_ref, all_alt)
  emp_rho <- '*'(all_ref, all_coverage^(-1))

  rho_title <- paste0(Individual, ' : Rho_hat across all treatements', '\n mean.Rho=', mean_rho_hat, ' | median.Rho=', median_rho_hat)

  pdf.file <- paste0(output.folder, '/plots/QC/', Individual, '_averageRho', '.pdf', sep='')
  pdf(file=pdf.file)
  hist(all_rho_hat[all_coverage>100], main=rho_title, xlim=c(0,1), breaks=seq(0,1,0.01), col='darkgrey', axes=FALSE)
  abline(v=mean_rho_hat, lty=2, col='red')
  axis(1, at=seq(0, 1, .1)); axis(2)
  dev.off()

  ##########################################
  ## 1.) QQ-plots for all treatments
  ##########################################
  for(ii in seq_along(1:length(inference.data))){
    ##ii <- 1
    treatment <- names(inference.data)[ii]
    pvals <- inference.data[[ii]]$dat$pval
    pi0 <- round(inference.data[[ii]]$meta.dat$pi0, 2)
    hets <- inference.data[[ii]]$meta.dat$hets
    qv.2 <- inference.data[[ii]]$meta.dat$qv.2
    coverage <- inference.data[[ii]]$dat$refCount + inference.data[[ii]]$dat$altCount
    avg.depth <- floor(mean(coverage))
    ##disp <- round(mean(inference.data[[ii]]$meta.dat$dispersion), 2)

    ## extract pvalues from only the high coverage loci
    pval_high <- inference.data[[1]]$dat$pval[which(coverage>100)]
    qqp <- qqplot(-log10(ppoints(length(pval_high))),-log10(pval_high), plot.it=F)

    pdf.file <- paste(output.folder, '/plots/QQ/', Individual, '_', treatment, "_cov", min.cov, '_', ii, '_QQ', '.pdf', sep='')
    title <- paste0(Individual, " | ", treatment, " | Pi0=", pi0, " | #hets=", hets, "\n #qv.2=", qv.2, " | avg.depth=", avg.depth) #, " | disp=", disp)
      pdf(file=pdf.file)
    qq(pvals)
    points(qqp,pch='.',cex=5,col='blue')
    title(main=title)
    legend("topleft", c(paste("all hets"), paste("hets with minCov=100")),
           fill=c("black","blue"))
    dev.off()
  }

  ##########################################
  ## 2.) Expression table across all treatments
  ##########################################
  asetable <- t(sapply(seq_along(1:length(inference.data)), FUN=function(ii){
    ##ii <- 1
    sapply(c(.01, .05, .1, .2), FUN=function(jj){sum(inference.data[[ii]]$dat$qv < jj)})
  }))

  rownames(asetable) <- names(inference.data)
  colnames(asetable) <- c('Q<.01', 'Q<.05', 'Q<.1', 'Q<.2')

  outfile <- paste(output.folder, '/',Individual, "_", Plate, '_Qhits.txt', sep='')
  write.table(asetable, file=outfile, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

##})

##                          ##
cat("###### THE END ######")
##                          ##
