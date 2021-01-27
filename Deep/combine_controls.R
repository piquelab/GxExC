# Anthony Findley
# 9/21/18

# Purpose: Combine control replicate pileups on the same plate
# homeDir <- "/nfs/rprdata/Anthony/GxE_iPSC/Deep/"

cargs <- commandArgs(trail=TRUE)
if(length(cargs)>=1)
  Plate <- cargs[ 1 ]

pair1 <- c("HT85", "HT87", "HT89", "HT91", "HT93", "HT95") # These 6 barcodes must be combined. pair1[x] goes with pair2[x], where x=1:6
pair2 <- c("HT86", "HT88", "HT90", "HT92", "HT94", "HT96")


for (i in 1:length(pair1)){

        # Load first replicate
        R1 <- read.table( paste0("/nfs/rprdata/Anthony/GxE_iPSC/Deep/", Plate, "/P1_R2/Pileups/", Plate, "-", pair1[i], ".pileup.clean.bed.gz"), header=F, stringsAsFactors=F, sep="\t")
        colnames(R1) <- c("chr", "pos0", "pos1", "ref", "alt", "rsID", "alt_AF", "ref_reads", "alt_reads", "error")

        R2 <- read.table( paste0("/nfs/rprdata/Anthony/GxE_iPSC/Deep/", Plate, "/P1_R2/Pileups/", Plate, "-", pair2[i], ".pileup.clean.bed.gz"), header=F, stringsAsFactors=F, sep="\t")
        colnames(R2) <- c("chr", "pos0", "pos1", "ref", "alt", "rsID", "alt_AF", "ref_reads", "alt_reads", "error")

        merged <- merge(R1, R2, by=c("chr", "pos0", "pos1", "ref", "alt", "rsID", "alt_AF"), all=T) # Inserts NA's for SNPs present in one file but not both

        # Change NA's for ref, alt, and error read counts to 0
        merged$ref_reads.x[which(is.na(merged$ref_reads.x))] <- 0
        merged$ref_reads.y[which(is.na(merged$ref_reads.y))] <- 0
        merged$alt_reads.x[which(is.na(merged$alt_reads.x))] <- 0
        merged$alt_reads.y[which(is.na(merged$alt_reads.y))] <- 0
        merged$error.x[which(is.na(merged$error.x))] <- 0
        merged$error.y[which(is.na(merged$error.y))] <- 0

        merged$total_ref <- merged$ref_reads.x + merged$ref_reads.y
        merged$total_alt <- merged$alt_reads.x + merged$alt_reads.y
        merged$total_error <- merged$error.x + merged$error.y

        final_tab <- merged[, c("chr", "pos0", "pos1", "ref", "alt", "rsID", "alt_AF", "total_ref", "total_alt", "total_error")]

        system( paste0( "mkdir -p ./combined_pileups/", Plate) )
        write.table(final_tab, file=paste0("./combined_pileups/", Plate, "/", Plate, "-", pair1[i], ".pileup.clean.bed"), quote=F, sep="\t", col.names=F, row.names=F)
        system( paste0("gzip ./combined_pileups/", Plate, "/", Plate, "-", pair1[i], ".pileup.clean.bed") )
}
