# Take bed files for each plate and turn into counts file.
# This iis an example using DCM1R1

files<-list.files()
DCM1R1files<-files[grep("_counts.bed.gz", files)]
DCM1R1 <- read.table(DCM1R1files[1], as.is=TRUE, header=FALSE, sep="\t")
DCM1R1counts <- as.data.frame(DCM1R1[,13], row.names = DCM1R1[,4])
for (i in 2:length(DCM1R1files)){
                print(i)
        DCM1R1 <- read.table(DCM1R1files[i], as.is=TRUE, header=FALSE, sep="\t")
        DCM1R1counts<-cbind(DCM1R1counts, DCM1R1[,13])
}
DCM1R1names <- strsplit(DCM1R1files, "_")
DCM1R1names<-as.data.frame(matrix(unlist(DCM1R1names), ncol =2, byrow=TRUE))
colnames(DCM1R1counts)<-DCM1R1names[,1]
write.table(DCM1R1counts, gzfile("DCM1R1.data.gz"),sep="\t",quote=F)
save(DCM1R1counts, file="DCM1R1.data.Rd")
