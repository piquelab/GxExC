#Take bed files for each plate and turn into counts file. Shown for plate CM1R1.
files<-list.files()
CM1R1files<-files[grep("CM1R1-HT", files)]
CM1R1 <- read.table(CM1R1files[1], as.is=TRUE, header=FALSE)
CM1R1counts <- as.data.frame(CM1R1[,15], row.names = CM1R1[,4])
for (i in 2:length(CM1R1files)){
      CM1R1 <- read.table(CM1R1files[i], as.is=TRUE, header=FALSE)
      CM1R1counts<-cbind(CM1R1counts, CM1R1[,15])
      print(i)
}
CM1R1names <- strsplit(CM1R1files, "_")
CM1R1names<-as.data.frame(matrix(unlist(CM1R1names), ncol =2, byrow=TRUE))
colnames(CM1R1counts)<-CM1R1names[,1]
write.table(CM1R1counts, gzfile("CM1R1.data.gz"),sep="\t",quote=F)
save(CM1R1counts, file="CM1R1.data.Rd")
