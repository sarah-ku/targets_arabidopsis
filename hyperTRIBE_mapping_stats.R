filenames <- list.files("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/output/",pattern = "*Log.final.out",full.names = T)
filenames3 <- list.files("/binf-isilon/PBgrp/qbp693/hyperTRIBE_ect3/output/",pattern = "*Log.final.out",full.names = T)

filenames <- c(filenames,filenames3)
#filenames <- gsub(".sort.bam","",filenames)
tur <- rep(0,length(filenames))
inp <- rep(0,length(filenames))
for(i in 1:length(filenames))
{
  inp[i] <- as.vector(read.table(filenames[i],sep="\t",fill = T)[5,]$V2)
  tur[i] <- as.vector(read.table(filenames[i],sep="\t",fill = T)[8,]$V2)
}

#pmr <- as.numeric(gsub("%","",tur))
pmr <- as.numeric(tur)
names(pmr) <- gsub(".+\\/(E[2|3]T_.+)_Log.final.out","\\1",filenames)
tir <- as.numeric(inp)
names(tir) <- gsub(".+\\/(E[2|3]T_.+)_Log.final.out","\\1",filenames)

library(RColorBrewer)
#pdf("input_reads_and_ADAR_expression.PDF",width=3,height=6)
par(mfrow=c(3,1))
par(mar=c(4,4,4,3))
barplot(tir,ylim=c(0,max(tir)),ylab="number of input reads \n (millions)",las=2,col=brewer.pal(4,"Set1")[rep(1:4,each=3)])
barplot(pmr,ylim=c(0,max(pmr)),ylab="number of input reads \n (millions)",las=2,col=brewer.pal(4,"Set1")[rep(1:4,each=3)])

read_stats <- data.frame("Input reads (millions)"=tir,"Uniquely mapped (millions)"=pmr)
read_stats <- read_stats[-grep("E2T_St5",row.names(read_stats)),]
row.names(read_stats)[grep("E2T_Se5",row.names(read_stats))] <- "E2T_St5"
read_stats <- read_stats[-grep("St|Rt",row.names(read_stats)),]

write.table(read_stats,"/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/ECT2_3_read_statistics.txt",quote=F,sep="\t")

