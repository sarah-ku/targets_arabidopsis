######NEW SET#######
E2_55k <- read.table("/binf-isilon/alab/projects/ECT2_TC/binding sites unfiltered/E2mCh_55kDa_1.bsites.bed",header=F)
E2_55kGR <- GRanges(Rle(as.vector(E2_55k$V1)),strand=Rle(as.vector(E2_55k$V6)),IRanges(E2_55k$V2+5,width=1),score=E2_55k$V5)

E2_55kW <- read.table("/binf-isilon/alab/projects/ECT2_TC/binding sites unfiltered/E2mChW_55kDa_1.bsites.bed",header=F)
E2_55kWGR <- GRanges(Rle(as.vector(E2_55kW$V1)),strand=Rle(as.vector(E2_55kW$V6)),IRanges(E2_55kW$V2+5,width=1),score=E2_55kW$V5)

E2_110k <- read.table("/binf-isilon/alab/projects/ECT2_TC/binding sites unfiltered/E2mCh_110kDa_1.bsites.bed",header=F)
E2_110kGR <- GRanges(Rle(as.vector(E2_110k$V1)),strand=Rle(as.vector(E2_110k$V6)),IRanges(E2_110k$V2+5,width=1),score=E2_110k$V5)

E2_110kW <- read.table("/binf-isilon/alab/projects/ECT2_TC/binding sites unfiltered/E2mChW_110kDa_1.bsites.bed",header=F)
E2_110kWGR <- GRanges(Rle(as.vector(E2_110kW$V1)),strand=Rle(as.vector(E2_110kW$V6)),IRanges(E2_110kW$V2+5,width=1),score=E2_110kW$V5)

seqlevels(E2_55kGR) <- gsub("Chr","",seqlevels(E2_55kGR))
seqlevels(E2_55kWGR) <- gsub("Chr","",seqlevels(E2_55kWGR))
seqlevels(E2_110kGR) <- gsub("Chr","",seqlevels(E2_110kGR))
seqlevels(E2_110kWGR) <- gsub("Chr","",seqlevels(E2_110kWGR))

seqlevels(E2_55kGR) <- c("1","2","3","4","5","Pt")
seqlevels(E2_55kWGR) <- c("1","2","3","4","5","Pt")
seqlevels(E2_110kGR) <- c("1","2","3","4","5","Pt")
seqlevels(E2_110kWGR) <- c("1","2","3","4","5","Pt")

names(E2_55kGR) <- paste(seqnames(E2_55kGR),start(E2_55kGR),sep="_")
names(E2_55kWGR) <- paste(seqnames(E2_55kWGR),start(E2_55kWGR),sep="_")
names(E2_110kGR) <- paste(seqnames(E2_110kGR),start(E2_110kGR),sep="_")
names(E2_110kWGR) <- paste(seqnames(E2_110kWGR),start(E2_110kWGR),sep="_")

iCLIPBS <- list()
iCLIPBS[["E2_55kGR"]] <- E2_55kGR
iCLIPBS[["E2_55kWGR"]] <- E2_55kWGR
iCLIPBS[["E2_110kGR"]] <- E2_110kGR
iCLIPBS[["E2_110kWGR"]] <- E2_110kWGR

#full set
length(iCLIP[[3]])
#binding site set
length(iCLIPBS[[3]])
#old binding site set
length(iclipBS[[1]])

#length of overlap with new binding site set
sum(countOverlaps(iCLIPBS[[3]],iCLIP[[3]]))
#length over overlap with old binding site set
sum(countOverlaps(iclipBS[[1]],iCLIP[[3]]))

#quantifications for hyperTRIBE - use WT samples (shoots+roots) as a proxy for expression in iCLIP
load("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/quantification/transcripts_matrix_araport11.Rdat")

quant.vec <- rowMeans(tpm.mat[,grep("c",colnames(tpm.mat))])

#annotate with pipeline
for(i in 1:4)
{
  iCLIPBS[[i]] <- addGenes(gtfGR = gtf,posGR = iCLIPBS[[i]],ncore=30,quant = quant.vec,assignStrand=T,geneids = ids)
}
length(unique(iCLIPBS[[3]]$gene))

save(iCLIPBS,file="/binf-isilon/alab/projects/ECT2_TC/binding sites unfiltered/iCLIPBS_with_meta.Rdat")
