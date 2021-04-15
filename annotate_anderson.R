source("/home/sarah/github/targets_arabidopsis/targets_paper_open.R")
library(RNAeditR)

load(file="/binf-isilon/alab/projects/ECT2_TC/parker19/parker19_nanoGR.Rdat")
load(file="/binf-isilon/alab/projects/ECT2_TC/parker19/parker19_micGR.Rdat")
names(nanoGR) <- paste0(as.vector(seqnames(nanoGR)),end(nanoGR))
names(micGR) <- paste0(as.vector(seqnames(micGR)),end(micGR))

quant.vec <- rowMeans(tpm.mat[,grep("c",colnames(tpm.mat))])

nanoGR <- addGenes(gtfGR = gtf,posGR=nanoGR,ncore=30,quant=quant.vec,assignStrand = T,geneids=ids)
micGR <- addGenes(gtfGR = gtf,posGR=micGR,ncore=30,quant=quant.vec,assignStrand = T,geneids=ids)

start(nanoGR) <- end(nanoGR)
start(micGR) <- end(micGR)

save(nanoGR,file="./parker19/parker19_nanoGR.Rdat")
save(micGR,file="./parker19/parker19_micGR.Rdat")

ander <- read.table("/binf-isilon/PBgrp/qbp693/hyperTRIBE_ect3/m6A_shen_andersson/m6A_andersson.bed")
andGR <- GRanges(seqnames=Rle(ander$V1),strand = Rle(ander$V6),ranges=IRanges(ander$V2,ander$V3),name=ander$V4)
save(andGR,file="/binf-isilon/alab/projects/ECT2_TC/anderson_et_al_andGR.Rdat")

load(file="/binf-isilon/alab/projects/ECT2_TC/anderson_et_al_andGR.Rdat")
quant.vec <- rowMeans(tpm.mat[,grep("Sc",colnames(tpm.mat))])
andGR <- addGenes(gtfGR = gtf,posGR=andGR,ncore=30,quant=quant.vec,assignStrand = T,geneids=ids)
save(andGR,file="/binf-isilon/alab/projects/ECT2_TC/anderson_et_al_andGR.Rdat")

plotGeneBody(myGR_list=list("ander"=andGR),list_style = T)
