source("/home/sarah/github/targets_arabidopsis/targets_paper_open.R")
source("/home/sarah/github/targets_arabidopsis/motifs_functions.R")

setwd("/binf-isilon/alab/people/sarah/meme/motif_plots")

library(doParallel)
registerDoParallel(cores=30)
library(ggplot2)
library(pheatmap)
library(rtracklayer)
library(GenomicRanges)
library(RColorBrewer)
library("BSgenome.Athaliana.TAIR.TAIR9")
seqnames(Athaliana) <- c("1","2","3","4","5","Mt","Pt")


load(file="/binf-isilon/alab/people/sarah/meme/matched_sets.Rdat")
load(file="/binf-isilon/alab/people/sarah/meme/motGR.Rdat")
load(file="/binf-isilon/alab/people/sarah/meme/randGR.Rdat")
load(file="/binf-isilon/alab/people/sarah/meme/ucontent.Rdat")

gtfE <- gtf

#############################################
########MATCH SETS PLOT BY POSITION##########
#############################################
rmg <- unique(unlist(lapply(pos.list.ECT2,function(x) x$gene)))
nontargs <- genes.list$genes.nontargets
summary(nontargs %in% rmg)
nontargs <- nontargs[!(nontargs %in% rmg)]

rGR <- randGR
rGR <- rGR[which(countOverlaps(rGR,gtfGR[gtfGR$gene_id %in% nontargs])>0)]
names(rGR) <- paste0(seqnames(rGR),":",start(rGR),",",strand(rGR))

pdf("motif_matched_location_sets_fg_bg_all_tol_25.pdf",width=144,height=20)
nset <- length(matchGR)-1
par(mfcol=c(nset,length(motGR)))
glmdat_list <- list()
rat_list <- list()
for(k in 1:length(motGR))
{
  getRun <- getCountDens(matchGR,motGR[[k]],tol = 25)
  fgbg_list <- getRun[[1]]
  glmdat_list[[k]] <- getRun[[3]]
  
  rat <- rep(0,nset)
  for(i in 1:nset)
  {
    tp <- fgbg_list[[i]][,1]
    tpx <- seq(0,3,by=0.1)
    tpy <- rep(0,length(tpx))
    names(tpy) <- paste(tpx)
    tpy[names(tp)] <- tp
    plot(tpx,tpy,type="l",col="black",lwd=2,main=paste(names(motGR)[k],names(fgbg_list)[[i]]),ylim=c(0,max(fgbg_list[[i]])),ylab="Number of motif",xlab="")
    
    abline(v=c(1,2),lty=2)
    
    tp <- fgbg_list[[i]][,2]
    tpx <- seq(0,3,by=0.1)
    tpy <- rep(0,length(tpx))
    names(tpy) <- paste(tpx)
    tpy[names(tp)] <- tp
    
    points(tpx,tpy,type="l",col="dark grey",lty=1,lwd=2)
    rat[i] <- sum(fgbg_list[[i]][,1])/sum(fgbg_list[[i]][,2])
  }
  rat_list[[k]] <- rat
}
dev.off()


######################################################################
###############RRACH+ and RRACH- enrichments     #####################
######################################################################

nmGR <- matchGR$nano
nmGR <- matchGR$random[matchGR$nano$match]
nmGR_list <- list(nanoGR,matchGR$random[matchGR$nano$match])
names(nmGR_list) <- c("nanopore","matched control sites")

ucontent_nr <- ucontent[-which(names(ucontent)=="RRACH")]
ucontent_nr[order(ucontent_nr)]
pdf("nanopore_RRACHp_vs_RRACHm_10.pdf",width=10,height=5)

par(mfrow=c(1,1))
rrn <- which(names(motGR)=="RRACH")
nmGR <- nmGR_list[[1]]
olcol <- do.call(cbind,lapply(motGR,function(x) (countOverlaps(nmGR+10,x))))
olcol[olcol>0] <- 1
olcol <- as.data.frame(olcol)
olcol[,"none"] <- 0
olcol[,"none"][rowSums(olcol[,c(1:ncol(olcol))[-rrn]])==0] <- 1
colSums(olcol)

#sum(rowSums(olcol[,c(1:6,8)][olcol[,rrn]==0,])==0)

rrNc <- colSums(olcol[olcol[,"RRACH"]==0,])[-rrn]/sum(olcol[,"RRACH"]==0)
rrPc <- colSums(olcol[olcol[,"RRACH"]==1,])[-rrn]/sum(olcol[,"RRACH"]==1)

sum(olcol[,"RRACH"]==0)
colSums(olcol[olcol[,"RRACH"]==0,])[-rrn]

sum(olcol[,"RRACH"]==1)
colSums(olcol[olcol[,"RRACH"]==1,])[-rrn]

nmGR <- nmGR_list[[2]]
olcol <- do.call(cbind,lapply(motGR,function(x) (countOverlaps(nmGR+10,x))))
olcol[olcol>0] <- 1
olcol <- as.data.frame(olcol)
olcol[,"none"] <- 0
olcol[,"none"][rowSums(olcol[,c(1:ncol(olcol))[-rrn]])==0] <- 1
rrNcbg <- colSums(olcol[olcol[,"RRACH"]==0,])[-rrn]/sum(olcol[,"RRACH"]==0)
rrPcbg <- colSums(olcol[olcol[,"RRACH"]==1,])[-rrn]/sum(olcol[,"RRACH"]==1)
# rrNc <- rrNc/sum(rrNc)
# rrPc <- rrPc/sum(rrPc)

nbg <- (rrNc/rrNcbg)
pbg <- (rrPc/rrPcbg)
tp <- rbind(nbg,pbg)
#tp[,1:(length(nbg)-1)] <- tp[,names(ucontent_nr)]

barplot(tp,las=2,ylim=c(0,max(nbg)),beside = T,col=brewer.pal(3,"Set1")[1:2],ylab="enrichment VS background count")
abline(h=1,lty=2)
legend("topleft",c("RRACH-","RRACH+"),col=brewer.pal(3,"Set1")[1:2],lwd=3)
#}

dev.off()

for(i in 1:length(rrNc))
{
  tab <- round(length(matchGR$nano)*rbind(c(rrNc[i],rrNcbg[i]),c(rrPc[i],rrPcbg[i])))
  print(paste(names(rrNc)[i],fisher.test(tab)$p.value))
}


######################################################################
###############RRACH+ and RRACH- motif locations #####################
######################################################################

nmGR <- matchGR$nano
nmGR <- matchGR$random[matchGR$nano$match]
nmGR_list <- list(matchGR$nano,matchGR$random[matchGR$nano$match])
names(nmGR_list) <- c("nanopore","matched control sites")

aSet1 <- nmGR_list[[1]][countOverlaps(nmGR_list[[1]]+10,motGR[["RRACH"]])>0]
aSet2 <- nmGR_list[[1]][countOverlaps(nmGR_list[[1]]+10,motGR[["RRACH"]])==0]
aSet3 <- nmGR_list[[2]][countOverlaps(nmGR_list[[2]]+10,motGR[["RRACH"]])>0]
aSet4 <- nmGR_list[[2]][countOverlaps(nmGR_list[[2]]+10,motGR[["RRACH"]])==0]

set.names <- c("RRACH_+","RRACH_-","bg_+","bg_-")

colset <- brewer.pal(8,"Set1")

pSets <- list(aSet1,aSet2,aSet3,aSet4)

pdf("motif_locations_relative_m6A_sets_RRACHP_RRACHM.pdf",width=3*length(motGR),height=2.6*2)
par(mfrow=c(2,length(motGR)))
for(i in 1:length(motGR))
{
  makeRelPlotSets(motGR[[i]],pSets,myn=200,my.col=colset,my.title=paste0(mot_name,"_nanopore"),use.para=TRUE)
  legend("topright",set.names,col=colset[1:4],lwd=2)
}
dev.off()

for(i in 1:length(iCLIP))
{
  makeRelPlotSets(iCLIP[[i]],pSets,myn=200,my.col=mycols,my.title=paste0(mot_name,"_nanopore"),use.para=TRUE)
}

############iCLIP

#pdf("iclip_locations_relative_m6A_sets_RRACHP_RRACHM.pdf",width=38,height=6)
#mvec.ud <- c(20,5,100,10,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
#par(mfrow=c(2,length(motGR)))
#for(i in 1:length(motGR))
#{
#  makeDBplot(aSet_list=list(aSet1,aSet2,aSet3,aSet4),motGR=iCLIP[[i]],colset=colset,my.name=names(iCLIP)[i],y.max=mvec.ud[i],my.set.names = set.names)
#}
#dev.off()

# pdf("raw_motif_overlap_iclip_status.pdf",width=10,height=8)
# par(mfrow=c(3,6))
# for(i in 1:18)
# {
#   wol <- (unlist(lapply(iCLIP,function(x) sum(countOverlaps(x+10,motGR[[i]])>0)/length(x))))
#   names(wol) <- c("55","55W","110","110W")
#   barplot(wol,las=2,col=colset[1:2],main=names(motGR)[i],ylab="prop. with motif")
# }
# dev.off()


############iCLIP
iCLIP_W <- lapply(iCLIP,function(x) x[countOverlaps(x+10,motGR[["WEI"]])>0])
iCLIP_NW <- lapply(iCLIP,function(x) x[countOverlaps(x+10,motGR[["WEI"]])==0])
GList <- c((iCLIP_W),(iCLIP_NW))
length(GList)

pdf("iclip_with_WEI_locations_relative_m6A_sets_RRACHP_RRACHM.pdf",width=38,height=6)
mvec.ud <- c(10,10,30,15,20,25,100,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
par(mfrow=c(2,length(motGR)))
for(i in 1:length(motGR))
{
  makeDBplot(aSet_list=list(aSet1,aSet2,aSet3,aSet4),motGR=GList[[i]],colset=colset,my.name=names(GList)[i],y.max=mvec.ud[i],my.set.names = set.names)
}
dev.off()

#######


# aSet1 <- nmGR_list[[1]][countOverlaps(nmGR_list[[1]]+10,motGRR[["RRACH"]])>0]
# aSet2 <- nmGR_list[[1]][countOverlaps(nmGR_list[[1]]+10,motGRR[["RRACH"]])==0]
# aSet3 <- nmGR_list[[2]][countOverlaps(nmGR_list[[2]]+10,motGRR[["RRACH"]])>0]
# aSet4 <- nmGR_list[[2]][countOverlaps(nmGR_list[[2]]+10,motGRR[["RRACH"]])==0]
# GLset <- list("R+"=aSet1,"R-"=aSet2,"bg+"=aSet3,"bg-"=aSet4)
# 
# 
# pdf("m6A_wei_locations_relative_iCLIP.pdf",width=38,height=6)
# mvec.ud <- c(100,100,20,15,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
# par(mfrow=c(2,length(motGR)))
# for(i in 1:length(motGR))
# {
#   makeDBplot(aSet_list=iCLIP,motGR=GLset[[i]],colset=colset,my.name=names(GLset)[i],y.max=mvec.ud[i],my.set.names = names(iCLIP))
# }
# dev.off()



##########
##########

mnames <- c("RRACH","AVAYU")
pdf("iclip_locations_relative_m6A_sets_RRACH_URUAY.pdf",width=10,height=4)
par(mfrow=c(1,3))
for(nmGR in nmGR_list)
{
  
  olcol <- do.call(cbind,lapply(motGR,function(x) (countOverlaps(nmGR+10,x))))
  olcol[olcol>0] <- 1
  olcol <- as.data.frame(olcol)
  olcol$empty <- rep(1,nrow(olcol))
  olcol$empty[rowSums(olcol)==1] <- 0
  
  
  m6aRR <- olcol[,mnames[1]]
  m6aRR[m6aRR>0] <- 1
  
  nanoPLUS <- nmGR[which(m6aRR==1)]
  nanoMINUS <- nmGR[m6aRR==0]
  
  m6aRR <- olcol[,mnames[2]]
  m6aRR[m6aRR>0] <- 1
  
  nanoPLUS_U <- nmGR[m6aRR==1]
  nanoMINUS_U <- nmGR[m6aRR==0]
  
  colset <- RColorBrewer::brewer.pal(4,"Set1")
  makeRelPlot(setGR[[1]],nanoMINUS_U,to.add=F,my.col=colset[4],my.ylim=c(0,30))
  makeRelPlot(setGR[[1]],nanoPLUS_U,to.add=T,my.col=colset[3])
  makeRelPlot(setGR[[1]],nanoMINUS,to.add=T,my.col=colset[2])
  makeRelPlot(setGR[[1]],nanoPLUS,to.add=T,my.col=colset[1])
  lengths_ta <- c(length(nanoPLUS),length(nanoMINUS),length(nanoPLUS_U),length(nanoMINUS_U))
  leg.names <- paste(rep(mnames,each=2),c("+","-","+","-"),lengths_ta)
  #leg.names
  legend("topright",leg.names,col=colset,lty=1)
  
}
dev.off()


setGR1 <- setGR[[1]]
setGR2 <- nanoMINUS



##########################################################################
################### motifs centered on m6A ###############################
##########################################################################
# colset <- brewer.pal(8,"Set1")
# 
# aSet <- matchGR$wei
# aSet1 <- aSet
# aSet2 <- matchGR$random[aSet$match]
# 
# set.names <- c("actual","background")
# 
# pdf("wei_VS_motif_matched_background.pdf",width=20,height=6)
# mvec.ud <- c(120,15,65,40,200,100,350,250)
# par(mfrow=c(2,8))
# for(i in 1:8)
# {
#   makeDBplot(aSet_list=list(aSet1,aSet2),motGR=motGRR[[i]],colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names = set.names)
# }
# dev.off()
# 
# aSet <- matchGR$mic
# aSet1 <- aSet[aSet$gene %in% genes.list$genes.union]
# aSet2 <- aSet[aSet$gene %in% genes.list$genes.nontargets]
# #aSet3 <- aSet[aSet$gene %in% genes.list$genes.intersect]
# #aSet4 <- aSet[(aSet$gene %in% genes.list$genes.nontargets) & (aSet$gene %in% up_genes)]
# 
# aSet5 <- matchGR$random[aSet$match]
# 
# pdf("m6A_sites_mic_VS_motif_targs_matched_background.pdf",width=20,height=6)
# 
# set.names <- c("actual-targs","act-nontargs","background")
# mvec.ud <- c(160,15,85,40,200,100,400,250)
# par(mfrow=c(2,8))
# for(i in 1:8)
# {
#   makeDBplot(aSet_list=list(aSet1,aSet2,aSet5),motGR=motGRR[[i]],colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names = set.names)
# }
# dev.off()

##########################################################################
################### open down and up regualted genes ###############################
##########################################################################


load(paste0("/binf-isilon/PBgrp/qbp693/ect2_protoplast/output_rnaseq_smartseq_analyses/pf2_de_output_smartseq.Rdata"))

down_genes <- resLFC_pf2_df_smartseq[resLFC_pf2_df_smartseq$log2FoldChange<(-0.5),]$gene_id
up_genes <- resLFC_pf2_df_smartseq[resLFC_pf2_df_smartseq$log2FoldChange>0.5,]$gene_id
nc_genes <- resLFC_pf2_df_smartseq[abs(resLFC_pf2_df_smartseq$log2FoldChange)<0.5,]$gene_id


down_expr <- (expr.list_ECT2$roots.expr[down_genes] + expr.list_ECT2$shoots.expr[down_genes])/2
up_expr <- (expr.list_ECT2$roots.expr[up_genes] + expr.list_ECT2$shoots.expr[up_genes])/2

t.test(log2(down_expr+1),log2(up_expr+1))
plot(density(log2(down_expr+1)))
points(density(log2(up_expr+1))$x,density(log2(up_expr+1))$y,type="l",col="blue")

#m6aRR <- (countOverlaps(nmGR+6,motG1))

##########################################################################
################### motifs centered on m6A ###############################
##########################################################################
# colset <- brewer.pal(8,"Set1")
# 
# #aSet <- nanoGR[countOverlaps(nanoGR+10,micGR)>0]
# aSet <- nanoGR
# aSet1 <- aSet[aSet$gene %in% up_genes]
# aSet2 <- aSet[aSet$gene %in% nc_genes]
# aSet3 <- aSet[aSet$gene %in% down_genes]
# set.names <- c("up_genes","no-change","down_genes")
# 
# 
# pdf("m6A_sites_VS_motifs_up_down.pdf",width=38,height=6)
# mvec.ud <- c(120,10,20,15,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
# #par(mfrow=c(2,2))
# par(mfrow=c(2,length(motGRR)))
# for(i in 1:length(motGRR))
# {
#   makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=motGR[[i]],colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names = set.names)
# }
# dev.off()
# 
# 
# #aSet <- nanoGR[countOverlaps(nanoGR+25,micGR)>0]
# aSet <- nanoGR
# aSet1 <- aSet[aSet$gene %in% genes.list$genes.nontargets]
# aSet2 <- aSet[aSet$gene %in% genes.list$genes.union[!(genes.list$genes.union %in% genes.list$genes.intersect)]]
# aSet3 <- aSet[aSet$gene %in% genes.list$genes.intersect]
# set.names <- c("non_targets","target(non-strict)","strict")
# 
# pdf("m6A_sites_VS_motifs_targets_non_targets.pdf",width=38,height=6)
# mvec.ud <- c(120,10,20,15,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
# par(mfrow=c(2,length(motGR)))
# for(i in 1:length(motGR))
# {
#   makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=motGR[[i]],colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names=set.names)
# }
# dev.off()
# 
# pdf("m6A_sites_VS_motifs_targets_non_targets.pdf",width=20,height=6)
# mvec.ud <- c(120,10,20,15,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
# par(mfrow=c(2,8))
# 
# aSet <- nanoGR[countOverlaps(nanoGR+25,micGR)>0]
# aSet1 <- aSet[aSet$gene %in% genes.list$genes.nontargets]
# aSet2 <- aSet[aSet$gene %in% genes.list$genes.union[!(genes.list$genes.union %in% genes.list$genes.intersect)]]
# aSet3 <- aSet[aSet$gene %in% genes.list$genes.intersect]
# set.names <- c("non_targets","target(non-strict)","strict")

#sum(countOverlaps(gtf[gtf$gene_id %in%],motGR[["WEI"]]))/length(aSet1)
# motGR[["WEI"]]
# library(RNAeditR)
# library(foreach)
# library(doParallel)
# #library()
# motWEIGR <- addGenes(gtfGR=gtfGR,posGR = motGR[["WEI"]],ncore=30,quant = rowMeans(tpm.mat[,grep("c",colnames(tpm.mat))]),geneids = ids,assignStrand = F)
# 
# ol <- countOverlaps(motGR[["WEI"]],gtfGR[gtfGR$gene_id %in% genes.list$genes.intersect])
# sum(ol)/length(ol)
# 
# gene.set <- sample(genes.list$genes.nontargets,length(genes.list$genes.intersect))
# ol <- countOverlaps(motGR[["WEI"]],gtfGR[gtfGR$gene_id %in% gene.set])
# sum(ol)/length(ol)
# 
# 
# makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=iCLIP[[3]],colset=colset,my.name="iCLIP",y.max=300,levels=set.names)
# makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=setGR[[4]],colset=colset,my.name="iCLIP",y.max=300,levels=set.names)
# 
# 
# aSet <- nanoGR[countOverlaps(nanoGR+25,micGR)>0]
# aSet1 <- aSet[aSet$gene %in% up_genes]
# aSet2 <- aSet[aSet$gene %in% nc_genes]
# aSet3 <- aSet[aSet$gene %in% down_genes]
# set.names <- c("up_genes","no-change","down_genes")
# 
# makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=iCLIP[[3]],colset=colset,my.name="iCLIP",y.max=300,levels=set.names)
# makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=setGR[[4]],colset=colset,my.name="iCLIP",y.max=300,levels=set.names)
# 
# dev.off()


####################################
olcol <- do.call(cbind,lapply(motGRR,function(x) (countOverlaps(aSet+50,x))))
olcol[olcol>0] <- 1
olcol <- as.data.frame(olcol)
olcol$resp <- c("none")
olcol$resp[aSet$gene %in% up_genes] <- "up"
olcol$resp[aSet$gene %in% down_genes] <- "down"
olcol2 <- olcol[olcol$resp %in% c("up","down"),]
olcol2$resp <- as.numeric(factor(olcol2$resp))-1
head(olcol2)
myfit <- (glm(resp~(.)^2,data=olcol2,family="binomial"))

plot(roc(olcol2$resp,myfit$fitted))
roc(olcol2$resp,myfit$fitted)

ep <- expand.grid(names(motGR),names(motGR))
resp.nos <- apply(ep,1,function(x) olcol2$resp[rowSums(olcol2[,x])==3])
rprop <- unlist(lapply(resp.nos,function(x) mean(x)))
names(rprop) <- apply(ep,1,function(x) paste(x,collapse=";"))
barplot((sort(rprop)),las=2)
plot(density(myfit$fitted.values[olcol2$resp==1]))
md <- density(myfit$fitted.values[olcol2$resp==0])
points(md$x,md$y,type="l",col="red")




#########Straight iCLIP VS nano etc.
length(unique(matchGR$iclip$gene))




makeDBplot(aSet_list=aSet_list,motGR=motGR[["WEI"]],colset=colset,my.name="WEI",y.max=100,my.set.names = set.names)

aSet_list <- c("iclip"=matchGR$iclip)
makeDBplot(aSet_list=aSet_list,motGR=matchGR$nano,colset=colset,my.name="nanoPore",y.max=100,my.set.names = set.names)
aSet_list <- c("nano"=matchGR$nano)
makeDBplot(aSet_list=aSet_list,motGR=matchGR$iclip,colset=colset,my.name="nanoPore",y.max=100,my.set.names = set.names)


aSet_list <- c("nano"=(matchGR$nano))
makeDBplot(aSet_list=aSet_list,motGR=motGR2[["RRACH"]],colset=colset,my.name="nano",y.max=1100,my.set.names = "",to.smooth = T)
makeDBplot(aSet_list=aSet_list,motGR=motGR2[["RRACH"]],colset=colset,my.name="nano",y.max=1100,my.set.names = "",to.smooth = F)

par(mfrow=c(2,2))
aSet_list <- c("nano"=(matchGR$nano))
makeDBplot(aSet_list=aSet_list,motGR=motGR[["RRACH"]],colset=colset,my.name="nano",y.max=1100,my.set.names = "",to.smooth = F)
makeDBplot(aSet_list=aSet_list,motGR=motGR2[["WEI"]],colset=colset,my.name="nano",y.max=350,my.set.names = "",to.smooth = F)

aSet_list <- c("nano"=(matchGR$iclip))
makeDBplot(aSet_list=aSet_list,motGR=motGR2[["RRACH"]],colset=colset,my.name="nano",y.max=400,my.set.names = "",to.smooth = F)
makeDBplot(aSet_list=aSet_list,motGR=motGR2[["WEI"]],colset=colset,my.name="nano",y.min = 100,y.max=350,my.set.names = "",to.smooth = F)



makeRelPlot(motGR[["RRACH"]],matchGR$nano,myn=200,to.add=F,my.col="black",my.ylim=c(0,1000),my.title="",to.smooth=F,use.para = TRUE)
makeRelPlot(motGR[["RRACH"]],matchGR$iclip,myn=200,to.add=F,my.col="black",my.ylim=c(0,400),my.title="",to.smooth=F,use.para = TRUE)

library(doParallel)
registerDoParallel(cores=40)
asizes <- c(10,8,8,8,8,5,7,15,10,8,12,35,20,40,40,50,20,10,50)
length(asizes)
names(asizes) <- names(motGR)
for(i in 1:length(motGR))
{
  pdf(paste0("binding_sites_check_",names(motGR)[i],".pdf"),width=8,height=8)
  par(mfrow=c(2,2))
  makeRelPlot(motGR[[i]],iCLIPBS[[3]],myn=200,to.add=F,my.col="black",my.ylim=c(0,asizes[names(motGR)[i]]),my.title="Binding sites",to.smooth=F,use.para = TRUE)
  makeRelPlot(motGR[[i]],iCLIP[[3]],myn=200,to.add=F,my.col="black",my.ylim=c(0,asizes[names(motGR)[i]]),my.title="Full set",to.smooth=F,use.para = TRUE)
  extraBS <- iCLIPBS[[3]][countOverlaps(iCLIPBS[[3]],iclipBS[[1]])==0]
  makeRelPlot(motGR[[i]],extraBS,myn=200,to.add=F,my.col="black",my.ylim=c(0,asizes[names(motGR)[i]]),my.title="Binding sites extra",to.smooth=F,use.para = TRUE)
  extraBS <- iCLIP[[3]][countOverlaps(iCLIP[[3]],iCLIPBS[[3]])==0]
  makeRelPlot(motGR[[i]],extraBS,myn=200,to.add=F,my.col="black",my.ylim=c(0,asizes[names(motGR)[i]]),my.title="Full set extra",to.smooth=F,use.para = TRUE)
  dev.off()
}


htl <- unique(unlist(GRangesList(pos.list.ECT2)))
setList <- list("nano"=nanoGR,"mic"=micGR,"ht"=htl,"wei"=weiGRC)
asizes <- c(250,1000,20,22)
for(i in 3)
{
  pdf(paste0("binding_sites_check_",names(setList)[i],".pdf"),width=8,height=8)
  par(mfrow=c(2,2))
  makeRelPlot(setList[[i]],iCLIPBS[[3]],myn=200,to.add=F,my.col="black",my.ylim=c(0,asizes[i]),my.title="Binding sites",to.smooth=F,use.para = TRUE)
  makeRelPlot(setList[[i]],iCLIP[[3]],myn=200,to.add=F,my.col="black",my.ylim=c(0,asizes[i]),my.title="Full set",to.smooth=F,use.para = TRUE)
  extraBS <- iCLIPBS[[3]][countOverlaps(iCLIPBS[[3]],iclipBS[[1]])==0]
  makeRelPlot(setList[[i]],extraBS,myn=200,to.add=F,my.col="black",my.ylim=c(0,asizes[i]),my.title="Binding sites extra",to.smooth=F,use.para = TRUE)
  extraBS <- iCLIP[[3]][countOverlaps(iCLIP[[3]],iCLIPBS[[3]])==0]
  makeRelPlot(setList[[i]],extraBS,myn=200,to.add=F,my.col="black",my.ylim=c(0,asizes[i]),my.title="Full set extra",to.smooth=F,use.para = TRUE)
  dev.off()
}






pdf(paste0("binding_sites_check_experiments.pdf"),width=8,height=8)
par(mfrow=c(2,2))
makeRelPlot(nanoGR,iCLIP[[3]],myn=200,to.add=F,my.col="black",my.ylim=c(0,250),my.title="Binding sites",to.smooth=F,use.para = TRUE)
makeRelPlot(weiGRC,iCLIP[[3]],myn=200,to.add=F,my.col="black",my.ylim=c(0,8),my.title="Full set",to.smooth=F,use.para = TRUE)
#extraBS <- iCLIPBS[[3]][countOverlaps(iCLIPBS[[3]],iclipBS[[1]])==0]
makeRelPlot(micGR,iCLIP[[3]],myn=200,to.add=F,my.col="black",my.ylim=c(0,900),my.title="Binding sites extra",to.smooth=F,use.para = TRUE)
makeRelPlot(pos.list.ECT2[[1]],iCLIP[[3]],myn=200,to.add=F,my.col="black",my.ylim=c(0,22),my.title="Full set extra",to.smooth=F,use.para = TRUE)
dev.off()

extraBS <- iCLIP[[3]][countOverlaps(iCLIP[[3]],iCLIPBS[[3]])==0]


sma <- tapply(start(iCLIP[[3]][strand(iCLIP[[3]])=="+"]),iCLIP[[3]][strand(iCLIP[[3]])=="+"]$gene,mean)
smab <- tapply(start(iCLIPBS[[3]][strand(iCLIPBS[[3]])=="+"]),iCLIPBS[[3]][strand(iCLIPBS[[3]])=="+"]$gene,mean)
table(sign(sma-smab))

sma <- tapply(start(iCLIP[[3]][strand(iCLIP[[3]])=="-"]),iCLIP[[3]][strand(iCLIP[[3]])=="-"]$gene,mean)
smab <- tapply(start(iCLIPBS[[3]][strand(iCLIPBS[[3]])=="-"]),iCLIPBS[[3]][strand(iCLIPBS[[3]])=="-"]$gene,mean)
table(sign(sma-smab))

head(sma)
head(smab)
summary( (sma-smab)/1000 )   
hist((sma-smab)/1000 ,breaks=1000,xlim=c(-2,2))
table(sign(sma-smab))

#if strand is negative, then smab >  sma (implies binding site is more upstream)
#if strand is positive, then sma > smab (implies binding site is more upstream)

extraBS <- iCLIP[[3]][countOverlaps(iCLIP[[3]],iCLIPBS[[3]])==0]
makeRelPlot(iCLIPBS[[3]],extraBS,myn=10,to.add=F,my.col="black",my.ylim=c(0,1000),my.title="Full set",to.smooth=F,use.para = TRUE)
makeRelPlot(iCLIPBS[[3]],iCLIP[[3]],myn=10,to.add=F,my.col="black",my.ylim=c(0,1000),my.title="Full set",to.smooth=F,use.para = TRUE)
extraBS <- iCLIPBS[[3]][countOverlaps(iCLIPBS[[3]],iclipBS[[1]])==0]
makeRelPlot(extraBS,iCLIP[[3]],myn=10,to.add=F,my.col="black",my.ylim=c(0,1000),my.title="Full set",to.smooth=F,use.para = TRUE)









pdf("enrichment_of_m6A_by_iCLIP_score.pdf",width=4,height=4)
par(mfrow=c(1,1))
makeRelPlot(setList[["nano"]],iCLIP[[3]][iCLIP[[3]]$score<1],myn=50,to.add=F,my.col="black",my.ylim=c(0,350),my.title="Binding sites",to.smooth=F,use.para = TRUE)
makeRelPlot(setList[["nano"]],iCLIP[[3]][iCLIP[[3]]$score>1 & iCLIP[[3]]$score<4],myn=50,to.add=T,my.col="blue",my.ylim=c(0,350),my.title="Full set",to.smooth=F,use.para = TRUE)
makeRelPlot(setList[["nano"]],iCLIP[[3]][iCLIP[[3]]$score>4 & iCLIP[[3]]$score<10],myn=50,to.add=T,my.col="red",my.ylim=c(0,350),my.title="Full set",to.smooth=F,use.para = TRUE)
makeRelPlot(setList[["nano"]],iCLIP[[3]][iCLIP[[3]]$score>4 & iCLIP[[3]]$score>10],myn=50,to.add=T,my.col="green",my.ylim=c(0,350),my.title="Full set",to.smooth=F,use.para = TRUE)
legend("topleft",c("score<1","1<score<4","4<score<10","score>10"),col=c("black","blue","red","green"),lwd=1)  
dev.off()


pdf("enrichment_of_m6A_by_iCLIP_peakset.pdf",width=4,height=4)
par(mfrow=c(1,1))
makeRelPlot(setList[["nano"]],iCLIPBS[[3]],myn=100,to.add=F,my.col="black",my.ylim=c(0,270),my.title="Binding sites",to.smooth=F,use.para = TRUE)
makeRelPlot(setList[["nano"]],iCLIP[[3]],myn=100,to.add=T,my.col="blue",my.ylim=c(0,270),my.title="Full set",to.smooth=F,use.para = TRUE)
extraBS <- iCLIP[[3]][countOverlaps(iCLIP[[3]],iCLIPBS[[3]])==0]
makeRelPlot(setList[["nano"]],extraBS,myn=100,to.add=T,my.col="red",my.ylim=c(0,270),my.title="Full set extra",to.smooth=F,use.para = TRUE)
legend("topleft",c("binding sites (8808)","full set (15960)","extra peaks (7156)"),col=c("black","blue","red"),lwd=1)
dev.off()

ms <- tapply(iCLIP[[3]]$score,iCLIP[[3]]$gene,mean)
msb <- tapply(iCLIPBS[[3]]$score,iCLIPBS[[3]]$gene,mean)

plot(ms,msb)  
abline(0,1)  



extraBS <- iCLIP$E2_110kGR[countOverlaps(iCLIP$E2_110kGR,iCLIPBS$E2_110kGR)==0]
oli <- countOverlaps(extraBS+5,motGR[["WEI"]]) #7156 extra peaks which are not binding sites
olibs <- countOverlaps(iCLIPBS$E2_110kGR+5,motGR[["WEI"]]) #8808 peaks which are binding sites

mean(oli)
mean(olibs)

mean(oli)/mean(olibs) #ratio
mean(oli>0)/mean(olibs>0) #ratio of peaks with at least 1 m6A






hist(oli[sample(1:length(oli),length(olibs),replace = T)],col="pink",ylim=c(0,5000))
hist(olibs,add=T,col="yellow")


gtfE <- gtf





setList <- list("nano"=nanoGR,"mic"=micGR,"ht"=htl,"wei"=weiGRC,"iclip"=matchGR$iclip)
asizes <- c(250,1000,20,22)

pdf(paste0("binding_sites_check_",names(setList)[i],".pdf"),width=8,height=8)
par(mfrow=c(2,2))
makeRelPlot(setList[[2]],setList[[1]],myn=100,to.add=F,my.col="black",my.ylim=c(0,asizes[i]),my.title="Binding sites",to.smooth=F,use.para = FALSE)
makeRelPlot(setList[[3]],setList[[1]],myn=100,to.add=T,my.col="black",my.ylim=c(0,asizes[i]),my.title="Full set",to.smooth=F,use.para = FALSE)
makeRelPlot(setList[[4]],setList[[1]],myn=100,to.add=T,my.col="black",my.ylim=c(0,asizes[i]),my.title="Binding sites extra",to.smooth=F,use.para = FALSE)
makeRelPlot(setList[[5]],setList[[1]],myn=100,to.add=T,my.col="black",my.ylim=c(0,asizes[i]),my.title="Full set extra",to.smooth=F,use.para = FALSE)
dev.off()

library(doParallel)
library(GenomicRanges)
registerDoParallel(cores=30)

par(mfrow=c(2,2))
makeRelPlot(motGR$RRACH,matchGR$iclip,myn=100,to.add=F,my.col="black",my.ylim=c(0,100),my.title="Binding sites",to.smooth=F,use.para = FALSE)
makeRelPlot(motGR$RRACH,matchGR$iclip$random,myn=100,to.add=T,my.col="blue",my.ylim=c(0,100),my.title="Full set",to.smooth=F,use.para = FALSE)


#################################################################################################
#################MOTIFS AROUND iCLIP AND m6A#####################################################
#################################################################################################



mycols <- brewer.pal(7,"Set1")


# pdf("enrichment_of_RRACH_around_iCLIP.pdf",width=3.5,height=4)
# par(mfrow=c(1,1))
# makeRelPlot(motGR$RRACH,iCLIP[[3]],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(15,120),my.title="Binding sites",to.smooth=F,use.para = TRUE)
# makeRelPlot(motGR$RRACH,matchGR$random[matchGR$iclip$match],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,50),my.title="Full set",to.smooth=F,use.para = TRUE)
# dev.off()


# yvals <- c("URUGUAY"=42,"WEI"=42,"UGUCUC"=30,"SCGKA"=8,"UUCCS"=6,"URCWC"=10,"AAGWC"=8,"UGACA"=5,"ACUCU"=10,"GGAUW"=20,"CUAUN"=20,"GMUAY"=18,"UGYAA"=18,"YUGUM"=60,"UUAKS"=35,"UHADG"=110,"URUAY"=60,"RRACH"=130,"AAUAA"=19,"AVAYU"=20)
# names(yvals) %in% names(motGR)
# yvals_add <- rep(50,length(names(motGR)[!(names(motGR) %in% yvals)]))
# names(yvals_add) <- names(motGR)[!(names(motGR) %in% yvals)]
# yvals <- c(yvals,yvals_add)
# 

names(motGR) %in% names(yvals)
for(i in 1:length(motGR))
{
  pdf(paste0("enrichment_of_",names(motGR)[i],"_around_iCLIP.pdf"),width=3.5,height=4)
  par(mfrow=c(1,1))
  makeRelPlot(motGR[[i]],iCLIP[[3]],myn=200,to.add=F,my.col=mycols[1],my.ylim=NULL,my.title="Binding sites",to.smooth=F,use.para = TRUE)
  #makeRelPlot(motGR[[i]],iCLIP[[3]],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,yvals[names(motGR)[i]]),my.title="Binding sites",to.smooth=F,use.para = TRUE)
  makeRelPlot(motGR[[i]],matchGR$random[matchGR$iclip$match],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,50),my.title="Full set",to.smooth=F,use.para = TRUE)
  dev.off()
}
# 
# yvals_nano <- c("WEI"=140,"SCGKA"=17,"UUCCS"=19,"URCWC"=22,"AAGWC"=75,"UGACA"=35,"ACUCU"=38,"GGAUW"=200,"CUAUN"=45,"GMUAY"=40,"UGYAA"=40,"YUGUM"=95,"UUAKS"=45,"UHADG"=180,"URUAY"=120,"RRACH"=400,"UUAKC"=20,"AAUAA"=30,"AVAYU"=300)
# names(yvals_nano) %in% names(motGR)
# yvals_add <- rep(50,length(names(motGR)[!(names(motGR) %in% names(yvals_nano))]))
# names(yvals_add) <- names(motGR)[!(names(motGR) %in% names(yvals_nano))]
# yvals_nano <- c(yvals_nano,yvals_add)

for(i in 1:length(motGR))
{
  pdf(paste0("enrichment_of_",names(motGR)[i],"_around_NANO_25bp.pdf"),width=3.5,height=4)
  par(mfrow=c(1,1))
  makeRelPlot(motGR[[i]],nanoGR,myn=25,to.add=F,my.col=mycols[1],my.ylim=NULL,my.title="Binding sites",to.smooth=F,use.para = TRUE)
  #makeRelPlot(motGR[[i]],nanoGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,yvals_nano[names(motGR)[i]]),my.title="Binding sites",to.smooth=F,use.para = TRUE)
  makeRelPlot(motGR[[i]],matchGR$random[matchGR$nano$match],myn=25,to.add=T,my.col=mycols[2],my.ylim=c(0,50),my.title="Full set",to.smooth=F,use.para = TRUE)
  dev.off()
}


#################GIANT###################

IDRGR <- iCLIP[[3]][countOverlaps(iCLIP[[3]]+10,iCLIP[[1]])==0]
NIDRGR <- iCLIP[[3]][countOverlaps(iCLIP[[3]]+10,iCLIP[[1]])>0]

st_idr <- table(as.vector(getSeq(Athaliana,IDRGR[as.vector(strand(IDRGR)) %in% c("-","+")])))
st_nidr <- table(as.vector(getSeq(Athaliana,NIDRGR[as.vector(strand(NIDRGR)) %in% c("-","+")])))
chisq.test(rbind(st_idr,st_nidr))
st_idr/sum(st_idr)
st_nidr/sum(st_nidr)

load(paste0("/binf-isilon/PBgrp/qbp693/ect2_protoplast/output_rnaseq_smartseq_analyses/pf2_de_output_smartseq.Rdata"))

down_genes <- resLFC_pf2_df_smartseq[resLFC_pf2_df_smartseq$log2FoldChange<(-0.5),]$gene_id
up_genes <- resLFC_pf2_df_smartseq[resLFC_pf2_df_smartseq$log2FoldChange>0.5,]$gene_id
nc_genes <- resLFC_pf2_df_smartseq[abs(resLFC_pf2_df_smartseq$log2FoldChange)<0.5,]$gene_id

nanoGR_i <- nanoGR[(countOverlaps(nanoGR+25,iCLIP[[3]])>0)]
nanoGR_ni <- nanoGR[(countOverlaps(nanoGR+25,iCLIP[[3]])==0)]

targ_genes <- genes.list$genes.union
targ_genes_strict <- genes.list$genes.intersect
non_targ_genes <- genes.list$genes.nontargets


mycols <- brewer.pal(7,"Set1")
#up genes, #ddown genes


names_ord <- names(ucontent)[order(ucontent)]
for(i in 1:length(motGR))
{
  
  mot_name <- names_ord[i]
  k <- which(names(motGR)==mot_name)
  print(mot_name)
  
  pdf(paste0("enrichment_of_motifs_around_iCLIP_nano_giant_",i,"_",mot_name,".pdf"),width=3.5,height=7*3.6)
  
  par(mfrow=c(7,1))
  
  ###m6A#############
  pSets <- list(nanoGR,matchGR$random[matchGR$nano$match])
  makeRelPlotSets(motGR[[k]],pSets,myn=200,my.col=mycols,my.title=paste0(mot_name,"_nanopore"),use.para=TRUE)
  
  #dev.off()
  ###m6A zoomed in########
  makeRelPlotSets(motGR[[k]],pSets,myn=25,my.col=mycols,my.title=paste0(mot_name,"_nanopore_zoom"),use.para=TRUE)
  
  ###ICLIP#############
  pSets <- list(iCLIP[[3]],matchGR$random[matchGR$iclip$match])
  makeRelPlotSets(motGR[[k]],pSets,myn=200,my.col=mycols,my.title=paste0(mot_name,"_iclip"),use.para=TRUE)
  
  ###ICLIP#############
  pSets <- list(nanoGR_i,nanoGR_ni)
  makeRelPlotSets(motGR[[k]],pSets,myn=200,my.col=mycols,my.title=paste0(mot_name,"_nano_near_clip"),use.para=TRUE)
  legend("topright",c("near","far"),col=mycols[1:2],lty=1,lwd=1)
  
  
  #####IDR###########
  pSets <- list(NIDRGR,IDRGR)
  makeRelPlotSets(motGR[[k]],pSets,myn=200,my.col=mycols,my.title=paste0(mot_name,"_IDR"),use.para=TRUE)
  legend("topright",c("idr","non-idr"),col=mycols[1:2],lty=1,lwd=1)
  
  
  ####UPDOWN##########
  pSets <- list(nanoGR[nanoGR$gene %in% up_genes],nanoGR[nanoGR$gene %in% down_genes],nanoGR[nanoGR$gene %in% nc_genes])
  makeRelPlotSets(motGR[[k]],pSets,myn=200,my.col=mycols,my.title=paste0(mot_name,"_nano_by_DE"),use.para=TRUE)
  legend("topright",c("up","down","nc"),col=mycols[1:3],lty=1,lwd=1)
  
  ###TARGETS########
  pSets <- list(nanoGR[nanoGR$gene %in% targ_genes],nanoGR[nanoGR$gene %in% targ_genes_strict],nanoGR[nanoGR$gene %in% non_targ_genes])
  makeRelPlotSets(motGR[[k]],pSets,myn=200,my.col=mycols,my.title=paste0(mot_name,"_nano_by_TARGs"),use.para=TRUE)
  legend("topright",c("targ","strict","non"),col=mycols[1:3],lty=1,lwd=1)
  
  dev.off()
}


#############50bp series of GIANT##################
start(nanoGR) <- start(nanoGR)-1

names_ord <- c("YYYYY","UUUUU","UNUNU","URURU","URUAY")
ymax <- c(350,350,350,40,40)
ymin <- c(0,0,0,5,5)
for(i in 1:5)
{
  
  mot_name <- names_ord[i]
  k <- which(names(motGR)==mot_name)
  print(mot_name)
  
  pdf(paste0("enrichment_of_motifs_around_iCLIP_nano_giant_75BP_selected_SERIES_",i,"_",mot_name,".pdf"),width=3.5,height=7*3.6)
  
  par(mfrow=c(7,1))
  
  ###m6A#############
  pSets <- list(nanoGR,matchGR$random[matchGR$nano$match])
  makeRelPlotSets(motGR[[k]],pSets,myn=75,my.ylim = c(ymin[i],ymax[i]),my.col=mycols,my.title=paste0(mot_name,"_nanopore"),use.para=TRUE)
  
  #dev.off()
  ###m6A zoomed in########
  makeRelPlotSets(motGR[[k]],pSets,myn=75,my.ylim = c(ymin[i],ymax[i]),my.col=mycols,my.title=paste0(mot_name,"_nanopore_zoom"),use.para=TRUE)
  
  ###ICLIP#############
  pSets <- list(iCLIP[[3]],matchGR$random[matchGR$iclip$match])
  makeRelPlotSets(motGR[[k]],pSets,myn=75,my.ylim = c(ymin[i],ymax[i]),my.col=mycols,my.title=paste0(mot_name,"_iclip"),use.para=TRUE)
  
  ###ICLIP#############
  pSets <- list(nanoGR_i,nanoGR_ni)
  makeRelPlotSets(motGR[[k]],pSets,myn=75,my.ylim = c(ymin[i],ymax[i]),my.col=mycols,my.title=paste0(mot_name,"_nano_near_clip"),use.para=TRUE)
  legend("topright",c("near","far"),col=mycols[1:2],lty=1,lwd=1)
  
  
  #####IDR###########
  pSets <- list(NIDRGR,IDRGR)
  makeRelPlotSets(motGR[[k]],pSets,myn=75,my.ylim = c(ymin[i],ymax[i]),my.col=mycols,my.title=paste0(mot_name,"_IDR"),use.para=TRUE)
  legend("topright",c("idr","non-idr"),col=mycols[1:2],lty=1,lwd=1)
  
  
  ####UPDOWN##########
  pSets <- list(nanoGR[nanoGR$gene %in% up_genes],nanoGR[nanoGR$gene %in% down_genes],nanoGR[nanoGR$gene %in% nc_genes])
  makeRelPlotSets(motGR[[k]],pSets,myn=75,my.ylim = c(ymin[i],ymax[i]),my.col=mycols,my.title=paste0(mot_name,"_nano_by_DE"),use.para=TRUE)
  legend("topright",c("up","down","nc"),col=mycols[1:3],lty=1,lwd=1)
  
  ###TARGETS########
  pSets <- list(nanoGR[nanoGR$gene %in% targ_genes],nanoGR[nanoGR$gene %in% targ_genes_strict],nanoGR[nanoGR$gene %in% non_targ_genes])
  makeRelPlotSets(motGR[[k]],pSets,myn=75,my.ylim = c(ymin[i],ymax[i]),my.col=mycols,my.title=paste0(mot_name,"_nano_by_TARGs"),use.para=TRUE)
  legend("topright",c("targ","strict","non"),col=mycols[1:3],lty=1,lwd=1)
  
  dev.off()
}





#############50bp series of GIANT##################
start(nanoGR) <- start(nanoGR)-1

names_ord <- names(ucontent)[order(ucontent)]
for(i in 1:length(motGR))
{
  
  mot_name <- names_ord[i]
  k <- which(names(motGR)==mot_name)
  print(mot_name)
  
  pdf(paste0("enrichment_of_motifs_around_iCLIP_nano_giant_50BP_SERIES_",i,"_",mot_name,".pdf"),width=3.5,height=7*3.6)
  
  par(mfrow=c(7,1))
  
  ###m6A#############
  pSets <- list(nanoGR,matchGR$random[matchGR$nano$match])
  makeRelPlotSets(motGR[[k]],pSets,myn=50,my.col=mycols,my.title=paste0(mot_name,"_nanopore"),use.para=TRUE)
  
  #dev.off()
  ###m6A zoomed in########
  makeRelPlotSets(motGR[[k]],pSets,myn=50,my.col=mycols,my.title=paste0(mot_name,"_nanopore_zoom"),use.para=TRUE)
  
  ###ICLIP#############
  pSets <- list(iCLIP[[3]],matchGR$random[matchGR$iclip$match])
  makeRelPlotSets(motGR[[k]],pSets,myn=50,my.col=mycols,my.title=paste0(mot_name,"_iclip"),use.para=TRUE)
  
  ###ICLIP#############
  pSets <- list(nanoGR_i,nanoGR_ni)
  makeRelPlotSets(motGR[[k]],pSets,myn=50,my.col=mycols,my.title=paste0(mot_name,"_nano_near_clip"),use.para=TRUE)
  legend("topright",c("near","far"),col=mycols[1:2],lty=1,lwd=1)
  
  
  #####IDR###########
  pSets <- list(NIDRGR,IDRGR)
  makeRelPlotSets(motGR[[k]],pSets,myn=50,my.col=mycols,my.title=paste0(mot_name,"_IDR"),use.para=TRUE)
  legend("topright",c("idr","non-idr"),col=mycols[1:2],lty=1,lwd=1)
  
  
  ####UPDOWN##########
  pSets <- list(nanoGR[nanoGR$gene %in% up_genes],nanoGR[nanoGR$gene %in% down_genes],nanoGR[nanoGR$gene %in% nc_genes])
  makeRelPlotSets(motGR[[k]],pSets,myn=50,my.col=mycols,my.title=paste0(mot_name,"_nano_by_DE"),use.para=TRUE)
  legend("topright",c("up","down","nc"),col=mycols[1:3],lty=1,lwd=1)
  
  ###TARGETS########
  pSets <- list(nanoGR[nanoGR$gene %in% targ_genes],nanoGR[nanoGR$gene %in% targ_genes_strict],nanoGR[nanoGR$gene %in% non_targ_genes])
  makeRelPlotSets(motGR[[k]],pSets,myn=50,my.col=mycols,my.title=paste0(mot_name,"_nano_by_TARGs"),use.para=TRUE)
  legend("topright",c("targ","strict","non"),col=mycols[1:3],lty=1,lwd=1)
  
  dev.off()
}




for(i in 1:length(motGR))
{
  #pdf(paste0("enrichment_of_",names(motGR)[i],"_around_NANO_25bp.pdf"),width=3.5,height=4)
  #par(mfrow=c(1,1))
   dev.off()
}

####################################################
#####IDR HITS####################################################
####################################################

mycols <- brewer.pal(7,"Set1")
#m6AGR <- unique(unlist(GRangesList(nanoGR,micGR)))
#IDRGR <- iCLIP[[3]][(countOverlaps(iCLIP[[3]]+10,m6AGR)==0)]
IDRGR <- iCLIP[[3]][countOverlaps(iCLIP[[3]]+10,iCLIP[[1]])==0]
NIDRGR <- iCLIP[[3]][countOverlaps(iCLIP[[3]]+10,iCLIP[[1]])>0]
length(IDRGR)
length(NIDRGR)
length(iCLIP[[1]])
length(iCLIP[[3]])
for(i in 1:length(motGR))
{
  pdf(paste0("enrichment_of_",names(motGR)[i],"_around_iCLIP_IDR.pdf"),width=3.5,height=4)
  par(mfrow=c(1,1))
  makeRelPlot(motGR[[i]],NIDRGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=NULL,my.title="Binding sites",to.smooth=F,use.para = TRUE)
  #makeRelPlot(motGR[[i]],iCLIP[[3]],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,yvals[names(motGR)[i]]),my.title="Binding sites",to.smooth=F,use.para = TRUE)
  makeRelPlot(motGR[[i]],IDRGR,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,50),my.title="Full set",to.smooth=F,use.para = TRUE)
  dev.off()
}

iclip_match <- matchGR$random[matchGR$iclip$match]
iclip_fg <- iCLIP[[3]][sample(1:length(iclip_match))]
i <- 2
odds_vec <- rep(0,length(motGR))
pval_vec <- rep(0,length(motGR))
names(pval_vec) <- names(motGR)
names(odds_vec) <- names(motGR)
for(i in 1:length(motGR))
{
  cfg <- (countOverlaps(iclip_fg+50,motGR[[i]]))
  cbg <- (countOverlaps(iclip_match+50,motGR[[i]]))
  sum(cfg)/sum(cbg)
  
  cfgo <- (countOverlaps(shift(iclip_fg,150)+75,motGR[[i]])) + (countOverlaps(shift(iclip_fg,-150)+75,motGR[[i]]))
  cbgo <-  (countOverlaps(shift(iclip_match,150)+75,motGR[[i]])) + (countOverlaps(shift(iclip_match,-150)+75,motGR[[i]]))
  sum(cfgo)/sum(cbgo)
  odds_vec[i] <-  (sum(cfg)/sum(cbg))/(sum(cfgo)/sum(cbgo))
  cbox <- rbind(c(sum(cfg),sum(cbg)),c(sum(cfgo),sum(cbgo)))
  pval_vec[i] <- fisher.test(cbox)$p.value
}
df <- cbind(pval_vec,odds_vec)
df[order(df[,2]),]


#################################################################################################
#################NANOPORE +/- RRACH AND GGAUW####################################################
#################################################################################################


pdf(paste0("enrichment_of_iCLIP_around_nano_RRACH_GGAUW.pdf"),width=4*3.5,height=4)
par(mfrow=c(1,4))
makeRelPlot(iCLIP[[3]],nanoGR[countOverlaps(nanoGR+5,motGR$RRACH)>0 & countOverlaps(nanoGR+5,motGR$GGAU)>0],myn=200,to.add=F,my.col=mycols[1],my.ylim=NULL,my.title="Binding sites",to.smooth=F,use.para = TRUE)
makeRelPlot(iCLIP[[3]],nanoGR[countOverlaps(nanoGR+5,motGR$RRACH)==0 & countOverlaps(nanoGR+5,motGR$GGAU)>0],myn=200,to.add=F,my.col=mycols[2],my.ylim=NULL,my.title="Binding sites",to.smooth=F,use.para = TRUE)
makeRelPlot(iCLIP[[3]],nanoGR[countOverlaps(nanoGR+5,motGR$RRACH)>0 & countOverlaps(nanoGR+5,motGR$GGAU)==0],myn=200,to.add=F,my.col=mycols[3],my.ylim=NULL,my.title="Binding sites",to.smooth=F,use.para = TRUE)
makeRelPlot(iCLIP[[3]],nanoGR[countOverlaps(nanoGR+5,motGR$RRACH)==0 & countOverlaps(nanoGR+5,motGR$GGAU)==0],myn=200,to.add=F,my.col=mycols[4],my.ylim=NULL,my.title="Binding sites",to.smooth=F,use.para = TRUE)
dev.off()


# pdf("enrichment_of_URUAY_around_iCLIP.pdf",width=3.5,height=4)
# par(mfrow=c(1,1))
# makeRelPlot(motGR$URUAY,iCLIP[[3]],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(8,65),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm = F)
# makeRelPlot(motGR$URUAY,matchGR$random[matchGR$iclip$match],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,135),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm = F)
# dev.off()
# 
# sum(countOverlaps(iCLIP[[3]]+50,motGR$URUAY))
# sum(countOverlaps(matchGR$random[matchGR$iclip$match]+50,motGR$URUAY))
# 
# 
# pdf("enrichment_of_RRACH_around_nano.pdf",width=3.5,height=4)
# par(mfrow=c(1,1))
# makeRelPlot(motGR$RRACH,nanoGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,750),my.title="Binding sites",to.smooth=F,use.para = TRUE)
# makeRelPlot(motGR$RRACH,matchGR$random[matchGR$nano$match],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,380),my.title="Full set",to.smooth=F,use.para = TRUE)
# dev.off()
# 
# pdf("enrichment_of_URUAY_around_nano.pdf",width=3.5,height=4)
# par(mfrow=c(1,1))
# makeRelPlot(motGR$URUAY,nanoGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(10,120),my.title="Binding sites",to.smooth=F,use.para = TRUE)
# makeRelPlot(motGR$URUAY,matchGR$random[matchGR$nano$match],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,60),my.title="Full set",to.smooth=F,use.para = TRUE)
# dev.off()

htl <- unique(unlist(GRangesList(pos.list.ECT2)))

# pdf("enrichment_of_sets_around_iCLIP.pdf",width=3.5,height=4)
# par(mfrow=c(1,1))
# makeRelPlot(nanoGR,iCLIP[[3]],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,.0115),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=T)
# makeRelPlot(micGR,iCLIP[[3]],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
# makeRelPlot(htl,iCLIP[[3]],myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
# makeRelPlot(weiGRC,iCLIP[[3]],myn=200,to.add=T,my.col=mycols[4],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
# legend("topright",c("nano","mic","ht","wei"),col=mycols[1:4],lty=1,lwd=1)
# dev.off()

# mycols <- brewer.pal(7,"Set1")
# pdf("enrichment_of_sets_around_m6A.pdf",width=3.5,height=4)
# par(mfrow=c(1,1))
# makeRelPlot(iCLIP[[3]],nanoGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,.013),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=T)
# makeRelPlot(micGR,nanoGR,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
# makeRelPlot(htl,nanoGR,myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
# makeRelPlot(weiGRC,nanoGR,myn=200,to.add=T,my.col=mycols[4],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
# legend("topright",c("iclip","mic","ht","wei"),col=mycols[1:4],lty=1,lwd=1)
# dev.off()

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_ht_around_iclip.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(htl,iCLIP[[3]],myn=750,to.add=F,my.col=mycols[1],my.ylim=c(0,55),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(htl,iCLIP$E2_110kWGR,myn=750,to.add=T,my.col=mycols[2],my.ylim=c(0,60),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
dev.off()

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_ht_around_iclip_all.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(htl,iCLIP[[3]],myn=750,to.add=F,my.col=mycols[1],my.ylim=c(0,60),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(htl,iCLIP$E2_110kWGR,myn=750,to.add=T,my.col=mycols[2],my.ylim=c(0,60),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(htl,iCLIP$E2_55kGR,myn=750,to.add=T,my.col=mycols[3],my.ylim=c(0,60),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(htl,iCLIP$E2_55kWGR,myn=750,to.add=T,my.col=mycols[4],my.ylim=c(0,60),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
dev.off()

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_uruay_around_iclip_all.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(motGR$URUAY,iCLIP$E2_110kGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,80),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
#makeRelPlot(motGR$URUAY,iCLIP$E2_110kWGR,myn=750,to.add=T,my.col=mycols[2],my.ylim=c(0,28),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$URUAY,iCLIP$E2_55kGR,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,28),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$URUAY,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,28),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
dev.off()



#################################################################################################
#################NANOPORE +/- RRACH##############################################################
#################################################################################################



nanoGRP <- nanoGR[countOverlaps(nanoGR+10,motGR$RRACH)>0]
nanoGRM <- nanoGR[countOverlaps(nanoGR+10,motGR$RRACH)==0]
dtn <- distanceToNearest(nanoGR,motGR$RRACH)
summary(as.data.frame(dtn)$distance)
length(nanoGRP)
length(nanoGRM)

ag <- c(unique(unlist(lapply(iCLIPBS,function(x) x$gene))),genes.list$genes.union,nanoGR$gene,weiGRC$gene,micGR$gene)
gtfE <- gtf[gtf$gene_id %in% ag]
#gtfE

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_around_iCLIP110_distances_55.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
iCLIP$E2_55kGR_not_110kGR <- iCLIP$E2_55kGR[countOverlaps(iCLIP$E2_55kGR+2,iCLIP$E2_110kGR)==0]
makeRelPlot(iCLIP$E2_55kGR_not_110kGR,iCLIP$E2_110kGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,0.02),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
iCLIP$E2_55kGR_not_110kGR <- iCLIP$E2_55kGR[countOverlaps(iCLIP$E2_55kGR+10,iCLIP$E2_110kGR)==0]
makeRelPlot(iCLIP$E2_55kGR_not_110kGR,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
iCLIP$E2_55kGR_not_110kGR <- iCLIP$E2_55kGR[countOverlaps(iCLIP$E2_55kGR+25,iCLIP$E2_110kGR)==0]
makeRelPlot(iCLIP$E2_55kGR_not_110kGR,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
iCLIP$E2_55kGR_not_110kGR <- iCLIP$E2_55kGR[countOverlaps(iCLIP$E2_55kGR+100,iCLIP$E2_110kGR)==0]
makeRelPlot(iCLIP$E2_55kGR_not_110kGR,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[4],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
legend("topleft",col=mycols[1:4],legend=c("2","10","25","100"))
dev.off()

iCLIP$E2_110kGR_not_55kGR <- iCLIP$E2_110kGR[countOverlaps(iCLIP$E2_110kGR+5,iCLIP$E2_55kGR)==0]
iCLIP$E2_55kGR_not_110kGR <- iCLIP$E2_55kGR[countOverlaps(iCLIP$E2_55kGR+5,iCLIP$E2_110kGR)==0]
iCLIP$E2_55kGR_not_110kGR_URUAY <- iCLIP$E2_55kGR_not_110kGR[countOverlaps(iCLIP$E2_55kGR_not_110kGR+20,motGR$URUAY)>0] 
iCLIP$E2_55kGR_not_110kGR_NO_URUAY <- iCLIP$E2_55kGR_not_110kGR[countOverlaps(iCLIP$E2_55kGR_not_110kGR+20,motGR$URUAY)==0] 

lapply(iCLIP,length)

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_around_iCLIP110_not55_URUAY_nn.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
#makeRelPlot(iCLIP$E2_55kGR,iCLIP$E2_110kGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,62),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_55kGR_not_110kGR_URUAY,iCLIP$E2_110kGR,myn=200,to.add=F,my.col=mycols[2],my.ylim=c(0,.018),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_55kGR_not_110kGR_NO_URUAY,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,.018),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_55kGR,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[4],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
dev.off()

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_around_iCLIP110_nn.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(iCLIP$E2_55kWGR,iCLIP$E2_110kGR,myn=200,to.add=F,my.col=mycols[3],my.ylim=c(0,0.025),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_110kWGR,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[4],my.ylim=c(0,110),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_55kGR,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[1],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_55kGR_not_110kGR,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
dev.off()

iCLIP$E2_110kGR_not_55kGR <- iCLIP$E2_110kGR[countOverlaps(iCLIP$E2_110kGR+5,iCLIP$E2_55kGR)==0]
summary(countOverlaps(iCLIP$E2_110kGR+5,iCLIP$E2_55kGR)==0)
iCLIP$E2_110kGR_URUAY <- iCLIP$E2_110kGR[countOverlaps(iCLIP$E2_110kGR+10,motGR$URUAY)>0]
iCLIP$E2_110kGR_NO_URUAY <- iCLIP$E2_110kGR[countOverlaps(iCLIP$E2_110kGR+10,motGR$URUAY)==0]
lapply(iCLIP,length)

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip110_around_iCLIP55_uruay.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(iCLIP$E2_110kGR,iCLIP$E2_55kGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,0.018),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_110kGR_URUAY,iCLIP$E2_55kGR,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,110),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_110kGR_NO_URUAY,iCLIP$E2_55kGR,myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,110),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
dev.off()


mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_around_iCLIP110_not55.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
#makeRelPlot(iCLIP$E2_55kGR,iCLIP$E2_110kGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,62),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_55kGR_not_110kGR,iCLIP$E2_110kGR,myn=200,to.add=F,my.col=mycols[2],my.ylim=c(0,16),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_55kWGR,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_110kWGR,iCLIP$E2_110kGR,myn=200,to.add=T,my.col=mycols[4],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
dev.off()


#################################################################################################
#################iCLIP sites in 55 but not 110KDa#################################################
#################################################################################################

iCLIP$E2_110kGR_not_55kGR <- iCLIP$E2_110kGR[countOverlaps(iCLIP$E2_110kGR+50,iCLIP$E2_55kGR)==0]
iCLIP$E2_55kGR_not_110kGR <- iCLIP$E2_55kGR[countOverlaps(iCLIP$E2_55kGR+5,iCLIP$E2_110kGR)==0]
lapply(iCLIP,length)

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_55_110_around_nanopore.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
#makeRelPlot(iCLIP$E2_55kGR,iCLIP$E2_110kGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,62),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_55kGR_not_110kGR,nanoGR,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,0.012),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_55kGR,nanoGR,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_110kGR,nanoGR,myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_110kGR_not_55kGR,nanoGR,myn=200,to.add=T,my.col=mycols[4],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
dev.off()


mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_NOT_55_around_rrach_pm_nano_all.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(iCLIP$E2_110kGR_not_55kGR,nanoGRP,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_110kGR_not_55kGR,nanoGRM,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,110),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
dev.off()


mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_around_rrach_pm_nano_all.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(iCLIP$E2_110kGR,nanoGRP,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,310),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_110kGR,nanoGRM,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,160),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
dev.off()



mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_110_55_around_rrach_pm_nano_all.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(iCLIP$E2_55kGR,nanoGRP,myn=200,to.add=F,my.col=mycols[3],my.ylim=c(0,300),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_55kGR,nanoGRM,myn=200,to.add=T,my.col=mycols[4],my.ylim=c(0,35),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_110kGR,nanoGRP,myn=200,to.add=T,my.col=mycols[1],my.ylim=c(0,180),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_110kGR,nanoGRM,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,160),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
dev.off()

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_110_55_around_rrach_pm_nano_all_NORM.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(iCLIP$E2_55kGR,nanoGRP,myn=200,to.add=F,my.col=mycols[3],my.ylim=c(0,0.012),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_55kGR,nanoGRM,myn=200,to.add=T,my.col=mycols[4],my.ylim=c(0,35),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_110kGR,nanoGRP,myn=200,to.add=T,my.col=mycols[1],my.ylim=c(0,180),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP$E2_110kGR,nanoGRM,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,160),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
dev.off()

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_55_around_rrach_pm_nano_all.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(iCLIP$E2_55kGR,nanoGRP,myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,80),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP$E2_55kGR,nanoGRM,myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,35),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
dev.off()

#################################################################################################
#################up, down and no changes#########################################
#################################################################################################



load(paste0("/binf-isilon/PBgrp/qbp693/ect2_protoplast/output_rnaseq_smartseq_analyses/pf2_de_output_smartseq.Rdata"))

down_genes <- resLFC_pf2_df_smartseq[resLFC_pf2_df_smartseq$log2FoldChange<(-0.5),]$gene_id
up_genes <- resLFC_pf2_df_smartseq[resLFC_pf2_df_smartseq$log2FoldChange>0.5,]$gene_id
nc_genes <- resLFC_pf2_df_smartseq[abs(resLFC_pf2_df_smartseq$log2FoldChange)<0.5,]$gene_id

#up genes, #ddown genes
mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_iclip_around_m6A_up_down_nc.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(iCLIP[[3]],nanoGR[nanoGR$gene %in% up_genes],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,300),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP[[3]],nanoGR[nanoGR$gene %in% down_genes],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,300),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(iCLIP[[3]],nanoGR[nanoGR$gene %in% nc_genes],myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,300),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
legend("topright",c("up","down","nc"),col=mycols[1:3],lty=1,lwd=1)
dev.off()


mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_wei_around_m6A_up_down_nc.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(motGR$WEI,nanoGR[nanoGR$gene %in% up_genes],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,120),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$WEI,nanoGR[nanoGR$gene %in% down_genes],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,120),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$WEI,nanoGR[nanoGR$gene %in% nc_genes],myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,120),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
legend("topright",c("up","down","nc"),col=mycols[1:3],lty=1,lwd=1)
dev.off()


mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_wei_around_iclip_up_down_nc.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(motGR$WEI,iCLIP[[3]][iCLIP[[3]]$gene %in% up_genes],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,70),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$WEI,iCLIP[[3]][iCLIP[[3]]$gene %in% down_genes],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$WEI,iCLIP[[3]][iCLIP[[3]]$gene %in% nc_genes],myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
legend("topright",c("up","down","nc"),col=mycols[1:3],lty=1,lwd=1)
dev.off()

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_rrach_around_iclip_up_down_nc.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(motGR$RRACH,iCLIP[[3]][iCLIP[[3]]$gene %in% up_genes],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,70),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$RRACH,iCLIP[[3]][iCLIP[[3]]$gene %in% down_genes],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$RRACH,iCLIP[[3]][iCLIP[[3]]$gene %in% nc_genes],myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
legend("topright",c("up","down","nc"),col=mycols[1:3],lty=1,lwd=1)
dev.off()

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_uuccs_around_iclip_up_down_nc.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(motGR$UUCCS,iCLIP[[3]][iCLIP[[3]]$gene %in% up_genes],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,25),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$UUCCS,iCLIP[[3]][iCLIP[[3]]$gene %in% down_genes],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,25),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$UUCCS,iCLIP[[3]][iCLIP[[3]]$gene %in% nc_genes],myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,25),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
legend("topright",c("up","down","nc"),col=mycols[1:3],lty=1,lwd=1)
dev.off()


#################################################################################################
#################targets and non-targets#########################################
#################################################################################################


targ_genes <- genes.list$genes.union
targ_genes_strict <- genes.list$genes.intersect
non_targ_genes <- genes.list$genes.nontargets


mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_wei_around_m6A_targ_nt.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(motGR$WEI,nanoGR[nanoGR$gene %in% targ_genes],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,150),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$WEI,nanoGR[nanoGR$gene %in% targ_genes_strict],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,120),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$WEI,nanoGR[nanoGR$gene %in% non_targ_genes],myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,120),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
legend("topright",c("targ","strict","non"),col=mycols[1:3],lty=1,lwd=1)
dev.off()

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_wei_around_iclip_targ_nt.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(motGR$WEI,iCLIP[[3]][iCLIP[[3]]$gene %in% targ_genes],myn=200,to.add=F,my.col=mycols[1],my.ylim=c(0,70),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$WEI,iCLIP[[3]][iCLIP[[3]]$gene %in% targ_genes_strict],myn=200,to.add=T,my.col=mycols[2],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(motGR$WEI,iCLIP[[3]][iCLIP[[3]]$gene %in% non_targ_genes],myn=200,to.add=T,my.col=mycols[3],my.ylim=c(0,70),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=F)
legend("topright",c("targ","strict","non"),col=mycols[1:3],lty=1,lwd=1)
dev.off()


#################################################################################################
#################m6A around iCLIP by score (martin/tino)#########################################
#################################################################################################


mycols <- brewer.pal(7,"Set1")
gtfE <- gtf

setGR <- list()
setGR[[1]] <- iCLIPBS[[3]]
setGR[[2]] <- iCLIP[[3]][iCLIP[[3]]$score>=1]
setGR[[3]] <- iCLIP[[3]][iCLIP[[3]]$score<1]
lapply(setGR,length)

unlist(lapply(setGR,length))
mean(countOverlaps(setGR[[1]],micGR))
mean(countOverlaps(setGR[[2]],micGR))
mean(countOverlaps(setGR[[3]],micGR))

pdf("miCLIP_around_iCLIP_sets.pdf",width=3.5,height=4)
makeRelPlot(micGR,setGR[[1]],myn=100,to.add=F,my.col=mycols[1],my.ylim=c(0,2000),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(micGR,setGR[[2]],myn=100,to.add=T,my.col=mycols[2],my.ylim=c(0,300),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(micGR,setGR[[3]],myn=100,to.add=T,my.col=mycols[3],my.ylim=c(0,70),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
legend("topright",col=mycols[1:3],legend=c("Binding sites (8804)","iCLIP peaks score>=1 (11275)","iCLIP peaks score<1 (4685)"))
dev.off()


pdf("nanopore_around_iCLIP_sets.pdf",width=3.5,height=4)
makeRelPlot(nanoGR,setGR[[1]],myn=100,to.add=F,my.col=mycols[1],my.ylim=c(0,500),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(nanoGR,setGR[[2]],myn=100,to.add=T,my.col=mycols[2],my.ylim=c(0,300),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
makeRelPlot(nanoGR,setGR[[3]],myn=100,to.add=T,my.col=mycols[3],my.ylim=c(0,70),my.title="Binding sites",to.smooth=F,use.para = TRUE,to.norm=F)
legend("topright",col=mycols[1:3],legend=c("Binding sites (8804)","iCLIP peaks score>=1 (11275)","iCLIP peaks score<1 (4685)"))
dev.off()

setGR <- list()
setGR[[1]] <- iCLIPBS[[3]]
setGR[[2]] <- iCLIPBS[[3]][iCLIPBS[[3]]$score>=1]
setGR[[3]] <- iCLIPBS[[3]][iCLIPBS[[3]]$score<1]

unlist(lapply(setGR,length))
mean(countOverlaps(setGR[[1]],micGR))
mean(countOverlaps(setGR[[2]],micGR))
mean(countOverlaps(setGR[[3]],micGR))



#################################################################################################
#################EXAMPLES OF iCLIP PEAKS WITH URUAY/m6A WITH RRACH ETC############################
#################################################################################################

names(nanoGR) <- paste(as.vector(seqnames(nanoGR)),end(nanoGR),sep="_")
start(nanoGR) <- end(nanoGR)
nanoGR

nanoGR_RRACH <- nanoGR[countOverlaps(nanoGR+3,motGR$RRACH)>0]
nanoGR_RRACH <- nanoGR[countOverlaps(nanoGR+3,motGR$GGAU)>0]
nanoGR_RRACH <- nanoGR[countOverlaps(nanoGR+3,motGR$ACUCU)>0]


iCLIP_URUAY <- iCLIP[[3]][countOverlaps( iCLIP[[3]]+3,motGR$URUAY)>0]  

hitsGR <- (iCLIP_URUAY[countOverlaps(iCLIP_URUAY+10, nanoGR_RRACH)>0])
hitsGR$near <- names(nanoGR_RRACH[nearest(hitsGR,nanoGR_RRACH)])

mydf <- data.frame(hitsGR[rev(order(hitsGR$score))]$gene,hitsGR[rev(order(hitsGR$score))]$near)
mydf <- cbind(mydf,ids[as.vector(mydf$hitsGR.rev.order.hitsGR.score....gene)])
head(mydf)
colnames(mydf) <- c("gene_id","m6A_loc","name")
head(mydf)
dim(mydf)
mydf <- unique(mydf)

write.table(mydf,"/binf-isilon/alab/people/sarah/example_genes_iCLIP_URUAY_m6A_RRACH.txt",quote=F,col.names = T,row.names = F)
write.table(mydf,"/binf-isilon/alab/people/sarah/example_genes_iCLIP_URUAY_m6A_GGAU.txt",quote=F,col.names = T,row.names = F)
write.table(mydf,"/binf-isilon/alab/people/sarah/example_genes_iCLIP_URUAY_m6A_ACUCU.txt",quote=F,col.names = T,row.names = F)

nanoGRT <- nanoGR[nanoGR$gene %in% genes.list$genes.union]
nanoGRT <- micGR

pvm <- rep(0,length(motGR))
for(i in 1:length(motGR))
{
  olsize <- 50
  uol <- countOverlaps(nanoGRT+olsize,motGR[[i]])
  #uol <- countOverlaps(nanoGRT+olsize,motGR$UUVUS)
  
  table(uol)
  uol[uol>4] <- 4
  iol <- countOverlaps(nanoGRT+olsize,iCLIP[[3]])
  table(iol)
  iol[iol>0] <- 1
  #boxplot(iol~uol,outline=F)
  barplot(tapply(iol,uol,mean))
  #t.test(iol[uol==0],iol[uol>2])
  
  pvm[i] <- (summary(glm(factor(iol)~uol,family="binomial")))$coefficients[2,4]
  
}
names(pvm) <- names(motGR)
rev(sort(pvm))


#olsize <- 25
nanoGRR <- matchGR$random[matchGR$nano$match]
nolu <- (countOverlaps(nanoGR+olsize,motGR$URUGUAY))
nrolu <- (countOverlaps(nanoGRR+olsize,motGR$URUGUAY))

nwin_U <- rep(0,100)
nwin_R <- rep(0,100)
for(i in 1:100)
{
  nwin_U[i] <- sum(countOverlaps(nanoGR+i-1,motGR$URUAY)>0)
  nwin_R[i] <- sum(countOverlaps(nanoGR+i-1,motGR$RRACH)>0)
}

plot(nwin_R,type="l",ylim=c(0,max(nwin_R)),ylab="Number of motif sites",xlab="Increasing window size (m6A +/- x)",lwd=2,col=brewer.pal(3,"Set2")[1],main="Number of motif sites at \n increasing intervals from m6A")
points(nwin_U,type="l",col=brewer.pal(3,"Set2")[2],lwd=2)
abline(v=5,lty=2)
nwin_U[4]
nwin_R[4]



nanoGR[nanoGR$score]

mn <- nanoGR[as.vector(strand(nanoGR)) %in% c("+","-")]
mytab <- table(as.vector(getSeq(Athaliana,shift(mn,+1))))
mytab/sum(mytab)

nanoGR_extra <- nanoGR[nanoGR$log_odds<(-1)]

mean(countOverlaps(nanoGR_filt+12,iCLIP[[3]]))
mean(countOverlaps(nanoGR_extra+12,iCLIP[[3]]))

mean(countOverlaps(nanoGR_filt+5,motGR$DRAC))
mean(countOverlaps(nanoGR_extra+5,motGR$DRAC))

mean(countOverlaps(nanoGR_filt+5,micGR))
mean(countOverlaps(nanoGR_extra+5,micGR))



as.data.frame(distanceToNearest(nanoGR_extra,nanoGR_filt))
(distanceToNearest(nanoGR_extra,nanoGR_filt))



hist(mn$score)


