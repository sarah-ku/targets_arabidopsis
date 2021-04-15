##########################################################################
#####RANDOM FOREST########################################################
##########################################################################
source("/home/sarah/github/targets_arabidopsis/targets_paper_open.R")
source("/home/sarah/github/targets_arabidopsis/motifs_functions.R")
setwd("/binf-isilon/alab/people/sarah/meme/motif_plots")

library(gbm)
require(pROC)
library(ggplot2)
#library(pheatmap)
library(rtracklayer)
library(GenomicRanges)
library(RColorBrewer)
library("BSgenome.Athaliana.TAIR.TAIR9")
seqnames(Athaliana) <- c("1","2","3","4","5","Mt","Pt")

load(file="/binf-isilon/alab/people/sarah/meme/matched_sets.Rdat")
#load("/binf-isilon/alab/people/sarah/meme/matched_sets_m6A_genes.Rdat")
load(file="/binf-isilon/alab/people/sarah/meme/motGR.Rdat")
load(file="/binf-isilon/alab/people/sarah/meme/randGR.Rdat")

#gtfE <- gtf

nanoGRme <- createMerged(nanoGR)
matches <- matchGR$nano[subjectHits(findOverlaps(nanoGRme,matchGR$nano))]$match
matchGRme <- matchGR$random[matches]
nanoMatchGR <- sort(matchGRme)
RFdata_nano <- makeRFdata(nanoGRme = nanoGRme,nanoMatchGR = nanoMatchGR,add_size=25,center_size=10)

nanoGRme <- createMerged(nanoGR)
nanoGRme_i <- nanoGRme[(countOverlaps(nanoGRme+25,iCLIP[[3]])>0)]
nanoGRme_ni <- nanoGRme[(countOverlaps(nanoGRme+25,iCLIP[[3]])==0)]
RFdata_nano_iCLIP <- makeRFdata(nanoGRme = nanoGRme_i,nanoMatchGR = nanoGRme_ni,add_size=25,center_size=10)


iCLIPGRme <- createMerged(iCLIP[[3]])
matches <- matchGR$iclip[subjectHits(findOverlaps(iCLIPGRme,matchGR$iclip))]$match
matchGRme <- matchGR$random[matches]
#matchGRme <- createMerged(matchGR$random[matchGR$iclip$match])
iclipMatchGR <- sort(matchGRme)
RFdata_iCLIP <- makeRFdata(nanoGRme = iCLIPGRme,nanoMatchGR = iclipMatchGR)


#####################################################################
################## m6A vs background ################################
#####################################################################

pdf("RF_results_nanopore_background_10bp_stages.pdf",width=10,height=24)
par(mfrow=c(5,2))

#####VERSION WITH ALL
resp <- RFdata_nano$resp
X <- RFdata_nano$X

boost_nano <- makeGRF(X=X,resp=resp)

dat_test <- getRFtest(resp,X)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_nano,X_test)

#####RESTRICT to 10
bres <- summary(boost_nano,plot=F)
top10 <- bres[1:10,]$var
X_10 <- X[,top10]

boost_nano_10 <- makeGRF(X=X_10,resp=resp)

dat_test <- getRFtest(resp,X_10)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_nano_10,X_test)

#####RESTRICT to UNUNU + DRACH
topTerms <- grep("UNUNU|DRACH",colnames(X))
X_UD <- X[,topTerms]

boost_nano_UD <- makeGRF(X=X_UD,resp=resp)

dat_test <- getRFtest(resp,X_UD)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_nano_UD,X_test)

#####RESTRICT to UNUNU 
topTerms <- grep("UNUNU",colnames(X))
X_U <- X[,topTerms]

boost_nano_U <- makeGRF(X=X_U,resp=resp)

dat_test <- getRFtest(resp,X_U)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_nano_U,X_test)

#####RESTRICT to  DRACH
topTerms <- grep("DRACH",colnames(X))
X_D <- X[,topTerms]

boost_nano_D <- makeGRF(X=X_D,resp=resp)

dat_test <- getRFtest(resp,X_D)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_nano_D,X_test)

dev.off()


#####################################################################
################## iCLIP vs background ##############################
#####################################################################

pdf("RF_results_iCLIP_background_10bp_stages.pdf",width=10,height=24)
par(mfrow=c(5,2))


#####VERSION WITH ALL
resp <- RFdata_iCLIP$resp
X <- RFdata_iCLIP$X

boost_iCLIP <- makeGRF(X=X,resp=resp)

dat_test <- getRFtest(resp,X)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP,X_test)

#####RESTRICT to 10
bres <- summary(boost_iCLIP,plot=F)
top10 <- bres[1:10,]$var
X_10 <- X[,top10]

boost_iCLIP_10 <- makeGRF(X=X_10,resp=resp)

dat_test <- getRFtest(resp,X_10)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP_10,X_test)

#####RESTRICT to UNUNU + DRACH
topTerms <- grep("UNUNU|DRACH",colnames(X))
X_UD <- X[,topTerms]

boost_iCLIP_UD <- makeGRF(X=X_UD,resp=resp)

dat_test <- getRFtest(resp,X_UD)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP_UD,X_test)

#####RESTRICT to UNUNU 
topTerms <- grep("UNUNU",colnames(X))
X_U <- X[,topTerms]

boost_iCLIP_U <- makeGRF(X=X_U,resp=resp)

dat_test <- getRFtest(resp,X_U)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP_U,X_test)

#####RESTRICT to  DRACH
topTerms <- grep("DRACH",colnames(X))
X_D <- X[,topTerms]

boost_iCLIP_D <- makeGRF(X=X_D,resp=resp)

dat_test <- getRFtest(resp,X_D)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP_D,X_test)

dev.off()


#####################################################################
################## m6A +/- iCLIP ####################################
#####################################################################

pdf("RF_results_iCLIP_m6A_10bp_stages.pdf",width=10,height=24)
par(mfrow=c(5,2))


#####VERSION WITH ALL
resp <- RFdata_nano_iCLIP$resp
X <- RFdata_nano_iCLIP$X

boost_iCLIP <- makeGRF(X=X,resp=resp)

dat_test <- getRFtest(resp,X)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP,X_test)

#####RESTRICT to 10
bres <- summary(boost_iCLIP,plot=F)
top10 <- bres[1:10,]$var
X_10 <- X[,top10]

boost_iCLIP_10 <- makeGRF(X=X_10,resp=resp)

dat_test <- getRFtest(resp,X_10)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP_10,X_test)

#####RESTRICT to UNUNU + DRACH
topTerms <- grep("UNUNU|DRACH",colnames(X))
X_UD <- X[,topTerms]

boost_iCLIP_UD <- makeGRF(X=X_UD,resp=resp)

dat_test <- getRFtest(resp,X_UD)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP_UD,X_test)

#####RESTRICT to UNUNU 
topTerms <- grep("UNUNU",colnames(X))
X_U <- X[,topTerms]

boost_iCLIP_U <- makeGRF(X=X_U,resp=resp)

dat_test <- getRFtest(resp,X_U)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP_U,X_test)

#####RESTRICT to  DRACH
topTerms <- grep("DRACH",colnames(X))
X_D <- X[,topTerms]

boost_iCLIP_D <- makeGRF(X=X_D,resp=resp)

dat_test <- getRFtest(resp,X_D)
resp_test <- dat_test$resp_test
X_test <- dat_test$X
makeSummaryPlot(boost_iCLIP_D,X_test)

dev.off()