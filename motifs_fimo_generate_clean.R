library(BSgenome.Athaliana.TAIR.TAIR9)
source("/home/sarah/github/targets_arabidopsis/targets_paper_open.R")
source("/home/sarah/github/targets_arabidopsis/motifs_functions.R")
#load(file="/binf-isilon/alab/projects/ECT2_TC/FA_CLIP/weiGRC_singlebp.Rdat")
setwd("/binf-isilon/alab/people/sarah/meme/motif_plots")
library(ggplot2)
library(pheatmap)
library(rtracklayer)
library(GenomicRanges)
library(RColorBrewer)
library("BSgenome.Athaliana.TAIR.TAIR9")
seqnames(Athaliana) <- c("1","2","3","4","5","Mt","Pt")


load("/binf-isilon/alab/projects/ECT2_TC/binding sites unfiltered/iCLIPBS_with_meta.Rdat")
##check that the genes match properly for binding site and peak set
summary(unique(iCLIPBS$E2_110kGR$gene) %in% unique(iCLIP$E2_110kGR$gene))
summary(unique(iCLIP$E2_110kGR$gene) %in% unique(iCLIPBS$E2_110kGR$gene))

#################################################
###########DEFINE MOTIFS#########################
#################################################

#DRACH (D=G/A/U, R=G/A, H=A/U/C) 
# R = G/A, H = A/C/U, D=A/G/U, V=C/A/G, S=C/G, Y=C/T, W=T/A, K=T/G, M=C/A
R <- c("G","A")
A <- "A"
C <- "C"
H <- c("A","C","T")
U <- "T"
D <- c("A","G","T")
G <- "G"
V <- c("C","A","G")
Y <- c("C","T")
W <- c("T","A")
K <- c("T","G")
S <- c("C","G")
N <- c("C","T","A","G")
M <- c("C","A")

motif.list <- list()

UUGAA <- expand.grid(U,U,G,A,A)
motif.list[["UUGAA"]] <- UUGAA

UGGAU <- expand.grid(U,G,G,A,U)
motif.list[["UGGAU"]] <- UGGAU

UGUCUC <- expand.grid(U,G,U,C,U,C)
motif.list[["UGUCUC"]] <- UGUCUC

UUVUS <- expand.grid(U,U,V,U,S)
motif.list[["UUVUS"]] <- UUVUS

UGASA <- expand.grid(U,G,A,S,A)
motif.list[["UGASA"]] <- UGASA

UUGUG <- expand.grid(U,U,G,U,G)
motif.list[["UUGUG"]] <- UUGUG

UGAAC <- expand.grid(U,G,A,A,C)
motif.list[["UGAAC"]] <- UGAAC

UACUS <- expand.grid(U,A,C,U,S)
motif.list[["UACUS"]] <- UACUS

UGUKC <- expand.grid(U,G,U,K,C)
motif.list[["UGUKC"]] <- UGUKC

CUVUV <- expand.grid(C,U,V,U,V)
motif.list[["CUVUV"]] <- CUVUV

RRACH <- expand.grid(R,R,A,C,H)
motif.list[["RRACH"]] <- RRACH

UHADG <- expand.grid(U,H,A,D,G)
motif.list[["UHADG"]] <- UHADG

AVAYU <- expand.grid(A,V,A,Y,U)
motif.list[["AVAYU"]] <- AVAYU

GGAUW <- expand.grid(G,G,A,U,W)
motif.list[["GGAUW"]] <- GGAUW

GGAU <- expand.grid(G,G,A,U)
motif.list[["GGAU"]] <- GGAU

UUAKS <- expand.grid(U,U,A,K,S)
motif.list[["UUAKS"]] <- UUAKS

UUCCS <- expand.grid(U,U,C,C,S) 
motif.list[["UUCCS"]] <- UUCCS

URUAY <- expand.grid(U,R,U,A,Y)
motif.list[["URUAY"]] <- URUAY

WEI <- expand.grid(U,N,U,R,U,R,W)
motif.list[["WEI"]] <- WEI

URUGUAY <- expand.grid(U,R,U,G,U,A,Y)
motif.list[["URUGUAY"]] <- URUGUAY

AAUAA <- expand.grid(A,A,U,A,A)
motif.list[["AAUAA"]] <- AAUAA

ACUCU <- expand.grid(A,C,U,C,U)
motif.list[["ACUCU"]] <- ACUCU

UGACA <- expand.grid(U,G,A,C,A)
motif.list[["UGACA"]] <- UGACA

AAGWC <- expand.grid(A,A,G,W,C)
motif.list[["AAGWC"]] <- AAGWC

SCGKA <- expand.grid(S,C,G,K,A)
motif.list[["SCGKA"]] <- SCGKA

YUGUM <- expand.grid(Y,U,G,U,M)
motif.list[["YUGUM"]] <- YUGUM

UGYAA <- expand.grid(U,G,Y,A,A)
motif.list[["UGYAA"]] <- UGYAA

GMUAY <- expand.grid(G,M,U,A,Y)
motif.list[["GMUAY"]] <- GMUAY

CUAUN <- expand.grid(C,U,A,U,N)
motif.list[["CUAUN"]] <- CUAUN

URCWC <- expand.grid(U,R,C,W,C)
motif.list[["URCWC"]] <- URCWC

GACUU <- expand.grid(G,A,C,U,U)
motif.list[["GACUU"]] <- GACUU

#U rich placeholders
UNUNU <- expand.grid(U,N,U,N,U)
motif.list[["UNUNU"]] <- UNUNU

UUNUU <- expand.grid(U,U,N,U,U)
motif.list[["UUNUU"]] <- UUNUU

UNUUU <- expand.grid(U,N,U,U,U)
motif.list[["UNUUU"]] <- UNUUU

UUUUU <- expand.grid(U,U,U,U,U)
motif.list[["UUUUU"]] <- UUUUU

#-----------------
DRAY <- expand.grid(D,R,A,Y)
motif.list[["DRAY"]] <- DRAY

YYYYY <- expand.grid(Y,Y,Y,Y,Y)
motif.list[["YYYYY"]] <- YYYYY

URURU <- expand.grid(U,R,U,R,U)
motif.list[["URURU"]] <- URURU

ACUCUCU <- expand.grid(A,C,U,C,U,C,U)
motif.list[["ACUCUCU"]] <- ACUCUCU

AUUUUU <- expand.grid(A,U,U,U,U,U)
motif.list[["AUUUUU"]] <- AUUUUU

ACUUCUU <- expand.grid(A,C,U,U,C,U,U)
motif.list[["ACUUCUU"]] <- ACUUCUU

DRACUCU <- expand.grid(D,R,A,C,U,C,U)
motif.list[["DRACUCU"]] <- DRACUCU

UUUUUUU <- expand.grid(U,U,U,U,U,U,U)
motif.list[["UUUUUUU"]] <- UUUUUUU

UUUUUU <- expand.grid(U,U,U,U,U,U)
motif.list[["UUUUUU"]] <- UUUUUU

DRACG <- expand.grid(D,R,A,C,G)
motif.list[["DRACG"]] <- DRACG

URACH <- expand.grid(U,R,A,C,H)
motif.list[["URACH"]] <- URACH

DRACH <- expand.grid(D,R,A,C,H)
motif.list[["DRACH"]] <- DRACH

RRAC <- expand.grid(R,R,A,C)
motif.list[["RRAC"]] <- RRAC

motif.list <- lapply(motif.list,function(x) apply(x,1,function(y) paste0(y,collapse="")))
motif.list <- lapply(motif.list,function(x) paste0(x,collapse="|"))

mot.lets <- lapply(motif.list,function(x) strsplit(paste(x,collapse=""),"")[[1]])

#sort by U content for plotting order later
ucontent <- unlist(lapply(mot.lets,function(x) length(x[x=="T"])/length(x)))
save(ucontent,file="/binf-isilon/alab/people/sarah/meme/ucontent.Rdat")

#################################################
###########FIND MOTIFS AT BINDING SITES##########
#################################################

rrfm_freq_list <- list()
motif_c_list <- list()
motif <- motif.list[[1]]
size <- c(10)
for(motif in motif.list)
{
  rrfm_freq_list[[motif]] <- list()
  motif_c_list[[motif]] <- list()
  for(size in c(10,25,50,100))
  {
    #regsGR <- iclipBS[[1]][iclipBS[[1]]$type=="dominantPeak"]+size
    regsGR <- iCLIP[[3]]+size
    regsGR <- regsGR[-which(as.vector(strand(regsGR))=="*")]
    oor <- which(end(ranges(regsGR)) > seqlengths(Athaliana)[as.vector(seqnames(regsGR))])
    if(length(oor)>0){
      regsGR <- regsGR[-oor]
    }
    oor <- which(start(ranges(regsGR)) < 1)
    if(length(oor)>0)
    {
      regsGR <- regsGR[-oor]
    }
    
    my.mat <- as.matrix(getSeq(Athaliana, regsGR))
    
    B <- table(as.vector(my.mat))[c("A","T","G","C")]
    B <- B/sum(B)
    seqvec <- apply(my.mat,1,function(x) paste(x,collapse=""))
    
    rrach_sep <- strsplit(motif,"\\|")[[1]]
    rrach_counts <- rep(0,length(rrach_sep))
    names(rrach_counts) <- rrach_sep
    for(i in 1:length(rrach_sep)){
      rrach_counts[i] <- (length(grep(rrach_sep[i],seqvec,value=F)))
    }
    
    motif_c_list[[motif]][[paste(size)]] <- ((rrach_counts))
    
    mlen <- length(strsplit(rrach_sep[1],"")[[1]])
    rrfm <- matrix(0,ncol=4,nrow=mlen)
    colnames(rrfm) <- c("A","T","G","C")
    drc <- do.call(rbind,strsplit(rrach_sep,""))
    for(i in 1:ncol(drc)){
      rrfm[i,"A"] <- sum(rrach_counts[which(drc[,i]=="A")])
      rrfm[i,"T"] <- sum(rrach_counts[which(drc[,i]=="T")])
      rrfm[i,"G"] <- sum(rrach_counts[which(drc[,i]=="G")])
      rrfm[i,"C"] <- sum(rrach_counts[which(drc[,i]=="C")])
    }
    
    PPM <- function(C) C / sum(C)
    PPMp <- function(C,n) (C + n/length(C) )/ (sum(C) + n)
    mut <- round(mean(rrfm*0.1))
    #ppm <- apply(rrfm,1,function(x) (PPMp(x,mut)/B))
    
    ppm <- apply(rrfm,1,function(x) log2(PPMp(x,mut)/B))
    fmat <- apply(rrfm,1,function(x) PPMp(x,mut))
    colnames(ppm) <- paste(1:mlen)
    
    rrfm_freq_list[[motif]][[paste(size)]] <- list("B"=B,"fmat"=fmat)
  }
}


#################################################
###########LOGOS#################################
#################################################

require(ggplot2)
require(ggseqlogo)

almat <- lapply(rrfm_freq_list,function(x) x[["50"]]$fmat)
names(almat) <- names(motif.list)

makeLogo <- function(amat)
{
  x <- amat
  row.names(x)[2] <- "U"
  a <- ggseqlogo( x )
  return(a)
}

logos <- lapply(almat,function(x) makeLogo(x))


for(i in 1:length(logos))
{
  logos[[i]] <- logos[[i]] + ggtitle(names(almat)[i])
}

library(gridExtra)
sg <- grid.arrange(grobs=logos)
sg
ggsave(plot=sg,filename = "PWM_fimo_all.pdf",height=20,width=20)


#################################################
###########FIMO FILES############################
#################################################

names(rrfm_freq_list) <- names(motif.list)
motif.lengths <- unlist(lapply(rrfm_freq_list,function(x) ncol(x[["50"]]$fmat)))
fmats <- lapply(rrfm_freq_list,function(x) x[["50"]]$fmat)

system("rm /binf-isilon/alab/people/sarah/meme/motif_files/motifs_all_homer.meme")
system("touch /binf-isilon/alab/people/sarah/meme/motif_files/motifs_all_homer.meme")
topt <- c("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.273 C 0.165 G 0.173 T 0.389")
write.table(topt,file="/binf-isilon/alab/people/sarah/meme/motif_files/tmp.txt",quote=F,col.names=F,row.names=F)
system(paste0("cat /binf-isilon/alab/people/sarah/meme/motif_files/tmp.txt >> /binf-isilon/alab/people/sarah/meme/motif_files/motifs_all_homer.meme"))

for(i in 1:length(rrfm_freq_list))
{
  headline <- paste0("\nMOTIF ",names(motif.list)[i],"\nletter-probability matrix: alength= 4 w= ",motif.lengths[i])
  write.table(headline,file="/binf-isilon/alab/people/sarah/meme/motif_files/tmp.txt",quote=F,col.names=F,row.names=F)
  system(paste0("cat /binf-isilon/alab/people/sarah/meme/motif_files/tmp.txt >> /binf-isilon/alab/people/sarah/meme/motif_files/motifs_all_homer.meme"))
  write.table(t(fmats[[i]][c("A","C","G","T"),]),file="/binf-isilon/alab/people/sarah/meme/motif_files/tmp.txt",quote=F,col.names=F,row.names=F)
  system(paste0("cat /binf-isilon/alab/people/sarah/meme/motif_files/tmp.txt >> /binf-isilon/alab/people/sarah/meme/motif_files/motifs_all_homer.meme"))
}

system("head /binf-isilon/alab/people/sarah/meme/motif_files/motifs_all_homer.meme")
system("tail /binf-isilon/alab/people/sarah/meme/motif_files/motifs_all_homer.meme")

#################################################################
###GENOMEWIDE####################################################
#################################################################

##run in bash
#export PATH=/binf-isilon/alab/people/sarah/meme/bin:/binf-isilon/alab/people/sarah/meme/libexec/meme-5.1.1:$PATH

## run from /binf-isilon/alab/people/sarah/meme/results/all_from_homer
#fimo --verbosity 4 --thresh 0.05 --max-stored-scores 1000000000 ../../motif_files/motifs_all_homer.meme /binf-isilon/alab/projects/ECT2_TC/annotations/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa

################################################################
#######################OPEN AND PROCESS#########################
#################################################################



####genomewide##################################################
#all 48 in one file
twrr <- read.table("/binf-isilon/alab/people/sarah/meme/results/fimo_all/fimo_out/fimo.tsv",sep="\t",header=T)

cmots <- as.vector(unique(twrr$motif_id))
cmots

motGR <- list()
for(mot in cmots)
{
  twrr_4 <- twrr[(twrr$score>4) & twrr$motif_id==mot,]
  mGR <- GRanges(seqnames=Rle(twrr_4$sequence_name),strand=Rle(twrr_4$strand),IRanges(twrr_4$start+2,width=1),match_sequence=twrr_4$matched_sequence)
  motGR[[mot]] <- sort(mGR)
}
save(motGR,file="/binf-isilon/alab/people/sarah/meme/motGR.Rdat")


for(mot in cmots)
{
  forbed <- twrr[twrr$motif_id==mot,][,c("sequence_name","start","stop","matched_sequence","score","strand")]
  #forbed$start <- as.numeric(forbed$start)-1
  #forbed[!forbed$start==0,]$start <- 1
  forbed <- forbed[forbed$score>4,]
  forbed <- forbed[order(forbed$strand),]
  forbed <- forbed[order(forbed$sequence_name),]
  #colnames(forbed) <- c("chr","strand","start","end","score")
  write.table(forbed,file=paste0("/binf-isilon/alab/projects/ECT2_TC/motif_bed/motif_bed_",mot,".bed"),quote=FALSE,col.names = F,row.names=F,sep="\t")
}

# slist <- list()
# for(mot in cmots)
# {
#   slist[[mot]] <- twrr$score[twrr$motif_id==mot]
# }

#####random background nucleotides in transcriptome (make it expressed genes)
#matched.set.list$union$non_targets
gtf_expr <- unique(gtfGR[gtfGR$gene_id %in% genes.list$genes.expressed])
gtf_expr_df <- as.data.frame(gtf_expr)
min_start <- tapply(as.vector(gtf_expr_df$start),gtf_expr_df$gene_id,min)
max_end <- tapply(as.vector(gtf_expr_df$end),gtf_expr_df$gene_id,max)
my_strand <- tapply(gtf_expr_df$strand,gtf_expr_df$gene_id,function(x) as.vector(x)[1])
my_chr <- tapply(gtf_expr_df$seqnames,gtf_expr_df$gene_id,function(x) as.vector(x)[1])
expr_newdf <- data.frame("Chrom"=my_chr,"Strand"=my_strand,"Start"=min_start,"End"=max_end,"id"=names(max_end))
all.positions <- apply(expr_newdf,1,function(x) c(x[3]:x[4]))
all.positions <- unlist(all.positions)
all.seqnames <- apply(expr_newdf,1,function(x) rep(x[1],length(x[3]:x[4])))
all.seqnames <- unlist(all.seqnames)
all.strands <- apply(expr_newdf,1,function(x) rep(x[2],length(x[3]:x[4])))
all.strands <- unlist(all.strands)

pGR <- makeGRangesFromDataFrame(data.frame("seqnames"=all.seqnames,"Strand"=all.strands,"Start"=all.positions,"End"=all.positions))
randGR <- pGR

save(randGR, file = "/binf-isilon/alab/people/sarah/meme/randGR.Rdat")

#############################################
########MOTIFS GENOME WIDE###################
#############################################

for(mot in names(motGR))
{
  twGR <- motGR[[mot]]
  dr <- returnDens(twGR,gtfGR)
  dr_rand <- returnDens(randGR,gtfGR)
  pdf(paste0("normalised_enrichments_plot_",mot,".pdf"),width=6.5,height=6)
  plotNE(bs=0.02,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=T,to.add=F,col.no=1,no.samps=25,to.smooth=F,mytitle = mot)
  dev.off()
}

#############################################
########MATCH SETS BY GENOMIC POS############
#############################################

rmg <- unique(unlist(lapply(pos.list.ECT2,function(x) x$gene)))
#nontargs <- Reduce(union,lapply(matched.set.list,function(x) x$non_targets))
#nontargs <- matched.set.list$intersect$non_targets
nontargs <- genes.list$genes.nontargets
summary(nontargs %in% rmg)
nontargs <- nontargs[!(nontargs %in% rmg)]

#nontargs <- genes.list$genes.expressed

rGR <- randGR
rGR <- rGR[which(countOverlaps(rGR,gtfGR[gtfGR$gene_id %in% nontargs])>0)]
names(rGR) <- paste0(seqnames(rGR),":",start(rGR),",",strand(rGR))

htGR <- unique(unlist(GRangesList(pos.list.ECT2)))
ranges(weiGRC) <- IRanges(start = start(weiGRC)-1,end=end(weiGRC))
weiGRC <- unique(weiGRC)
setGR <- list("iclip"=iCLIPBS[[3]],"iclip55"=iCLIPBS$E2_55kGR,"iclip110W"=iCLIPBS$E2_110kWGR,"nano"=unique(nanoGR),"mic"=unique(micGR),"wei"=unique(weiGRC),"ht"=htGR)
for(i in 1:length(setGR)){ 
  ranges(setGR[[i]]) <- IRanges(end(unique(setGR[[i]])),width=1)
  names(setGR[[i]]) <- paste0(seqnames(setGR[[i]]),":",start(setGR[[i]]),",",strand(setGR[[i]])) 
}

matchGR <- getMatchedSets(tGR=setGR,rGR = rGR,mdb = 0.1)
nset <- length(matchGR)-1


save(matchGR,file="/binf-isilon/alab/people/sarah/meme/matched_sets.Rdat")
#save(matchGR,file="/binf-isilon/alab/people/sarah/meme/matched_sets_m6A_genes.Rdat")


###########################################
###########################################
###########################################


load(file="/binf-isilon/alab/people/sarah/meme/matched_sets.Rdat")
load(file="/binf-isilon/alab/people/sarah/meme/motGR.Rdat")


############################################
###########MAKE BED FILES FOR HOMER#########
############################################

fg <- matchGR$iclip
bg <- matchGR$random[matchGR$iclip$match]
bg <- unique(bg)

saveHomerMatch(fg+100,bg+100,my.name="iclip_100bp_homer")
saveHomerMatch(fg+50,bg+50,my.name="iclip_50bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="iclip_20bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="iclip_10bp_homer")
saveHomerMatch(fg+5,bg+5,my.name="iclip_5bp_homer")



fg <- matchGR$wei
bg <- matchGR$random[matchGR$wei$match]
bg <- unique(bg)

saveHomerMatch(fg+100,bg+100,my.name="wei_100bp_homer")
saveHomerMatch(fg+50,bg+50,my.name="wei_50bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="wei_20bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="wei_10bp_homer")
saveHomerMatch(fg+5,bg+5,my.name="wei_5bp_homer")


fg <- matchGR$nano
bg <- matchGR$random[matchGR$nano$match]
bg <- unique(bg)

saveHomerMatch(fg+100,bg+100,my.name="nano_100bp_homer")
saveHomerMatch(fg+50,bg+50,my.name="nano_50bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="nano_20bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="nano_10bp_homer")
saveHomerMatch(fg+5,bg+5,my.name="nano_5bp_homer")


fg <- matchGR$ht
bg <- matchGR$random[matchGR$ht$match]
bg <- unique(bg)

saveHomerMatch(fg+100,bg+100,my.name="ht_100bp_homer")
saveHomerMatch(fg+50,bg+50,my.name="ht_50bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="ht_20bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="ht_10bp_homer")
saveHomerMatch(fg+5,bg+5,my.name="ht_5bp_homer")

fg <- matchGR$mic
bg <- matchGR$random[matchGR$mic$match]
bg <- unique(bg)

saveHomerMatch(fg+100,bg+100,my.name="mic_100bp_homer")
saveHomerMatch(fg+50,bg+50,my.name="mic_50bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="mic_20bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="mic_10bp_homer")
saveHomerMatch(fg+5,bg+5,my.name="mic_5bp_homer")


fg <- iclipBS$E2mCh_110kDa_1.bsites.bed
bg <- iclipBS$E2mChW_110kDa_1.bsites.bed

saveHomerMatch(fg+100,bg+100,my.name="iclip_cage_110_100bp_homer")
saveHomerMatch(fg+50,bg+50,my.name="iclip_cage_110_50bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="iclip_cage_110_20bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="iclip_cage_110_10bp_homer")
saveHomerMatch(fg+5,bg+5,my.name="iclip_cage_110_5bp_homer")

fg <- iclipBS$E2mCh_55kDa_1.bsites.bed
bg <- iclipBS$E2mChW_55kDa_1.bsites.bed

saveHomerMatch(fg+100,bg+100,my.name="iclip_cage_55_100bp_homer")
saveHomerMatch(fg+50,bg+50,my.name="iclip_cage_55_50bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="iclip_cage_55_20bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="iclip_cage_55_10bp_homer")
saveHomerMatch(fg+5,bg+5,my.name="iclip_cage_55_5bp_homer")


fg <- iclipBS$E2mCh_110kDa_1.bsites.bed
bg <- iclipBS$E2mCh_55kDa_1.bsites.bed

saveHomerMatch(fg+5,bg+5,my.name="iclip_110_vs_55_5bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="iclip_110_vs_55_10bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="iclip_110_vs_55_20bp_homer")
saveHomerMatch(fg+50,bg+50,my.name="iclip_110_vs_55_50bp_homer")
saveHomerMatch(fg+100,bg+100,my.name="iclip_110_vs_55_100bp_homer")



saveHomerMatch <- function(fg,bg,my.name="iclip_20")
{
  #seqlevels(fg) <- gsub("Chr","",seqlevels(fg))
  my.bed <- as.data.frame(fg)[,c(1,2,3,5)]
  my.bed$id <- paste0(as.vector(seqnames(fg)),"_",start(fg))
  my.bed$nu <- "*"
  my.bed <- my.bed[,c(1:3,5,6,4)]
  head(my.bed)
  write.table(my.bed,file=paste0("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/motifs/results_new/",my.name,"_fg.bed"),quote=F,sep="\t",col.names=F,row.names=F)
  
  my.bed <- as.data.frame(bg)[,c(1,2,3,5)]
  my.bed$id <- paste0(as.vector(seqnames(bg)),"_",start(bg))
  my.bed$nu <- "*"
  my.bed <- my.bed[,c(1:3,5,6,4)]
  head(my.bed)
  write.table(my.bed,file=paste0("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/motifs/results_new/",my.name,"_bg.bed"),quote=F,sep="\t",col.names=F,row.names=F)
}


#PATH=$PATH:/binf-isilon/alab/people/sarah/homer/bin/
  
# findMotifsGenome.pl iclip_5bp_homer_fg.bed tair10 resultsENWSH_iclip_5bp_homer_fg -bg iclip_5bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_10bp_homer_fg.bed tair10 resultsENWSH_iclip_10bp_homer_fg -bg iclip_10bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_20bp_homer_fg.bed tair10 resultsENWSH_iclip_20bp_homer_fg -bg iclip_20bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_50bp_homer_fg.bed tair10 resultsENWSH_iclip_50bp_homer_fg -bg iclip_50bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_100bp_homer_fg.bed tair10 resultsENWSH_iclip_100bp_homer_fg -bg iclip_100bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
#   
# findMotifsGenome.pl wei_5bp_homer_fg.bed tair10 resultsENWSH_wei_5bp_homer_fg -bg wei_5bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl wei_10bp_homer_fg.bed tair10 resultsENWSH_wei_10bp_homer_fg -bg wei_10bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl wei_20bp_homer_fg.bed tair10 resultsENWSH_wei_20bp_homer_fg -bg wei_20bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl wei_50bp_homer_fg.bed tair10 resultsENWSH_wei_50bp_homer_fg -bg wei_50bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl wei_100bp_homer_fg.bed tair10 resultsENWSH_wei_100bp_homer_fg -bg wei_100bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
#   
# findMotifsGenome.pl nano_5bp_homer_fg.bed tair10 resultsENWSH_nano_5bp_homer_fg -bg nano_5bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl nano_10bp_homer_fg.bed tair10 resultsENWSH_nano_10bp_homer_fg -bg nano_10bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl nano_20bp_homer_fg.bed tair10 resultsENWSH_nano_20bp_homer_fg -bg nano_20bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl nano_50bp_homer_fg.bed tair10 resultsENWSH_nano_50bp_homer_fg -bg nano_50bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl nano_100bp_homer_fg.bed tair10 resultsENWSH_nano_100bp_homer_fg -bg nano_100bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
#  
# findMotifsGenome.pl ht_5bp_homer_fg.bed tair10 resultsENWSH_ht_5bp_homer_fg -bg ht_5bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl ht_10bp_homer_fg.bed tair10 resultsENWSH_ht_10bp_homer_fg -bg ht_10bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl ht_20bp_homer_fg.bed tair10 resultsENWSH_ht_20bp_homer_fg -bg ht_20bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl mic_5bp_homer_fg.bed tair10 resultsENWSH_mic_5bp_homer_fg -bg mic_5bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl mic_10bp_homer_fg.bed tair10 resultsENWSH_mic_10bp_homer_fg -bg mic_10bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl mic_20bp_homer_fg.bed tair10 resultsENWSH_mic_20bp_homer_fg -bg mic_20bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
#   
# findMotifsGenome.pl iclip_cage_55_5bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_55_5bp_homer_fg -bg iclip_cage_55_5bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_cage_55_10bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_55_10bp_homer_fg -bg iclip_cage_55_10bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_cage_55_20bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_55_20bp_homer_fg -bg iclip_cage_55_20bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_cage_55_50bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_55_50bp_homer_fg -bg iclip_cage_55_50bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_cage_55_100bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_55_100bp_homer_fg -bg iclip_cage_55_100bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
#   
# findMotifsGenome.pl iclip_cage_110_5bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_110_5bp_homer_fg -bg iclip_cage_110_5bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_cage_110_10bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_110_10bp_homer_fg -bg iclip_cage_110_10bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_cage_110_10bp_homer_fg.bed tair10 resultsESH_iclip_cage_110_10bp_homer_fg -bg iclip_cage_110_10bp_homer_bg.bed -size given -rna -h -len 5,6,7,8 &
# findMotifsGenome.pl iclip_cage_110_20bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_110_20bp_homer_fg -bg iclip_cage_110_20bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_cage_110_50bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_110_50bp_homer_fg -bg iclip_cage_110_50bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_cage_110_100bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_110_100bp_homer_fg -bg iclip_cage_110_100bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# 
# findMotifsGenome.pl iclip_110_vs_55_5bp_homer_fg.bed tair10 resultsENWSH_iclip_110_vs_55_5bp_homer_fg -bg iclip_110_vs_55_5bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_110_vs_55_10bp_homer_fg.bed tair10 resultsENWSH_iclip_110_vs_55_10bp_homer_fg -bg iclip_110_vs_55_10bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_110_vs_55_20bp_homer_fg.bed tair10 resultsENWSH_iclip_110_vs_55_20bp_homer_fg -bg iclip_110_vs_55_20bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_110_vs_55_50bp_homer_fg.bed tair10 resultsENWSH_iclip_110_vs_55_50bp_homer_fg -bg iclip_110_vs_55_50bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
# findMotifsGenome.pl iclip_110_vs_55_100bp_homer_fg.bed tair10 resultsENWSH_iclip_110_vs_55_100bp_homer_fg -bg iclip_110_vs_55_100bp_homer_bg.bed -size given -rna -h -nlen 1 -len 5,6,7,8 &
  