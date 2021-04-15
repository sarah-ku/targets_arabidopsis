source("/home/sarah/github/targets_arabidopsis/targets_paper_open.R")
source("/home/sarah/github/targets_arabidopsis/motifs_functions.R")
load(file="/binf-isilon/alab/projects/ECT2_TC/FA_CLIP/weiGRC_singlebp.Rdat")

setwd("/binf-isilon/alab/people/sarah/meme/motif_plots")


library(ggplot2)
library(pheatmap)
library(rtracklayer)
library(GenomicRanges)
library(RColorBrewer)
library("BSgenome.Athaliana.TAIR.TAIR9")
seqnames(Athaliana) <- c("1","2","3","4","5","Mt","Pt")

#################################################
###########DEFINE MOTIFS#########################
#################################################

#DRACH (D=G/A/U, R=G/A, H=A/U/C) 
# R = G/A, H = A/C/U, D=A/G/U, V=C/A/G, S=C/G, Y=C/T, W=T/A, K=T/G, M=C/A


RRACH <- expand.grid(c("G","A"),c("G","A"),c("A"),c("C"),c("A","T","C"))
UHADG <- expand.grid(c("T"),c("A","C","T"),"A",c("A","G","T"),"G")
AVAYU <- expand.grid(c("A"),c("C","A","G"),c("A"),c("C","T"),c("T"))
GGAUW <- expand.grid("G","G","A","T",c("T","A"))
UUAKS <- expand.grid("T","T","A",c("T","G"),c("G","C"))
UUCCS <- expand.grid("T","T","C","C",c("C","G")) 
URUAY <- expand.grid(c("T"),c("G","A"),c("T"),c("A"),c("T","C"))
WEI <- expand.grid(c("T"),c("G","T","A","C"),c("T"),c("G","A"),c("T"),c("A","G"),c("T","A"))

AAUAA <- expand.grid(c("A"),c("A"),c("T"),c("A"),c("A"))
ACUCU <- expand.grid(c("A"),c("C"),c("T"),c("C"),c("T"))
UGACA <- expand.grid(c("T"),c("G"),c("A"),c("C"),c("A"))
AAGWC <- expand.grid(c("A"),c("A"),c("G"),c("T","A"),c("C"))
SCGKA <- expand.grid(c("C","G"),c("C"),c("G"),c("T","G"),c("A"))
YUGUM <- expand.grid(c("C","T"),c("T"),c("G"),"T",c("T","A"))
UGYAA <- expand.grid(c("T"),c("G"),c("C","T"),"A","A")
GMUAY <- expand.grid(c("G"),c("C","A"),"T","A",c("C","T"))
CUAUN <- expand.grid(c("C"),"T","A","T",c("A","C","T","G"))
URCWC <- expand.grid(c("T"),c("G","A"),"C",c("T","A"),"C")

motif.list <- list("RRACH"=RRACH,"URUAY"=URUAY,"WEI"=WEI,"UHADG"=UHADG,"AVAYU"=AVAYU,"GGAUW"=GGAUW,"UUAKS"=UUAKS,"UUCCS"=UUCCS,
                   "AAUAA"=AAUAA,"ACUCU"=ACUCU,"UGACA"=UGACA,"AAGWC"=AAGWC,"SCGKA"=SCGKA,"YUGUM"=YUGUM,"UGYAA"=UGYAA,"GMUAY"=GMUAY,"CUAUN"=CUAUN,"URCWC"=URCWC)
motif.list <- lapply(motif.list,function(x) apply(x,1,function(y) paste0(y,collapse="")))
motif.list <- lapply(motif.list,function(x) paste0(x,collapse="|"))

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
  for(size in c(10,25,50,100,250))
  {
    #regsGR <- iclipBS[[1]][iclipBS[[1]]$type=="dominantPeak"]+size
    regsGR <- iclipBS[[1]]+size
    regsGR <- regsGR[-which(as.vector(strand(regsGR))=="*")]
    #regsGR
    #regsGR <- shift(regsGR,100)
    #seqnames(regsGR)
    oor <- which(end(ranges(regsGR)) > seqlengths(Athaliana)[as.vector(seqnames(regsGR))])
    if(length(oor)>0){
      regsGR <- regsGR[-oor]
      #end(ranges(regsGR))[oor] <- seqlengths(Athaliana)[as.vector(seqnames(regsGR[oor]))]
    }
    oor <- which(start(ranges(regsGR)) < 1)
    if(length(oor)>0)
    {
      regsGR <- regsGR[-oor]
      #start(ranges(regsGR))[oor] <- 1
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
    #sort(rrach_counts)
    
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

#devtools::install_github("omarwagih/ggseqlogo")
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
ggsave(plot=sg,filename = "PWM_fimo_all.pdf",height=12,width=12)


#################################################
###########FIMO FILES############################
#################################################

names(rrfm_freq_list) <- names(motif.list)
motif.lengths <- unlist(lapply(rrfm_freq_list,function(x) ncol(x[["10"]]$fmat)))
fmats <- lapply(rrfm_freq_list,function(x) x[["50"]]$fmat)

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
###GENOMEWIDE#################################################################
#################################################################

export PATH=/binf-isilon/alab/people/sarah/meme/bin:/binf-isilon/alab/people/sarah/meme/libexec/meme-5.1.1:$PATH

# run from /binf-isilon/alab/people/sarah/meme/results/all_from_homer 
fimo --verbosity 4 --thresh 0.05 --max-stored-scores 1000000000 ../../motif_files/motifs_all_homer.meme /binf-isilon/alab/projects/ECT2_TC/annotations/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa

################################################################
#######################OPEN AND PROCESS#########################
#################################################################

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

####genomewide##################################################
twrr <- read.table("/binf-isilon/alab/people/sarah/meme/results/all_from_homer/fimo_out/fimo.tsv",sep="\t",header=T)
cmots <- as.vector(unique(twrr$motif_id))
cmots

motGR <- list()
for(mot in cmots)
{
  twrr_4 <- twrr[(twrr$score>4) & twrr$motif_id==mot,]
  mGR <- GRanges(seqnames=Rle(twrr_4$sequence_name),strand=Rle(twrr_4$strand),IRanges(twrr_4$start+2,width=1),match_sequence=twrr_4$matched_sequence)
  motGR[[mot]] <- sort(mGR)
}

slist <- list()
for(mot in cmots)
{
  slist[[mot]] <- twrr$score[twrr$motif_id==mot]
}

motGen <- function(motG1)
{
  motG1[as.vector(strand(motG1))=="+"] <- shift(motG1[as.vector(strand(motG1))=="+"],2) + 2
  motG1[as.vector(strand(motG1))=="-"] <- shift(motG1[as.vector(strand(motG1))=="-"],-2) + 2
  return(motG1)
}
motGRR <- lapply(motGR,motGen)


#############################################
########MOTIFS GENOME WIDE###################
#############################################

for(mot in names(motGR))
{
  twGR <- motGR[[mot]]
  dr <- returnDens(twGR,gtfGR)
  dr_rand <- returnDens(randGR,gtfGR)
  pdf(paste0("normalised_enrichment_plot_",mot,".pdf"),width=6.5,height=6)
  plotNE(bs=0.02,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=T,to.add=F,col.no=1,no.samps=25,to.smooth=F,mytitle = mot)
  dev.off()
}


for(mot in names(motGR))
{
  twGR <- motGR[[mot]]
  dr <- returnDens(twGR,gtfGR)
  dr_rand <- returnDens(randGR,gtfGR)
  pdf(paste0("normalised_enrichment_plot_",mot,".pdf"),width=6.5,height=6)
  plotNE(bs=0.1,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=T,to.add=F,col.no=1,no.samps=25,to.smooth=F,mytitle = mot)
  dr <- returnDens(twGR,gtfGR[gtfGR$gene_id %in% t])
  plotNE(bs=0.1,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=T,to.add=T,col.no=1,no.samps=25,to.smooth=F,mytitle = mot)
  
  dev.off()
}


#############################################
########MATCH SETS BY GENOMIC POS############
#############################################

rmg <- union(unlist(lapply(iclipBS,function(x) x$gene)),unlist(lapply(pos.list.ECT3,function(x) x$gene)))
#nontargs <- Reduce(union,lapply(matched.set.list,function(x) x$non_targets))
nontargs <- genes.list$genes.nontargets
summary(nontargs %in% rmg)
nontargs <- nontargs[!(nontargs %in% rmg)]

nontargs <- genes.list$genes.expressed

rGR <- randGR
rGR <- rGR[which(countOverlaps(rGR,gtfGR[gtfGR$gene_id %in% nontargs])>0)]
names(rGR) <- paste0(seqnames(rGR),":",start(rGR),",",strand(rGR))


ranges(weiGRC) <- IRanges(start = start(weiGRC)-1,end=end(weiGRC))
weiGRC <- unique(weiGRC)
setGR <- list("iclip"=iclipBS[[1]],"nano"=unique(nanoGR),"mic"=unique(micGR),"wei"=unique(weiGRC),"ht"=pos.list.ECT2[["roots"]])
for(i in 1:length(setGR)){ 
  ranges(setGR[[i]]) <- IRanges(end(unique(setGR[[i]])),width=1)
  names(setGR[[i]]) <- paste0(seqnames(setGR[[i]]),":",start(setGR[[i]]),",",strand(setGR[[i]])) 
}

matchGR <- getMatchedSets(tGR=setGR,rGR = rGR,mdb = 0.1)
nset <- length(matchGR)-1

############################################
###########MAKE BED FILES FOR HOMER#########
############################################

fg <- matchGR$iclip
bg <- matchGR$random[matchGR$iclip$match]
bg <- unique(bg)

saveHomerMatch(fg+20,bg+20,my.name="iclip_20bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="iclip_10bp_homer")
saveHomerMatch(fg+5,bg+5,my.name="iclip_5bp_homer")



fg <- matchGR$wei
bg <- matchGR$random[matchGR$wei$match]
bg <- unique(bg)

saveHomerMatch(fg+20,bg+20,my.name="wei_20bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="wei_10bp_homer")
saveHomerMatch(fg+5,bg+5,my.name="wei_5bp_homer")



fg <- matchGR$nano
bg <- matchGR$random[matchGR$nano$match]
bg <- unique(bg)

saveHomerMatch(fg+5,bg+5,my.name="nano_5bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="nano_10bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="nano_20bp_homer")


fg <- matchGR$ht
bg <- matchGR$random[matchGR$ht$match]
bg <- unique(bg)

saveHomerMatch(fg+5,bg+5,my.name="ht_5bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="ht_10bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="ht_20bp_homer")

fg <- matchGR$mic
bg <- matchGR$random[matchGR$mic$match]
bg <- unique(bg)

saveHomerMatch(fg+5,bg+5,my.name="mic_5bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="mic_10bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="mic_20bp_homer")


fg <- iclipBS$E2mCh_110kDa_1.bsites.bed
bg <- iclipBS$E2mChW_110kDa_1.bsites.bed

saveHomerMatch(fg+5,bg+5,my.name="iclip_cage_110_5bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="iclip_cage_110_10bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="iclip_cage_110_20bp_homer")

fg <- iclipBS$E2mCh_55kDa_1.bsites.bed
bg <- iclipBS$E2mChW_55kDa_1.bsites.bed

saveHomerMatch(fg+5,bg+5,my.name="iclip_cage_55_5bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="iclip_cage_55_10bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="iclip_cage_55_20bp_homer")


fg <- iclipBS$E2mCh_110kDa_1.bsites.bed
bg <- iclipBS$E2mCh_55kDa_1.bsites.bed

saveHomerMatch(fg+5,bg+5,my.name="iclip_110_vs_55_5bp_homer")
saveHomerMatch(fg+10,bg+10,my.name="iclip_110_vs_55_10bp_homer")
saveHomerMatch(fg+20,bg+20,my.name="iclip_110_vs_55_20bp_homer")



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


PATH=$PATH:/binf-isilon/alab/people/sarah/homer/bin/

  
findMotifsGenome.pl iclip_5bp_homer_fg.bed tair10 resultsENWSH_iclip_5bp_homer_fg -bg iclip_5bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
    
findMotifsGenome.pl iclip_10bp_homer_fg.bed tair10 resultsENWSH_iclip_10bp_homer_fg -bg iclip_10bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
findMotifsGenome.pl iclip_20bp_homer_fg.bed tair10 resultsENWSH_iclip_20bp_homer_fg -bg iclip_20bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &

  
  findMotifsGenome.pl wei_5bp_homer_fg.bed tair10 resultsENWSH_wei_5bp_homer_fg -bg wei_5bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl wei_10bp_homer_fg.bed tair10 resultsENWSH_wei_10bp_homer_fg -bg wei_10bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl wei_20bp_homer_fg.bed tair10 resultsENWSH_wei_20bp_homer_fg -bg wei_20bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  

findMotifsGenome.pl nano_5bp_homer_fg.bed tair10 resultsENWSH_nano_5bp_homer_fg -bg nano_5bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
findMotifsGenome.pl nano_10bp_homer_fg.bed tair10 resultsENWSH_nano_10bp_homer_fg -bg nano_10bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
findMotifsGenome.pl nano_20bp_homer_fg.bed tair10 resultsENWSH_nano_20bp_homer_fg -bg nano_20bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &


  findMotifsGenome.pl ht_5bp_homer_fg.bed tair10 resultsENWSH_ht_5bp_homer_fg -bg ht_5bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl ht_10bp_homer_fg.bed tair10 resultsENWSH_ht_10bp_homer_fg -bg ht_10bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl ht_20bp_homer_fg.bed tair10 resultsENWSH_ht_20bp_homer_fg -bg ht_20bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  
  findMotifsGenome.pl mic_5bp_homer_fg.bed tair10 resultsENWSH_mic_5bp_homer_fg -bg mic_5bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl mic_10bp_homer_fg.bed tair10 resultsENWSH_mic_10bp_homer_fg -bg mic_10bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl mic_20bp_homer_fg.bed tair10 resultsENWSH_mic_20bp_homer_fg -bg mic_20bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  
  



findMotifsGenome.pl iclip_cage_55_5bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_55_5bp_homer_fg -bg iclip_cage_55_5bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl iclip_cage_55_10bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_55_10bp_homer_fg -bg iclip_cage_55_10bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl iclip_cage_55_20bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_55_20bp_homer_fg -bg iclip_cage_55_20bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  

  findMotifsGenome.pl iclip_cage_110_5bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_110_5bp_homer_fg -bg iclip_cage_110_5bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl iclip_cage_110_10bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_110_10bp_homer_fg -bg iclip_cage_110_10bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl iclip_cage_110_20bp_homer_fg.bed tair10 resultsENWSH_iclip_cage_110_20bp_homer_fg -bg iclip_cage_110_20bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  

  findMotifsGenome.pl iclip_110_vs_55_5bp_homer_fg.bed tair10 resultsENWSH_iclip_110_vs_55_5bp_homer_fg -bg iclip_110_vs_55_5bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl iclip_110_vs_55_10bp_homer_fg.bed tair10 resultsENWSH_iclip_110_vs_55_10bp_homer_fg -bg iclip_110_vs_55_10bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  findMotifsGenome.pl iclip_110_vs_55_20bp_homer_fg.bed tair10 resultsENWSH_iclip_110_vs_55_20bp_homer_fg -bg iclip_110_vs_55_20bp_homer_bg.bed -size given -rna -noweight -h -nlen 0 -len 5,6,7,8 &
  
  
#############################################
########MATCH SETS PLOT BY POSITION##########
#############################################

pdf("motif_matched_location_sets_fg_bg_all.pdf",width=60,height=12)
par(mfcol=c(nset,length(motGR)))
glmdat_list <- list()
rat_list <- list()
for(k in 1:length(motGR))
{
  getRun <- getCountDens(matchGR,motGR[[k]],tol = 10)
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


#############################################
########MATCH SETS BY GENOMIC POS############
#############################################

pvd <- list()
for(k in 1:length(motGR))
{
  pvd[[k]] <- getCountDens(matchGR,motGR[[k]])[[2]]
}
pvd <- do.call(rbind,pvd)
#pvd <- matrix(p.adjust(pvd),ncol=nset)
colnames(pvd) <- names(matchGR)[1:nset]
row.names(pvd) <- names(motGR)

library(pheatmap)
mat <- do.call(rbind,rat_list)
row.names(mat) <- names(motGR)
colnames(mat) <- names(matchGR)[1:nset]

p <- pheatmap(log2(t(mat)),display_numbers = round(t(pvd),3))

min(log2(t(mat)))
max(log2(t(mat)))

breaksList = seq(-1.5,1.5,by=0.1)
colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "PRGn")))(length(breaksList))
#colors<-colorRampPalette(rev(brewer.pal(n=7,name="PRGn")))(255)
p <- pheatmap(log2(t(mat)),col=colors,breaks=breaksList,border_color = "black",cellwidth = 15,cellheight = 15)
p

ggsave(p[[4]],file="pheatmap_motif_matched_fg_bg_enrichment_wei.pdf",width=7,height=5)


glmdat_list <- list()
for(k in 1:8)
{
  getRun <- getCountDens(matchGR,motGR[[k]],tol=25)
  glmdat_list[[k]] <- getRun[[3]]
}


length(glmdat_list)

mca <- do.call(cbind,lapply(glmdat_list,function(x) x[[4]][,1]))
head(mca)

mca <- as.data.frame(mca)
colnames(mca) <- names(motGR)
mca$resp <- as.numeric(glmdat_list[[1]][[4]][,2])-1

mnames <- expand.grid(names(motGR),names(motGR))
mnames <- mnames[!(mnames[,1]==mnames[,2]),]

pvc <- rep(0,nrow(mnames))
coev <- rep(0,nrow(mnames))
for(i in 1:nrow(mnames))
{
  mca_short <- mca[,paste(c(as.vector(mnames[i,1]),as.vector(mnames[i,2]),"resp"))]
  myres <- (glm((resp) ~ (.)^2,data=mca_short,family=binomial(link="logit")))
  pvc[i] <- summary(myres)$coefficients[4,4]
  coev[i] <- summary(myres)$coefficients[4,3]
}
pvcmat <- matrix(NA,8,8)
colnames(pvcmat) <- names(motGR)
row.names(pvcmat) <- names(motGR)
pvcmat[as.matrix(mnames)] <- pvc

coefmat <- matrix(NA,8,8)
colnames(coefmat) <- names(motGR)
row.names(coefmat) <- names(motGR)
coefmat[as.matrix(mnames)] <- coev

starmat <- pvcmat
starmat[pvcmat>=0.05] <- "ns"
starmat[pvcmat<0.05] <- c("*")
starmat[pvcmat<0.01] <- c("**")
starmat[pvcmat<0.001] <- c("***")
diag(starmat) <- ""

p <- pheatmap((coefmat),display_numbers = starmat)

p
ggsave(p[[4]],file="pheatmap_motif_pairs_glm_ht.pdf",width=4,height=3)



#############################################
########MOTIFS DE VS NON-DE genes############
#############################################
load(paste0("/binf-isilon/PBgrp/qbp693/ect2_protoplast/output_rnaseq_smartseq_analyses/pf2_de_output_smartseq.Rdata"))


myres <- (glm((resp) ~ (.)^2,data=mca,family=binomial(link="logit")))
summary(myres)



library(caret)
mdf <- data.frame(factor(mca$resp),factor(as.numeric(myres$fitted.values>0.5)),1-myres$fitted.values,myres$fitted.values)
colnames(mdf) <- c("obs","pred","0","1")
twoClassSummary(data=mdf,lev=levels(factor(mca$resp)))

library(pROC)
roc_obj <- roc(mdf$obs, mdf$`1`)
auc(roc_obj)
par(mfrow=c(1,1))
plot(roc_obj)

summary(glm((resp) ~ GGAUW * URUAY,data=mca))
table(mca$URUAY,mca$WEI)

s1 <- sum(mca$URUAY[mca$resp==1]==1 & (mca$RRACH)[mca$resp==1]==1)
s2 <- sum(mca$URUAY[mca$resp==0]==1 & (mca$RRACH)[mca$resp==0]==1)
s1/s2

g1 <- (glm((resp) ~ RRACH * URUAY,data=mca,family=binomial(link="logit")))
g2 <- (glm((resp) ~ RRACH + URUAY,data=mca,family=binomial(link="logit")))
(anova(g2,g1))
summary(g1)
summary(g2)
mytab <- table(mca$URUAY[mca$resp==1],(mca$RRACH)[mca$resp==1])
fisher.test(mytab)

###############
plot(matchGR[["random"]][matchGR[[1]]$match]$pos,matchGR[[1]]$pos)

hist(matchGR[[1]]$pos)

hist(matchGR[[2]]$pos)

hist(matchGR[["random"]]$pos)






###################


cl <- which(as.data.frame(distanceToNearest(setGR[[1]],setGR[[2]]))[,3]<50)

tGR <- setGR[[1]][cl]
motG1 <- motGR[["WEI"]]
motG2 <- motGR[["RRACH"]]

calcCo(tGR=tGR,motG1=motG1,motG2=motG2,shift.left=0,shift.right=0,window=50)


calcCo <- function(motG1,motG2,tGR,shift.left=0,shift.right=0,window=0)
{
  
  tGR[as.vector(strand(tGR))=="+"] <- shift(tGR[as.vector(strand(tGR))=="+"], + shift.right)
  tGR[as.vector(strand(tGR))=="-"] <- shift(tGR[as.vector(strand(tGR))=="-"], - shift.right)
  
  tGR[as.vector(strand(tGR))=="+"] <- shift(tGR[as.vector(strand(tGR))=="+"], - shift.left)
  tGR[as.vector(strand(tGR))=="-"] <- shift(tGR[as.vector(strand(tGR))=="-"], + shift.left)
  
  tGR <- tGR + window
  
  motG1[as.vector(strand(motG1))=="+"] <- shift(motG1[as.vector(strand(motG1))=="+"],2) + 2
  motG1[as.vector(strand(motG1))=="-"] <- shift(motG1[as.vector(strand(motG1))=="-"],-2) + 2
  
  motG2[as.vector(strand(motG2))=="+"] <- shift(motG2[as.vector(strand(motG2))=="+"],2) + 2
  motG2[as.vector(strand(motG2))=="-"] <- shift(motG2[as.vector(strand(motG2))=="-"],-2) + 2
  
  nol_1 <- countOverlaps(tGR,motG1)  
  nol_1[nol_1>0] <- 1
  
  nol_2 <- countOverlaps(tGR,motG2)
  nol_2[nol_2>0] <- 1
  
  tab <- table(nol_1 , nol_2)
  test <- fisher.test(tab)
  
  return(list("pval"=test$p.val,"odds"=test$estimate,"table"=tab))
}


tGR1 <- setGR[[2]]
tGR2 <- setGR[[1]]
motG1 <- motGR[["RRACH"]]

calcMS(motG1,tGR1,tGR2,shift.left=00,shift.right=00,window=25,tol=02)

tGR1 <- setGR[[1]]
tGR2 <- setGR[[2]]
motG1 <- motGR[["RRACH"]]

calcMS(motG1,tGR1,tGR2,shift.left=00,shift.right=00,window=25,tol=10)


calcMS <- function(motG1,tGR1,tGR2,shift.left=0,shift.right=0,window=0,tol=0)
{
  
  tGR1M <- tGR1
  tGR1M[as.vector(strand(tGR1M))=="+"] <- shift(tGR1M[as.vector(strand(tGR1M))=="+"], + shift.right)
  tGR1M[as.vector(strand(tGR1M))=="-"] <- shift(tGR1M[as.vector(strand(tGR1M))=="-"], - shift.right)
  
  tGR1M[as.vector(strand(tGR1M))=="+"] <- shift(tGR1M[as.vector(strand(tGR1M))=="+"], - shift.left)
  tGR1M[as.vector(strand(tGR1M))=="-"] <- shift(tGR1M[as.vector(strand(tGR1M))=="-"], + shift.left)
  
  tGR1M <- tGR1M + window
  
  motG1[as.vector(strand(motG1))=="+"] <- shift(motG1[as.vector(strand(motG1))=="+"],2) + 2
  motG1[as.vector(strand(motG1))=="-"] <- shift(motG1[as.vector(strand(motG1))=="-"],-2) + 2
  
  nol_1 <- countOverlaps(tGR1M,tGR2)  
  nol_1[nol_1>0] <- 1
  
  nol_2 <- countOverlaps(tGR1+tol,motG2)
  nol_2[nol_2>0] <- 1
  
  motif_ol <- nol_2
  m6A_ol <- nol_1
  
  tab <- table(m6A_ol , motif_ol)
  colnames(tab) <- c("no_mot","yes_mot")
  row.names(tab) <- c("no_m6A","yes_m6A")
  test <- fisher.test(tab)
  
  return(list("pval"=test$p.val,"odds"=test$estimate,"table"=tab,"m6A"=table(m6A_ol)))
}

motG1 <- motGR[["RRACH"]]
motG1[as.vector(strand(motG1))=="+"] <- shift(motG1[as.vector(strand(motG1))=="+"],2) + 2
motG1[as.vector(strand(motG1))=="-"] <- shift(motG1[as.vector(strand(motG1))=="-"],-2) + 2
names(setGR)

m6aRR <- (countOverlaps(setGR[[2]]+5,motG1))
m6aRR[m6aRR>0] <- 1
table(m6aRR>0)
t1 <- (countOverlaps(setGR[[1]],setGR[[2]][m6aRR==0]+10)>0)
t2 <- (countOverlaps(setGR[[1]],setGR[[2]][m6aRR==1]+10)>0)
table(t1)
table(t2)

fisher.test(table(t1,t2))

sum(t1==1)/length(t1)
sum(t2==1)/length(t2)


me <- (expr.list_ECT2$shoots.expr +expr.list_ECT2$roots.expr)/2
t.test(log2(me[setGR[[2]][m6aRR==0]$gene]+1),log2(me[setGR[[2]][m6aRR==0]$gene]+0))



#nanopore alone
dtnrr <- as.data.frame(distanceToNearest(nanoGR,motG1))[,3]
summary(dtnrr<5)

#miCLIP alone
dtnrr <- as.data.frame(distanceToNearest(micGR,motG1))[,3]
summary(dtnrr<5)

#intersection
nmGR <- nanoGR[countOverlaps(setGR[[2]],setGR[[3]])>0]
dtnrr <- as.data.frame(distanceToNearest(nmGR,motG1))[,3]
summary(dtnrr<5)


m6aRR <- (countOverlaps(nmGR+5,motG1))
m6aRR[m6aRR>0] <- 1
table(m6aRR>0)
t1 <- (countOverlaps(setGR[[1]],setGR[[2]][m6aRR==0]+10)>0)
t2 <- (countOverlaps(setGR[[1]],setGR[[2]][m6aRR==1]+10)>0)
table(t1)
table(t2)

fisher.test(table(t1,t2))


hist((dtnrr),breaks=300,col"red")
abline(v=5,lty=2)

############################################
nmGR_list <- list()
nmGR_list[[1]] <- setGR[[2]]
nmGR_list[[2]] <- setGR[[3]]
#nmGR_list[[3]] <- unique(intersect(setGR[[2]],setGR[[3]]))
nmGR_list[[3]] <- setGR$wei



par(mfrow=c(1,1))
#nmGR <- nmGR_list[[3]]

nmGR <- matchGR$nano

olcol <- do.call(cbind,lapply(motGRR,function(x) (countOverlaps(nmGR+6,x))))
olcol[olcol>0] <- 1
mpca <- princomp(olcol)
plot(mpca$loadings[,1],mpca$loadings[,2],type="n")
text(mpca$loadings[,1],mpca$loadings[,2],names(motGR))

which(rowSums(olcol)==0)
barplot(sort(colSums(olcol)),las=2,col="light blue")
mmat <- ((as.matrix(olcol[,c("URUAY","RRACH")])))
fisher.test(table(mmat[,1],(mmat[,2])))

nmGR <- nmGR_list[[3]]

##########
##########


rwins <- matchGR$nano + 50
rwW <- countOverlaps(rwins,motGR[["WEI"]])>0 
rwU <- countOverlaps(rwins,motGR[["URUAY"]])>0 

rwW[rwW>0] <- 1
rwU[rwU>0] <- 1

rwinsR <- matchGR$random[matchGR$nano$match] + 50
rwWR <- countOverlaps(rwinsR,motGR[["WEI"]])>0 
rwUR <- countOverlaps(rwinsR,motGR[["URUAY"]])>0 

rwWR[rwWR>0] <- 1
rwUR[rwUR>0] <- 1


sum(rwW*rwU)
sum(rwWR*rwUR)


fisher.test(table(rwW,rwU))

######################################################################
###############RRACH+ and RRACH- enrichments     #####################
######################################################################

nmGR <- matchGR$nano
nmGR <- matchGR$random[matchGR$nano$match]
nmGR_list <- list(matchGR$nano,matchGR$random[matchGR$nano$match])
names(nmGR_list) <- c("nanopore","matched control sites")

pdf("nanopore_RRACHp_vs_RRACHm_10.pdf",width=8,height=5)

rrn <- 16
  nmGR <- nmGR_list[[1]]
  olcol <- do.call(cbind,lapply(motGRR,function(x) (countOverlaps(nmGR+10,x))))
  olcol[olcol>0] <- 1
  olcol <- as.data.frame(olcol)
  olcol[,"none"] <- 0
  olcol[,"none"][rowSums(olcol[,c(1:ncol(olcol))[-rrn]])==0] <- 1
  colSums(olcol)
  
  sum(rowSums(olcol[,c(1:6,8)][olcol[,rrn]==0,])==0)
  
  rrNc <- colSums(olcol[olcol[,"RRACH"]==0,])[-rrn]/sum(olcol[,"RRACH"]==0)
  rrPc <- colSums(olcol[olcol[,"RRACH"]==1,])[-rrn]/sum(olcol[,"RRACH"]==1)
  
  sum(olcol[,"RRACH"]==0)
  colSums(olcol[olcol[,"RRACH"]==0,])[-rrn]
  
  sum(olcol[,"RRACH"]==1)
  colSums(olcol[olcol[,"RRACH"]==1,])[-rrn]
  
  
  nmGR <- nmGR_list[[2]]
  olcol <- do.call(cbind,lapply(motGRR,function(x) (countOverlaps(nmGR+10,x))))
  olcol[olcol>0] <- 1
  olcol <- as.data.frame(olcol)
  olcol[,"none"] <- 0
  olcol[,"none"][rowSums(olcol[,c(1:ncol(olcol))[-rrn]])==0] <- 1
  rrNcbg <- colSums(olcol[olcol[,"RRACH"]==0,])[-rrn]/sum(olcol[,"RRACH"]==0)
  rrPcbg <- colSums(olcol[olcol[,"RRACH"]==1,])[-rrn]/sum(olcol[,"RRACH"]==1)
 # rrNc <- rrNc/sum(rrNc)
 # rrPc <- rrPc/sum(rrPc)
  
  nbg <- rrNc/rrNcbg
  pbg <- rrPc/rrPcbg
  
  barplot(rbind(nbg,pbg),las=2,ylim=c(0,max(nbg)),beside = T,main=names(nmGR_list)[k],col=brewer.pal(3,"Set1")[1:2],ylab="enrichment VS background count")
  abline(h=1,lty=2)
  legend("topleft",c("RRACH-","RRACH+"),col=brewer.pal(3,"Set1")[1:2],lwd=3)
#}

dev.off()

for(i in 1:length(rrNc))
{
  tab <- round(length(matchGR$nano)*rbind(c(rrNc[i],rrNcbg[i]),c(rrPc[i],rrPcbg[i])))
  print(paste(names(rrNc)[i],fisher.test(tab)$p.value))
}


nmGR <- matchGR$nano
nmGR <- matchGR$random[matchGR$nano$match]
nmGR_list <- list(matchGR$nano,matchGR$random[matchGR$nano$match])
names(nmGR_list) <- c("nanopore","matched control sites")

aSet1 <- nmGR_list[[1]][countOverlaps(nmGR_list[[1]]+10,motGRR[["RRACH"]])>0]
aSet2 <- nmGR_list[[1]][countOverlaps(nmGR_list[[1]]+10,motGRR[["RRACH"]])==0]
aSet3 <- nmGR_list[[2]][countOverlaps(nmGR_list[[2]]+10,motGRR[["RRACH"]])>0]
aSet4 <- nmGR_list[[2]][countOverlaps(nmGR_list[[2]]+10,motGRR[["RRACH"]])==0]

set.names <- c("RRACH_+","RRACH_-","bg_+","bg_-")

colset <- brewer.pal(8,"Set1")

pdf("motif_locations_relative_m6A_sets_RRACHP_RRACHM.pdf",width=38,height=6)
mvec.ud <- c(120,10,20,15,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
par(mfrow=c(2,length(motGR)))
for(i in 1:length(motGR))
{
  makeDBplot(aSet_list=list(aSet1,aSet2,aSet3,aSet4),motGR=motGR[[i]],colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names = set.names)
}
dev.off()


############iCLIP

pdf("iclip_locations_relative_m6A_sets_RRACHP_RRACHM.pdf",width=38,height=6)
mvec.ud <- c(20,5,100,10,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
par(mfrow=c(2,length(motGR)))
for(i in 1:length(motGR))
{
  makeDBplot(aSet_list=list(aSet1,aSet2,aSet3,aSet4),motGR=iCLIP[[i]],colset=colset,my.name=names(iCLIP)[i],y.max=mvec.ud[i],my.set.names = set.names)
}
dev.off()

pdf("raw_motif_overlap_iclip_status.pdf",width=10,height=8)
par(mfrow=c(3,6))
for(i in 1:18)
{
  wol <- (unlist(lapply(iCLIP,function(x) sum(countOverlaps(x+10,motGR[[i]])>0)/length(x))))
  names(wol) <- c("55","55W","110","110W")
  barplot(wol,las=2,col=colset[1:2],main=names(motGR)[i],ylab="prop. with motif")
}
dev.off()


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


aSet1 <- nmGR_list[[1]][countOverlaps(nmGR_list[[1]]+10,motGRR[["RRACH"]])>0]
aSet2 <- nmGR_list[[1]][countOverlaps(nmGR_list[[1]]+10,motGRR[["RRACH"]])==0]
aSet3 <- nmGR_list[[2]][countOverlaps(nmGR_list[[2]]+10,motGRR[["RRACH"]])>0]
aSet4 <- nmGR_list[[2]][countOverlaps(nmGR_list[[2]]+10,motGRR[["RRACH"]])==0]
GLset <- list("R+"=aSet1,"R-"=aSet2,"bg+"=aSet3,"bg-"=aSet4)


pdf("m6A_wei_locations_relative_iCLIP.pdf",width=38,height=6)
mvec.ud <- c(100,100,20,15,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
par(mfrow=c(2,length(motGR)))
for(i in 1:length(motGR))
{
  makeDBplot(aSet_list=iCLIP,motGR=GLset[[i]],colset=colset,my.name=names(GLset)[i],y.max=mvec.ud[i],my.set.names = names(iCLIP))
}
dev.off()



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
colset <- brewer.pal(8,"Set1")

aSet <- matchGR$wei
aSet1 <- aSet
aSet2 <- matchGR$random[aSet$match]

set.names <- c("actual","background")

pdf("wei_VS_motif_matched_background.pdf",width=20,height=6)
mvec.ud <- c(120,15,65,40,200,100,350,250)
par(mfrow=c(2,8))
for(i in 1:8)
{
  makeDBplot(aSet_list=list(aSet1,aSet2),motGR=motGRR[[i]],colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names = set.names)
}
dev.off()

aSet <- matchGR$mic
aSet1 <- aSet[aSet$gene %in% genes.list$genes.union]
aSet2 <- aSet[aSet$gene %in% genes.list$genes.nontargets]
#aSet3 <- aSet[aSet$gene %in% genes.list$genes.intersect]
#aSet4 <- aSet[(aSet$gene %in% genes.list$genes.nontargets) & (aSet$gene %in% up_genes)]

aSet5 <- matchGR$random[aSet$match]

pdf("m6A_sites_mic_VS_motif_targs_matched_background.pdf",width=20,height=6)

set.names <- c("actual-targs","act-nontargs","background")
mvec.ud <- c(160,15,85,40,200,100,400,250)
par(mfrow=c(2,8))
for(i in 1:8)
{
  makeDBplot(aSet_list=list(aSet1,aSet2,aSet5),motGR=motGRR[[i]],colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names = set.names)
}
dev.off()

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
colset <- brewer.pal(8,"Set1")

#aSet <- nanoGR[countOverlaps(nanoGR+10,micGR)>0]
aSet <- nanoGR
aSet1 <- aSet[aSet$gene %in% up_genes]
aSet2 <- aSet[aSet$gene %in% nc_genes]
aSet3 <- aSet[aSet$gene %in% down_genes]
set.names <- c("up_genes","no-change","down_genes")


pdf("m6A_sites_VS_motifs_up_down.pdf",width=38,height=6)
mvec.ud <- c(120,10,20,15,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
#par(mfrow=c(2,2))
par(mfrow=c(2,length(motGRR)))
for(i in 1:length(motGRR))
{
  makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=motGR[[i]],colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names = set.names)
}
dev.off()


#aSet <- nanoGR[countOverlaps(nanoGR+25,micGR)>0]
aSet <- nanoGR
aSet1 <- aSet[aSet$gene %in% genes.list$genes.nontargets]
aSet2 <- aSet[aSet$gene %in% genes.list$genes.union[!(genes.list$genes.union %in% genes.list$genes.intersect)]]
aSet3 <- aSet[aSet$gene %in% genes.list$genes.intersect]
set.names <- c("non_targets","target(non-strict)","strict")

pdf("m6A_sites_VS_motifs_targets_non_targets.pdf",width=38,height=6)
mvec.ud <- c(120,10,20,15,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
par(mfrow=c(2,length(motGR)))
for(i in 1:length(motGR))
{
  makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=motGR[[i]],colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names=set.names)
}
dev.off()

pdf("m6A_sites_VS_motifs_targets_non_targets.pdf",width=20,height=6)
mvec.ud <- c(120,10,20,15,50,25,50,150,50,40,50,75,40,120,100,500,25,25,250,rep(300,30))
par(mfrow=c(2,8))

aSet <- nanoGR[countOverlaps(nanoGR+25,micGR)>0]
aSet1 <- aSet[aSet$gene %in% genes.list$genes.nontargets]
aSet2 <- aSet[aSet$gene %in% genes.list$genes.union[!(genes.list$genes.union %in% genes.list$genes.intersect)]]
aSet3 <- aSet[aSet$gene %in% genes.list$genes.intersect]
set.names <- c("non_targets","target(non-strict)","strict")

#sum(countOverlaps(gtf[gtf$gene_id %in%],motGR[["WEI"]]))/length(aSet1)
motGR[["WEI"]]
library(RNAeditR)
library(foreach)
library(doParallel)
#library()
motWEIGR <- addGenes(gtfGR=gtfGR,posGR = motGR[["WEI"]],ncore=30,quant = rowMeans(tpm.mat[,grep("c",colnames(tpm.mat))]),geneids = ids,assignStrand = F)

ol <- countOverlaps(motGR[["WEI"]],gtfGR[gtfGR$gene_id %in% genes.list$genes.intersect])
sum(ol)/length(ol)

gene.set <- sample(genes.list$genes.nontargets,length(genes.list$genes.intersect))
ol <- countOverlaps(motGR[["WEI"]],gtfGR[gtfGR$gene_id %in% gene.set])
sum(ol)/length(ol)


makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=iCLIP[[3]],colset=colset,my.name="iCLIP",y.max=300,levels=set.names)
makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=setGR[[4]],colset=colset,my.name="iCLIP",y.max=300,levels=set.names)


aSet <- nanoGR[countOverlaps(nanoGR+25,micGR)>0]
aSet1 <- aSet[aSet$gene %in% up_genes]
aSet2 <- aSet[aSet$gene %in% nc_genes]
aSet3 <- aSet[aSet$gene %in% down_genes]
set.names <- c("up_genes","no-change","down_genes")

makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=iCLIP[[3]],colset=colset,my.name="iCLIP",y.max=300,levels=set.names)
makeDBplot(aSet_list=list(aSet1,aSet2,aSet3),motGR=setGR[[4]],colset=colset,my.name="iCLIP",y.max=300,levels=set.names)

dev.off()


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



makeRelPlot(setGR[[1]],aSet1,to.add=F,my.col=colset[1],my.ylim=c(10,80),my.title="iCLIP")
makeRelPlot(setGR[[1]],aSet2,to.add=T,my.col=colset[2],my.ylim=c(8,25))
makeRelPlot(setGR[[3]],aSet3,to.add=T,my.col=colset[3],my.ylim=c(8,25))



makeRelPlot(setGR[[4]],aSet1,to.add=F,my.col=colset[2],my.ylim=c(0,40),my.title="hyperTRIBE")
makeRelPlot(setGR[[4]],aSet2,to.add=T,my.col=colset[1],my.ylim=c(8,25))


npg <- tapply(nanoGR$gene,nanoGR$gene,length)
dg <- npg[names(npg) %in% down_genes]
ug <- npg[names(npg) %in% up_genes]
dens <- density(dg)
plot(dens)
dens <- density(ug)



down_genes <- resLFC_pf2_df_smartseq[resLFC_pf2_df_smartseq$log2FoldChange<0,]$gene_id
up_genes <- resLFC_pf2_df_smartseq[resLFC_pf2_df_smartseq$log2FoldChange>0,]$gene_id
resLFC_pf2_df_smartseq$target <- "non_target"
resLFC_pf2_df_smartseq$target[resLFC_pf2_df_smartseq$gene_id %in% genes.list$genes.union] <- c("target")
#resLFC_pf2_df_smartseq$target[resLFC_pf2_df_smartseq$gene_id %in% genes.list$genes.intersect] <- c("strict")


groups <- paste(resLFC_pf2_df_smartseq$log2FoldChange>0.25,resLFC_pf2_df_smartseq$target,sep="-")
names(groups) <- resLFC_pf2_df_smartseq$gene_id
groups <- factor(groups)
levels(groups) <- c("down-non_target","down-target","up-non_target","up-target")

colset <- brewer.pal(length(unique(groups)),"Set1")

par(mfrow=c(1,1))
ug <- tapply(nanoGR$gene,nanoGR$gene,length)
pot <- cumsum(sort(ug))/length(ug)
pot <- quantile(pot,seq(0,1,by=0.001))
plot(pot,col=colset[c(2)],ylim=c(0,5),xlab="Quantile",ylab="No. m6A sites per gene",type="n")

for(i in 1:length(levels(groups)))
{
  dg <- npg[names(npg) %in% names(groups[groups==levels(groups)[i]])]
  pot <- cumsum(sort(dg))/length(dg)
  pot <- quantile(pot,seq(0,1,by=0.001))
  points(pot,type="l",col=colset[c(i)],lwd=2)
}
legend("topleft",levels(groups),col=colset,lwd=2)














aSet <- matchGR$mic
aSet1 <- aSet[aSet$gene %in% genes.list$genes.union]
aSet2 <- aSet[aSet$gene %in% genes.list$genes.nontargets]
#aSet3 <- aSet[aSet$gene %in% genes.list$genes.intersect]
#aSet4 <- aSet[(aSet$gene %in% genes.list$genes.nontargets) & (aSet$gene %in% up_genes)]

aSet5 <- matchGR$random[aSet$match]

weimotGR <- motGRR[["WEI"]]

rmatch <- matchGR$random[matchGR$iclip$match]
rseqs <- randGR[sample(1:length(randGR),length(iCLIP[[3]]))]
aSets <- c(iCLIP,rmatch,rseqs)

set.names <- c(names(iCLIP),"matched_random","any_random")
mvec.ud <- c(180,15,85,40,200,100,400,250)
par(mfrow=c(2,8))
for(i in 1)
{
  par(mfrow=c(1,2))
  makeDBplot(aSet_list=aSets,motGR=weimotGR,colset=colset,my.name=names(motGRR)[i],y.max=mvec.ud[i],my.set.names = set.names)
}


nr <- countOverlaps(matchGR$nano+10,motGR$RRACH)
nr[nr>0] <- 1
nu <- countOverlaps(matchGR$nano+10,motGR$GGAUW)
nu[nu>0] <- 1

barplot(as.vector(ftable(evid,nu,nr)))
rat <- (as.vector(table(nu[evid==1],nr[evid==1]))/as.vector(table(nu[evid==0],nr[evid==0])))/(sum(evid==1)/sum(evid==0))
rat
barplot(rat,ylab="Enrichment in 'targets'")

abline(h=1,lty=2)
RRACH==0 and GGAUW==1 = more likely to be a target than non-target

olmat <- do.call(cbind, lapply(motGR,function(x) countOverlaps(matchGR$nano+10,x)))
olmat[olmat>0] <- 1

evid <- as.numeric(matchGR$nano$gene %in% genes.list$genes.union)
table(evid)
evid <- as.numeric(matchGR$nano$gene %in% genes.list$genes.intersect)

ga <- which( (matchGR$nano$gene %in% genes.list$genes.union) & !(matchGR$nano$gene %in% genes.list$genes.intersect))
#evid <- evid[-ga]

summary(evid)


tab0 <- table(nu[evid==0],nr[evid==0])
fisher.test(tab0)
tab1 <- table(nu[evid==1],nr[evid==1])
fisher.test(tab1)

round(tab0/sum(tab0),3)
round(tab1/sum(tab1),3)

#proportions of RRACH and GGAUW


data <- as.data.frame(olmat)
data <- as.data.frame(apply(olmat,2,as.factor))
resp <- evid==1

barplot(sort(tapply(evid,paste(data$RRACH,data$GGAUW,data$URUAY),mean)))
myglm <- glm(evid~(data$RRACH+data$GGAUW+data$URUAY)^2,family="binomial")
summary(myglm)

tapply(rowSums(olmat),evid,mean)

tapply(evid,paste(data$RRACH,data$URUAY),mean)

#data <- data[-ga,]

#msamp <- sample(which(resp==0),length(resp[resp==1]),replace=T)
#resp2 <- c(resp[resp==1],resp[msamp])
#data2 <- data[c(which(resp==1),msamp),]

summary(resp)
mygm <- (glm(resp[data$RRACH==1] ~ (.),data=data[data$RRACH==1,-which(colnames(data)=="RRACH")],family="binomial"))
summary(mygm)
anova(mygm, test="Chisq")
roc_obj <- roc(resp[data$RRACH==1], mygm$fitted.values)
auc(roc_obj)
par(mfrow=c(1,1))
plot(roc_obj)

#RRACH <- data$RRACH
mygm <- (glm(resp==TRUE ~ RRACH,data=data,family="binomial"))
summary(mygm)



mygm <- (glm(resp[data$RRACH==0] ~ (.),data=data[data$RRACH==0,-which(colnames(data)=="RRACH")],family="binomial"))
anova(mygm, test="Chisq")

roc_obj <- roc(resp, mygm$fitted.values)
auc(roc_obj)
par(mfrow=c(1,1))
plot(roc_obj)




library(dominanceanalysis)
dom <- dominanceAnalysis(mygm)
dominanceMatrix(dom, type="complete",fit.functions = "r2.m", ordered=TRUE)


summary(mygm)
mygm$coefficients
table(mygm$fitted)

library(randomForest)
myrf <- randomForest(y=as.factor(as.numeric(resp2)),x=data2,trees=500,importance=TRUE,norm.votes=T)


library(pROC)
roc_obj <- roc(factor(evid), mygm$fitted.values)
#roc_obj <- roc(factor(evid), myrf$votes[,1])
auc(roc_obj)
par(mfrow=c(1,1))
plot(roc_obj)

#myrf$importance
