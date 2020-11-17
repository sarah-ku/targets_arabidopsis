setwd("/binf-isilon/alab/people/sarah/meme/motif_plots")
load("/home/sarah/github/targets_arabidopsis/targets_paper_open.R")

library(rtracklayer)
library(GenomicRanges)
library(RColorBrewer)
library("BSgenome.Athaliana.TAIR.TAIR9")
seqnames(Athaliana) <- c("1","2","3","4","5","Mt","Pt")

#################################################
###########DEFINE MOTIFS#########################
#################################################

#DRACH (D=G/A/U, R=G/A, H=A/U/C) 
# R = G/A, H = A/C/U, D=A/G/U, V=C/A/G, S=C/G, Y=C/T, W=T/A, K=T/G

RRACH <- expand.grid(c("G","A"),c("G","A"),c("A"),c("C"),c("A","T","C"))
UHADG <- expand.grid(c("T"),c("A","C","T"),"A",c("A","G","T"),"G")
AVAYU <- expand.grid(c("A"),c("C","A","G"),c("A"),c("C","T"),c("T"))
GGAUW <- expand.grid("G","G","A","T",c("T","A"))
UUAKC <- expand.grid("T","T","A",c("T","G"),"C")
UUCCS <- expand.grid("T","T","C","C",c("C","G")) 
URUAY <- expand.grid(c("T"),c("G","A"),c("T"),c("A"),c("T","C"))
WEI <- expand.grid(c("T"),c("G","T","A","C"),c("T"),c("G","A"),c("T"),c("A","G"),c("T","A"))

motif.list <- list("RRACH"=RRACH,"URUAY"=URUAY,"WEI"=WEI,"UHADG"=UHADG,"AVAYU"=AVAYU,"GGAUW"=GGAUW,"UUAKC"=UUAKC,"UUCCS"=UUCCS)
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
    regsGR <- iclipBS[[1]][iclipBS[[1]]$type=="dominantPeak"]+size
    #regsGR <- iclipBS[[1]]+size
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

almat <- lapply(rrfm_freq_list,function(x) x[["25"]]$fmat)
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
fmats <- lapply(rrfm_freq_list,function(x) x[["25"]]$fmat)

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
gtf_expr <- unique(gtfGR[gtfGR$gene_id %in% genes.list$genes.intersect])
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
  mGR <- GRanges(seqnames=Rle(twrr_4$sequence_name),strand=Rle(twrr_4$strand),IRanges(twrr_4$start,width=1),match_sequence=twrr_4$matched_sequence)
  motGR[[mot]] <- sort(mGR)
}

slist <- list()
for(mot in cmots)
{
  slist[[mot]] <- twrr$score[twrr$motif_id==mot]
}


names(motGR)
twGR <- motGR[["AVAYU"]]
dr <- returnDens(twGR,gtfGR)
dr_rand <- returnDens(randGR,gtfGR)

utab <- rev(sort(table(as.vector(twGR$match_sequence))))

dr_ind <- list()
for(i in 1:6)
{
  dr_ind[[i]] <- returnDens(twGR[twGR$match_sequence==names(utab)[i]],gtfGR)
}
names(dr_ind) <- names(utab)[1:6]
dr <- returnDens(twGR,gtfGR)
dr_rand <- returnDens(randGR,gtfGR)

plotDens(bs=0.01,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=T,to.add=F,col.no=1,no.samps=10)



for(mot in names(motGR))
{
  twGR <- motGR[[mot]]
  dr <- returnDens(twGR,gtfGR)
  dr_rand <- returnDens(randGR,gtfGR)
  pdf(paste0("normalised_enrichment_plot_",mot,".pdf"),width=6.5,height=6)
  plotNE(bs=0.02,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=T,to.add=F,col.no=1,no.samps=25,to.smooth=F,mytitle = mot)
  dev.off()
}


