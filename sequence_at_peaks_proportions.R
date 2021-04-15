setGR <- list("55kDa"=iCLIP[[1]],"55kWDa"=iCLIP[[2]],"110kDa"=iCLIP[[3]],"110kWDa"=iCLIP[[4]],"55kDa_BS"=iCLIPBS[[1]],"55kWDa_BS"=iCLIPBS[[2]],"110kDa_BS"=iCLIPBS[[3]],"110kWDa_BS"=iCLIPBS[[4]],"miCLIP"=micGR,"nanopore"=nanoGR)

stab <- matrix(NA,nrow=length(setGR),ncol=4)
row.names(stab) <- names(setGR)
colnames(stab) <- c("A","T","C","G")


for(i in 1:length(setGR))
{
  mtab <- table(as.vector(getSeq(Athaliana,setGR[[i]][as.vector(strand(setGR[[i]])) %in% c("-","+")])))
  stab[i,] <- mtab[c("A","T","C","G")]
}

#write.table(stab,file="/binf-isilon/alab/people/sarah/iCLIP_nt_counts.txt",quote=FALSE,sep="\t")

proptab <- apply(stab,1,function(x) x/sum(x))
barplot(proptab,beside = T,las=2)
colnames(proptab) <- names(setGR)

mydf <- melt(proptab)
colnames(mydf) <- c("base","exper","prop")
myp <- ggplot(mydf,aes(x=base,y=prop,fill=exper)) + geom_bar(stat = "identity",col="black") + facet_grid(~exper) + theme_bw() + ylab("Proportion of peaks at reference nucleotide")
myp
ggsave(myp,file="/binf-isilon/alab/projects/ECT2_TC/iCLIP/proportion_peaks_bases.pdf",width=8,height=4)