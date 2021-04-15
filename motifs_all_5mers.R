N <- c("C","T","A","G")
M <- c("C","A")

motif.list <- list()

all5 <- expand.grid(N,N,N,N,N)
all5 <- apply(all5,1,function(x) paste(x,collapse=""))
length(all5)

bs <- iCLIPBS[[3]]
bs <- bs[as.vector(strand(bs)) %in% c("+","-")]
names(bs) <- paste(as.vector(seqnames(bs)),":",start(bs)-1,",",as.vector(strand(bs)),sep="")
cnames <- intersect(names(matchGR$iclip),names(bs))
bs <- bs[cnames]
bgp <- matchGR$iclip[cnames]$match
bgp <- bgp[which(bgp %in% names(matchGR$random))]
inames <- names(matchGR$iclip[cnames])[which(bgp %in% names(matchGR$random))]

bs_bg <- matchGR$random[bgp]
bs <- (bs[inames])

length(bs)
length(bs_bg)

fs <- as.vector(getSeq(Athaliana,bs+20))
fs_bg <- as.vector(getSeq(Athaliana,bs_bg+20))

library(stringr)
glm.res <- matrix(NA,ncol=5,nrow=length(all5))
for(i in 1:length(all5))
{
  resp <- factor(c(rep(1,length(fs)),rep(0,length(fs_bg))))
  mot <- c(str_count(fs, all5[i]),str_count(fs_bg, all5[i]))
  if(sum(mot)>0)
  {
    myglm <- (glm(resp~mot,family="binomial"))
    counts <- c(sum(mot[resp==0]),sum(mot[resp==1]))
    glm.res[i,] <- c(summary(myglm)$coefficients[2,3:4],myglm$coefficients[2],counts)
  }
}

#glm.res[,3] <- exp(glm.res[,3])
res.names <- all5
res.names <- gsub("T","U",res.names)
row.names(glm.res) <- res.names

glm.res[order(glm.res[,1]),][1:20,]
glm.res[rev(order(glm.res[,1],na.last = F)),][1:20,]

row.names(glm.res[rev(order(glm.res[,1]))[1:20],])
row.names(glm.res[(order(glm.res[,1]))[1:20],])


pdf("/binf-isilon/alab/projects/ECT2_TC/motif_plots/most_frequent_least_frequent_bg_fg_iCLIP_20.pdf",width=14,height=7)
par(mfrow=c(2,1))
barplot(glm.res[rev(order(glm.res[,4])),][1:100,4],las=2,cex.names =.75,col=rgb(1,0,0,alpha=0.5),ylim=c(0,max(glm.res[,4:5])),main="Most common (ordered by sites at matched BG (red))")
barplot(glm.res[rev(order(glm.res[,4])),][1:100,5],las=2,cex.names =.75,add=T,col=rgb(0,0,1,alpha=0.5))

barplot(glm.res[(order(glm.res[,4])),][1:100,4],las=2,cex.names=.75,col=rgb(1,0,0,alpha=0.5),main="Least common (ordered by sites at matched BG (red))")
barplot(glm.res[(order(glm.res[,4])),][1:100,5],las=2,cex.names =.75,add=T,col=rgb(0,0,1,alpha=0.5))
dev.off()


