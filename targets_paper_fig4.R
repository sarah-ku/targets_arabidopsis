source("/home/sarah/github/targets_arabidopsis/targets_paper_open.R")
source("/home/sarah/github/targets_arabidopsis/motifs_functions.R")
setwd("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/")
# pl <- pos.list.ECT2[["roots"]]
# expressed_genes_roots <- genes.list$genes.roots.expressed
# 
# mycuts <- cut2(log2(c(expr.list_ECT2$roots.expr,expr.list_ECT2$shoots.expr)+1),g=9,onlycuts = T)
# 
# g10 <- as.numeric(cut2(log(expr.list_ECT2$roots.expr+1),cuts = mycuts))
# names(g10) <- names(expr.list_ECT2$roots.expr)
# table(g10)

###################################################################
##############MAKE DATA.FRAME######################################
###################################################################

boxdf <- data.frame("gene"=genes.list$genes.expressed)
row.names(boxdf) <- boxdf$gene

#append data from other sources
for(i in 1:length(genes.list))
{
  boxdf[,names(genes.list)[i]] <- 0
  boxdf[intersect(genes.list[[i]],row.names(boxdf)),names(genes.list)[i]] <- 1
}
boxdf[,"nanopore"] <- 0
boxdf[,"miclip"] <- 0
boxdf[,"faclip"] <- 0
boxdf[,"shen"] <- 0
boxdf[intersect(row.names(boxdf),unique(as.vector(nanoGR$gene))),"nanopore"] <- 1
boxdf[intersect(row.names(boxdf),unique(as.vector(micGR$gene))),"miclip"] <- 1
boxdf[intersect(row.names(boxdf),unique(as.vector(weiGRC$gene))),"faclip"] <- 1
boxdf[intersect(row.names(boxdf),unique(as.vector(m6AGR$genes[,"genes"]))),"shen"] <- 1


boxdf[,"ht_supported"] <- 0
boxdf[,"ht_unsupported"] <- 0
boxdf[,"faclip_supported"] <- 0
boxdf[,"faclip_unsupported"] <- 0

summary(unique(unique(m6AGR$genes[,"genes"])) %in% c(nanoGR$gene,micGR$gene))

m6A_genes <- unique(c(nanoGR$gene,micGR$gene,unique(m6AGR$genes[,"genes"])))
ht_supported <- genes.list$genes.union[genes.list$genes.union %in% m6A_genes]
ht_unsupported <- genes.list$genes.union[!(genes.list$genes.union %in% m6A_genes)]

fa_supported <- unique(weiGRC$gene[weiGRC$gene %in% m6A_genes])
fa_unsupported <- unique(weiGRC$gene[!(weiGRC$gene %in% m6A_genes)])
length(fa_unsupported)

boxdf[intersect(row.names(boxdf),ht_supported),"ht_supported"] <- 1
boxdf[intersect(row.names(boxdf),ht_unsupported),"ht_unsupported"] <- 1
boxdf[intersect(row.names(boxdf),fa_supported),"faclip_supported"] <- 1
boxdf[intersect(row.names(boxdf),fa_unsupported),"faclip_unsupported"] <- 1

for(i in 1:length(iCLIP))
{
  boxdf[,names(iCLIP)[i]] <- 0
  boxdf[intersect(row.names(boxdf),unique(as.vector(iCLIP[[i]]$gene))),names(iCLIP)[i]] <- 1
}

###################################################################
##############APPEND EXPRESSION####################################
###################################################################

boxdf$ebin <- NA
boxdf$ebin.roots <- NA
boxdf$ebin.shoots <- NA

roots_expr <- log(log(expr.list_ECT2$roots.expr+1)+.1)
hist(roots_expr)
roots_expr <- roots_expr[genes.list$genes.roots.expressed]
my_cut <- cut(roots_expr,9)
#table(my_cut)
my_cut <- as.numeric(my_cut)
names(my_cut) <- genes.list$genes.roots.expressed
cgenes <- intersect(row.names(boxdf),genes.list$genes.roots.expressed)
boxdf[cgenes,"ebin.roots"] <- my_cut[cgenes]

shoots_expr <- log(log(expr.list_ECT2$shoots.expr+1)+.1)
shoots_expr <- shoots_expr[genes.list$genes.shoots.expressed]
my_cut <- cut(shoots_expr, 9)
my_cut <- as.numeric(my_cut)
names(my_cut) <- genes.list$genes.shoots.expressed
cgenes <- intersect(row.names(boxdf),genes.list$genes.shoots.expressed)
boxdf[cgenes,"ebin.shoots"] <- my_cut[cgenes]

expr <- log(log(c(expr.list_ECT2$roots.expr+expr.list_ECT2$shoots.expr)/2+1)+.1)
expr <- expr[genes.list$genes.expressed]
my_cut <- cut(expr, 9)
levels(my_cut)
my_cut <- as.numeric(my_cut)
names(my_cut) <- genes.list$genes.expressed
cgenes <- intersect(row.names(boxdf),genes.list$genes.expressed)
boxdf[cgenes,"ebin"] <- my_cut[cgenes]

table(boxdf$ebin)
table(boxdf$ebin.roots)
table(boxdf$ebin.shoots)

boxdf$class_dummy <- 1

###################################################################
##############MAKE GGPLOTS#########################################
###################################################################

#df.roots <- boxdf[,c("genes.roots","genes.roots.intersect","genes.nontargets.roots")]
df.roots <- boxdf[,c("class_dummy","genes.roots.intersect","genes.nontargets.roots")]
df.sum.roots <- apply(df.roots,2,function(x) tapply(x,boxdf$ebin.roots,function(y) sum(y)/sum(x)))
df.sum.roots <- apply(df.roots,2,function(x) tapply(x,boxdf$ebin.roots,function(y) sum(y)))
df.sum.roots <- melt(df.sum.roots)
colnames(df.sum.roots) <- c("expr","name","prop")
levels(df.sum.roots$name)[1] <- "class_dummy.roots"


df.shoots <- boxdf[,c("class_dummy","genes.shoots")]
df.sum.shoots <- apply(df.shoots,2,function(x) tapply(x,boxdf$ebin.shoots,function(y) sum(y)/sum(x)))
df.sum.shoots <- apply(df.shoots,2,function(x) tapply(x,boxdf$ebin.shoots,function(y) sum(y)))
df.sum.shoots <- melt(df.sum.shoots)
colnames(df.sum.shoots) <- c("expr","name","prop")
levels(df.sum.shoots$name)[1] <- "class_dummy.shoots"

##########################################
colnames(boxdf)
boxdf$faclip_union <- boxdf$faclip*boxdf$genes.union
boxdf$union_faclip <- boxdf$faclip*boxdf$genes.union
boxdf$nano_mic_shen <- boxdf$nanopore+boxdf$miclip+boxdf$shen 
boxdf$nano_mic_shen[boxdf$nano_mic_shen>0] <- 1
boxdf$parker <- boxdf$nanopore+boxdf$miclip
boxdf$parker[boxdf$parker>0] <- 1


df <- boxdf[,c("class_dummy","parker","miclip","nanopore","genes.union","faclip","nano_mic_shen","shen","E2_110kGR","E2_110kWGR")]
#df <- boxdf[,c("class_dummy","ht_supported","ht_unsupported","faclip_supported","genes.union","faclip","union_faclip","faclip_union","E2_110kGR","E2_110kWGR")]
#df <- boxdf[,c("class_dummy","genes.union","genes.nontargets","nanopore","miclip","faclip","E2_55kGR","E2_55kWGR","E2_110kGR","E2_110kWGR")]
#df <- boxdf[,c("genes.single","genes.union","genes.nontargets","nanopore","miclip","faclip","E2_55kGR","E2_55kWGR","E2_110kGR","E2_110kWGR")]
#df.sum <- apply(df,2,function(x) tapply(x,boxdf$ebin,function(y) sum(y)/sum(x)))
#df.sum <- apply(df,2,function(x) tapply(x,boxdf$ebin,function(y) sum(y)))
df.sum <- apply(df,2,function(x) tapply(x,boxdf$ebin,function(y) sum(y)/length(y)))

# for(i in 2:9)
# {
#   df.sum[i,"faclip_union"] <- sum(boxdf$faclip_union[boxdf$ebin==i])/sum(boxdf$faclip[boxdf$ebin==i])
# }
# 
# for(i in 2:9)
# {
#   df.sum[i,"union_faclip"] <- sum(boxdf$faclip_union[boxdf$ebin==i])/sum(boxdf$genes.union[boxdf$ebin==i])
# }
# 


df.sum <- melt(df.sum)
colnames(df.sum) <- c("expr","name","prop")
levels(df.sum$name)[1] <- "class_dummy.all"



df.sum <- rbind(df.sum.shoots,df.sum.roots,df.sum)

cols <- rev(heat.colors(10))
names(cols) <- 1:10
df.sum$col <- cols[df.sum$expr]

p <- ggplot(data = df.sum,aes(x=factor(expr),y=prop,fill=factor(expr))) + geom_bar(stat="identity",colour="black") + facet_grid(rows=vars(df.sum$name),scales = "free_y")
p <- p + theme_classic() + scale_fill_manual(values=cols) + ylab("Raw proportion of genes supported by method") + xlab("Expression bin")
p


ggsave(p,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/proportion_of_class_by_expression_NANOPORE_miCLIP_SHEN.pdf",width=3.75,height=16)


###################################################################################
###################################################################################

##############################################
#########ECT2 roots###########################
##############################################
roots_expr_HT <- (log2(expr.list_ECT2$roots.expr+1))
roots_expr_HT <- roots_expr_HT[genes.list$genes.roots]
hist(roots_expr_HT)
my_cut <- cut(roots_expr_HT,9)
table(my_cut)
my_cut <- as.numeric(my_cut)
names(my_cut) <- genes.list$genes.roots
cgenes <- intersect(row.names(boxdf),genes.list$genes.roots)
boxdf[cgenes,"ebin.roots.HT"] <- my_cut[cgenes]
pl <- pos.list.ECT2[["roots"]]

q <- createBinByHits(pl,boxdf)
q
ggsave(q,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/expression_bins_by_hits_pval_prop_ECT2_roots.pdf",width=7,height=8)

##############################################
#########ECT2 roots###########################
##############################################
shoots_expr_HT <- (log2(expr.list_ECT2$shoots.expr+1))
shoots_expr_HT <- shoots_expr_HT[genes.list$genes.shoots]
hist(shoots_expr_HT)
my_cut <- cut(shoots_expr_HT,9)
table(my_cut)
my_cut <- as.numeric(my_cut)
names(my_cut) <- genes.list$genes.shoots
cgenes <- intersect(row.names(boxdf),genes.list$genes.shoots)
boxdf[cgenes,"ebin.shoots.HT"] <- my_cut[cgenes]
pl <- pos.list.ECT2[["shoots"]]

q <- createBinByHits(pl,boxdf)
ggsave(q,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/expression_bins_by_hits_pval_prop_ECT2_shoots.pdf",width=7,height=8)

##############################################
#########ECT3 roots###########################
##############################################
roots_expr_HT <- (log2(expr.list_ECT3$roots.expr+1))
roots_expr_HT <- roots_expr_HT[genes.list$genes.roots]
hist(roots_expr_HT)
my_cut <- cut(roots_expr_HT,9)
table(my_cut)
my_cut <- as.numeric(my_cut)
names(my_cut) <- genes.list$genes.roots
cgenes <- intersect(row.names(boxdf),genes.list$genes.roots)
boxdf[cgenes,"ebin.roots.HT"] <- my_cut[cgenes]
pl <- pos.list.ECT3[["roots"]]

q <- createBinByHits(pl,boxdf)
ggsave(q,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/expression_bins_by_hits_pval_prop_ECT3_roots.pdf",width=7,height=8)


##############################################
#########ECT3 shoots###########################
##############################################
shoots_expr_HT <- (log2(expr.list_ECT3$shoots.expr+1))
shoots_expr_HT <- shoots_expr_HT[genes.list$genes.shoots]
hist(shoots_expr_HT)
my_cut <- cut(shoots_expr_HT,9)
table(my_cut)
my_cut <- as.numeric(my_cut)
names(my_cut) <- genes.list$genes.shoots
cgenes <- intersect(row.names(boxdf),genes.list$genes.shoots)
boxdf[cgenes,"ebin.shoots.HT"] <- my_cut[cgenes]
pl <- pos.list.ECT3[["shoots"]]

q <- createBinByHits(pl,boxdf)
ggsave(q,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/expression_bins_by_hits_pval_prop_ECT3_shoots.pdf",width=7,height=8)




createBinByHits <- function(pl,boxdf)
{
  min.p <- tapply(pl$padj,pl$gene,function(x) -log10(min(x,na.rm=T)))
  mean.p <- tapply(pl$padj,pl$gene,function(x) -log10(mean(x,na.rm=T)))
  n.hits <- tapply(pl$prop,pl$gene,function(x) length(x))
  max.prop <- tapply(pl$prop,pl$gene,function(x) max(x))
  mean.prop <- tapply(pl$prop,pl$gene,function(x) mean(x))
  
  boxdf$n.hits <- NA
  boxdf[names(n.hits),"n.hits"] <- n.hits
  
  boxdf$mean.prop <- NA
  boxdf[names(mean.prop),"mean.prop"] <- mean.prop
  
  boxdf$max.prop <- NA
  boxdf[names(max.prop),"max.prop"] <- max.prop
  
  boxdf$min.p <- NA
  boxdf[names(min.p),"min.p"] <- min.p
  
  boxdf$mean.p <- NA
  boxdf[names(mean.p),"mean.p"] <- mean.p
  
  
  df <- boxdf[!is.na(boxdf$ebin.roots.HT),c("ebin.roots.HT","n.hits","max.prop","mean.prop","min.p","mean.p")]
  df.p <- melt(df[,-1])
  df.p$expr <- rep(df$ebin.roots.HT,5)
  
  #df.p <- df.p[-which(df.p$variable=="n.hits" & df.p$value<30),]
  p1 <- ggplot(data = df.p[df.p$variable=="n.hits",],aes(x=factor(expr),y=value,fill=factor(expr))) + geom_boxplot(outlier.shape=NA) + ylim(0,30)
  p1 <- p1 + theme_classic() + scale_fill_manual(values=cols) + ylab("number of hits in gene") + xlab("Expression bin")
  
  p2 <- ggplot(data = df.p[df.p$variable=="max.prop",],aes(x=factor(expr),y=value,fill=factor(expr))) + geom_boxplot(outlier.shape=NA) + ylim(0,0.6)
  p2 <- p2 + theme_classic() + scale_fill_manual(values=cols) + ylab("editing proportion (max in gene)") + xlab("Expression bin")
  
  p5 <- ggplot(data = df.p[df.p$variable=="mean.prop",],aes(x=factor(expr),y=value,fill=factor(expr))) + geom_boxplot(outlier.shape=NA) + ylim(0,0.6)
  p5 <- p5 + theme_classic() + scale_fill_manual(values=cols) + ylab("editing proportion (mean in gene)") + xlab("Expression bin")
  
  p3 <- ggplot(data = df.p[df.p$variable=="min.p",],aes(x=factor(expr),y=value,fill=factor(expr))) + geom_boxplot(outlier.shape=NA) + ylim(0,175)
  p3 <- p3 + theme_classic() + scale_fill_manual(values=cols) + ylab("-log10(min adj p-value)") + xlab("Expression bin")
  
  p4 <- ggplot(data = df.p[df.p$variable=="mean.p",],aes(x=factor(expr),y=value,fill=factor(expr))) + geom_boxplot(outlier.shape=NA) + ylim(0,10)
  p4 <- p4 + theme_classic() + scale_fill_manual(values=cols) + ylab("-log10(mean adj p-value)") + xlab("Expression bin") 
  
  q <- grid.arrange(p1,p2,p5,p3,p4)
  return(q)
}



####################################################################################
######HYPERTRIBE VS ICLIP PEAKS#####################################################
####################################################################################

mycols <- brewer.pal(7,"Set1")
htl <- unique(unlist(GRangesList(pos.list.ECT2[1:2])))

dtn <- (mcols(distanceToNearest(htl,iCLIP[[3]]))$distance)
dtn <- dtn[dtn<750]
which.max(table(dtn))
md <- (density(dtn))
which.max(md$y)

mycols <- brewer.pal(7,"Set1")
pdf("enrichment_of_ht_around_iclip.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
pSets <- iCLIP
makeRelPlotSets(htl,pSets,myn=1000,my.col=mycols,my.title=paste0("hyperTRIBE vs iCLIP"),use.para=TRUE)
legend("topright",c(names(iCLIP)),col=mycols[1:4],lty=1,lwd=1)
abline(v=1001-80,lty=2)
dev.off()


#percentage of iCLIP peaks covered by hyperTRIBE within distance X
iclip_cov_110 <- rep(NA,1001)
iclip_cov_110W <- rep(NA,1001)
iclip_cov_55 <- rep(NA,1001)
iclip_cov_55W <- rep(NA,1001)
for(i in 0:1000)
{
  iclip_cov_110[i+1] <- sum(countOverlaps(iCLIP[[3]]+i,htl)>0)/length(iCLIP[[3]])
  iclip_cov_110W[i+1] <- sum(countOverlaps(iCLIP[[4]]+i,htl)>0)/length(iCLIP[[4]])
  iclip_cov_55[i+1] <- sum(countOverlaps(iCLIP[[1]]+i,htl)>0)/length(iCLIP[[1]])
  iclip_cov_55W[i+1] <- sum(countOverlaps(iCLIP[[2]]+i,htl)>0)/length(iCLIP[[2]])
}

#getwd()
pdf("hyperTRIBE_cumulative_near_iCLIP_peaks_proportions.pdf",width=5,height=5.5)
plot(iclip_cov_110,type="l",col=mycols[1],lwd=2,ylab="Proportion of iCLIP peaks supported by hyperTRIBE",xlab="Distance of hyperTRIBE peaks from iCLIP peaks (nt)")
points(iclip_cov_110W,type="l",lty=2,col=mycols[1],lwd=2)
points(iclip_cov_55,type="l",col=mycols[2],lwd=2)
points(iclip_cov_55W,type="l",lty=2,col=mycols[2],lwd=2)
abline(v=200,lty=3)
abline(h=iclip_cov_110[201],lty=3)
text(5,iclip_cov_110[201],round(iclip_cov_110[201],3))
abline(h=iclip_cov_110W[201],lty=3)
text(5,iclip_cov_110W[201],round(iclip_cov_110W[201],3))
abline(h=iclip_cov_55[201],lty=3)
text(5,iclip_cov_55[201],round(iclip_cov_55[201],3))
abline(h=iclip_cov_55W[201],lty=3)
text(5,iclip_cov_55W[201],round(iclip_cov_55W[201],3))
legend("right",names(iCLIP),col=c(mycols[c(2,2,1,1)]),lty=c(1,2,1,2),lwd=2)
dev.off()

####################################################################################
######peaks sets around other peak sets#############################################
####################################################################################

andGRC <- andGR
start(andGRC) <- round(start(andGR) + 0.5*(end(andGR)-start(andGR)))
end(andGRC) <- round(start(andGR) + 0.5*(end(andGR)-start(andGR)))

pdf("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/enrichment_of_sets_around_m6A.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(iCLIP[[3]],nanoGR,myn=250,to.add=F,my.col=mycols[1],my.ylim=c(0,0.013),my.title="Nanopore",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(micGR,nanoGR,myn=250,to.add=T,my.col=mycols[2],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(htl,nanoGR,myn=250,to.add=T,my.col=mycols[3],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(weiGRC,nanoGR,myn=250,to.add=T,my.col=mycols[4],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(andGRC,nanoGR,myn=250,to.add=T,my.col=mycols[5],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
legend("topright",c("iclip","mic","ht","wei","anderson"),col=mycols[1:5],lty=1,lwd=1)
dev.off()

pdf("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/enrichment_of_sets_around_iCLIP.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(nanoGR,iCLIP[[3]],myn=250,to.add=F,my.col=mycols[1],my.ylim=c(0,0.011),my.title="Nanopore",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(micGR,iCLIP[[3]],myn=250,to.add=T,my.col=mycols[2],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(htl,iCLIP[[3]],myn=250,to.add=T,my.col=mycols[3],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(weiGRC,iCLIP[[3]],myn=250,to.add=T,my.col=mycols[4],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(andGRC,iCLIP[[3]],myn=250,to.add=T,my.col=mycols[5],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
legend("topright",c("nano","mic","ht","wei","anderson"),col=mycols[1:5],lty=1,lwd=1)
dev.off()

pdf("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/enrichment_of_sets_around_miclip.pdf",width=3.5,height=4)
par(mfrow=c(1,1))
makeRelPlot(nanoGR,micGR,myn=250,to.add=F,my.col=mycols[1],my.ylim=c(0,0.013),my.title="miclip",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(iCLIP[[3]],micGR,myn=250,to.add=T,my.col=mycols[2],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(htl,micGR,myn=250,to.add=T,my.col=mycols[3],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(weiGRC,micGR,myn=250,to.add=T,my.col=mycols[4],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
makeRelPlot(andGRC,micGR,myn=250,to.add=T,my.col=mycols[5],my.ylim=c(0,.013),my.title="Full set",to.smooth=F,use.para = TRUE,to.norm=T)
legend("topright",c("nano","iclip","ht","wei","anderson"),col=mycols[1:5],lty=1,lwd=1)
dev.off()


####################################################################################
######geneset enrichment plots######################################################
####################################################################################

dr_rand <- returnDens(randGR,gtfGR)


par(mfrow=c(1,1))
twGR <- iCLIP[[3]]
dr <- returnDens(twGR,gtfGR)
pdf(paste0("normalised_enrichments_plot_iCLIP.pdf"),width=6.5,height=6)
plotNE(bs=0.01,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=F,to.add=F,col.no=3,no.samps=25,to.smooth=T,mytitle = "iclip")
for(i in c(1,2,4))
{
  twGR <- iCLIP[[i]]
  dr <- returnDens(twGR,gtfGR)
  plotNE(bs=0.01,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=F,to.add=T,col.no=i,no.samps=25,to.smooth=T,mytitle = "iclip")
}
dev.off()



pdf(paste0("normalised_enrichments_plot_experiments.pdf"),width=6.5,height=6)
mybs <- 0.02
dr <- returnDens(iCLIP[[3]],gtfGR)
plotNE(bs=mybs,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=F,to.add=F,col.no=1,no.samps=25,to.smooth=T,mytitle = "iclip")

dr <- returnDens(micGR,gtfGR)
plotNE(bs=mybs,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=F,to.add=T,col.no=2,no.samps=25,to.smooth=T,mytitle = "iclip")

dr <- returnDens(htl,gtfGR)
plotNE(bs=mybs,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=F,to.add=T,col.no=3,no.samps=25,to.smooth=T,mytitle = "iclip")


dr <- returnDens(weiGRC,gtfGR)
plotNE(bs=mybs,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=F,to.add=T,col.no=4,no.samps=25,to.smooth=T,mytitle = "iclip")

dr <- returnDens(nanoGR,gtfGR)
plotNE(bs=mybs,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=F,to.add=T,col.no=5,no.samps=25,to.smooth=T,mytitle = "iclip")


dr <- returnDens(andGRC,gtfGR)
plotNE(bs=mybs,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=F,to.add=T,col.no=6,no.samps=25,to.smooth=T,mytitle = "iclip")


dr <- returnDens(m6AGR,gtfGR)
plotNE(bs=mybs,dr_rand=dr_rand,dr=dr,dr_ind=dr_ind,inc_groups=F,inc_int=F,to.add=T,col.no=7,no.samps=25,to.smooth=T,mytitle = "iclip")

legend("topleft",c("iclip","miclip","hyperTRIBE","wei","nanopore","anderson","shen"),col=brewer.pal("Set1",n = 9)[1:7],lwd=1)
dev.off()


