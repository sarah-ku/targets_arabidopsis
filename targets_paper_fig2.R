source("/home/sarah/github/targets_arabidopsis/targets_paper_open.R")
source("/home/sarah/github/targets_arabidopsis/ECT2_hyperTRIBE_model_data.R")

unlist(lapply(pos.list.ECT2,length))
unlist(lapply(pos.list.ECT3,length))

save_dir_r2 <-  "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/model_roots_single/saved_output/"
save_dir_s2 <-  "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/model_shoots_single/saved_output/"
save_dir_r3 <-  "/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/pipeline/model_roots_single/saved_output/"
save_dir_s3 <-  "/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/pipeline/model_shoots_single/saved_output/"

#overlap positions
#install.packages("eulerr")
library(eulerr)
#install.packages("GeneOverlap")
library(GeneOverlap)
library(systemPipeR)
library(gridExtra)
library(ggplot2)
library(Biostrings)
library(BSgenome.Athaliana.TAIR.TAIR9)
require(ggplot2)
require(ggseqlogo)
library(cluster)
setwd("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/")
getwd()

################################################################################################
#############ADAR EXPRESSION IN ROOTS AND SHOOTS########################
################################################################################################


###
adar <- tpm.mat["ADARclone",]
adar <- adar[grep("c|e",names(adar))]
treat <- rep("ECT2:ADAR",length(adar))
treat[grep("c",names(adar))] <- "WT"
my.df <- data.frame(names(adar),treat,adar)
colnames(my.df) <- c("sample","treat","adar")
my.df <- my.df[c(1:5,7:10,6,11:15,17:20,16),]
my.df$sample <- factor(as.vector(my.df$sample),levels=as.vector(my.df$sample))

p1 <- ggplot(my.df,aes(x=sample,y=adar,fill=treat),col=treat) + geom_bar(stat="identity") + theme_bw() +   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
  ggtitle("ADAR (TPM)")
p1
ggsave(filename="adar_expression_over_samples_roots_single.pdf",plot = p1,width = 6,height=4)
#t.test(adar[1:5],adar[6:10])


################################################################################################
#############VENN DIAGRAMS ROOTS AND SHOOTS EXPRESSED BOTH########################
################################################################################################

names(pos.list)
ht.shoots <- names(pos.list.ECT2[["shoots"]])
ht.roots <- names(pos.list.ECT2[["roots"]])

expressed.both <- intersect(genes.list$genes.shoots.expressed,genes.list$genes.roots.expressed)
genes.shoots.expressed <- genes.list$genes.shoots[genes.list$genes.shoots %in% expressed.both]
genes.roots.expressed <- genes.list$genes.roots[genes.list$genes.roots %in% expressed.both]

length(intersect(genes.shoots.expressed,genes.roots.expressed))

p1 <- makeEulerDiagram(list("ht.shoots"=genes.shoots.expressed, "ht.roots"=genes.roots.expressed), 
                 color_palette = c("#283148", "#913535", "#bbbbbb"),
                 plot_name = "overlaps_venn_all",
                 plot_title = "")

p2 <- makeEulerDiagram(list("ht.shoots"=ht.shoots, "ht.roots"=ht.roots), 
                 color_palette = c("#283148", "#913535", "#bbbbbb"),
                 plot_name = "overlaps_venn_all",
                 plot_title = "") 

p1
p3 <- grid.arrange(p1,p2,ncol=2)
p3
ggsave(p3,filename=paste0("HT_single_overlaps_expressed_common.pdf"),width=10,height=5)


################################################################################################
#############hyperTRIBE sequence preferences####################################################
################################################################################################

newseqs <- c("1","2","3","4","5","Mt","Pt")
names(newseqs) <- seqnames(seqinfo(Athaliana))
Athaliana <- renameSeqlevels(Athaliana, newseqs)
seqnames(seqinfo(Athaliana))

posGRseq <- pos.list.ECT2[["roots"]]
mcols(posGRseq) <- NULL
posGRseq <- posGRseq[as.vector(strand(posGRseq)) %in% c("+","-")]
to.add <- 2
my.seqs <- (getSeq(Athaliana, posGRseq+to.add))
my.mat <- as.matrix(my.seqs)
center <- to.add+1
apply(my.mat,2,table)

my.seqs <- my.seqs[my.mat[,center]=="A"]
my.mat <- my.mat[my.mat[,center]=="A",]

my.seqs <- as.vector(my.seqs)
gsub("T","U",my.seqs)

p1 <- ggseqlogo(as.vector(my.seqs))
p1

project_id <- "ECT2_roots"
ggsave(plot = p1,file=paste0("ECT2_hyperTRIBE_sequence_",project_id,".pdf"),width=3,height=2)

################################################################################################
##DENSITY PLOT of of EDITING PROPORTIONS of ect2-1 ECT2-FLAG-ADAR in roots – low!! (shoots in EXT.2A)#
################################################################################################





################################################################################################
###SCATTER PLOT with editing proportions in roots vs. shoots of targets
################################################################################################


prop_scatter <- function(GR1,GR2,name1="shoots_single",name2="roots_single")
{
  common.pos <- intersect(names(GR1),names(GR2))
  props <- data.frame("r1"=GR1[common.pos]$prop,"r2"=GR2[common.pos]$prop)

  scatterPlot <- ggplot(props,aes(x=sqrt(r1), y=sqrt(r2))) + 
    geom_bin2d(bins=200) + xlim(0,1) + ylim(0,1) +
    scale_color_manual(values = c('#999999','#E69F00')) + 
    geom_abline(lty=2) + xlab(paste0("SQRT editing proportion ",name1)) + ylab(paste0("SQRT editing proportion ",name2)) +
    theme_classic()
  return(scatterPlot)
}

scatterPlot <- prop_scatter(pos.list.ECT2[["roots"]],pos.list.ECT2[["shoots"]],"roots_single_ECT2","shoots_single_ECT2")
ggsave(plot=scatterPlot,file="ECT2_hyperTRIBE_roots_vs_shoots_single_proportions.pdf",width=4,height=3.2)

scatterPlot <- prop_scatter(pos.list.ECT3[["roots"]],pos.list.ECT3[["shoots"]],"roots_single_ECT3","shoots_single_ECT3")
ggsave(plot=scatterPlot,file="ECT3_hyperTRIBE_roots_vs_shoots_single_proportions.pdf",width=4,height=3.2)

scatterPlot <- prop_scatter(pos.list.ECT2[["roots"]],pos.list.ECT3[["roots"]],"roots_single_ECT2","roots_single_ECT3")
ggsave(plot=scatterPlot,file="hyperTRIBE_roots_ECT2_vs_roots_ECT3_single_proportions.pdf",width=4,height=3.2)

scatterPlot <- prop_scatter(pos.list.ECT2[["shoots"]],pos.list.ECT3[["shoots"]],"shoots_single_ECT2","shoots_single_ECT3")
ggsave(plot=scatterPlot,file="hyperTRIBE_shoots_ECT2_vs_shoots_ECT3_single_proportions.pdf",width=4,height=3.2)


################################################################################################
###SCATTER PLOT with editing proportions in roots vs. shoots of targets
################################################################################################


prop_scatter <- function(GR1,name="ECT2_root")
{
  #common.pos <- intersect(names(GR1),names(GR2))
  props <- data.frame("r1"=GR1$prop,"r2"=GR1$prop_ctrl)
  
  scatterPlot <- ggplot(props,aes(x=sqrt(r1), y=sqrt(r2))) + 
    geom_bin2d(bins=200) + xlim(0,1) + ylim(0,1) +
    scale_color_manual(values = c('#999999','#E69F00')) + 
    geom_abline(lty=2) + xlab(paste0("SQRT editing proportion in fusion ",name)) + ylab(paste0("SQRT editing proportion in WT ",name)) +
    theme_classic()
  return(scatterPlot)
}

scatterPlot <- prop_scatter(pos.list.ECT2[["roots"]],"ECT2_root")
ggsave(plot=scatterPlot,file="hyperTRIBE_roots_ECT2_WT_vs_fusion_proportions.pdf",width=4,height=3.2)

scatterPlot <- prop_scatter(pos.list.ECT2[["shoots"]],"ECT2_shoots")
ggsave(plot=scatterPlot,file="hyperTRIBE_shoots_ECT2_WT_vs_fusion_proportions.pdf",width=4,height=3.2)

scatterPlot <- prop_scatter(pos.list.ECT3[["roots"]],"ECT3_roots")
ggsave(plot=scatterPlot,file="hyperTRIBE_roots_ECT3_WT_vs_fusion_proportions.pdf",width=4,height=3.2)

scatterPlot <- prop_scatter(pos.list.ECT3[["shoots"]],"ECT3_shoots")
ggsave(plot=scatterPlot,file="hyperTRIBE_shoots_ECT3_WT_vs_fusion_proportions.pdf",width=4,height=3.2)


#################################################################################
#################################################################################
#PCA of significant edit sites
#################################################################################
#################################################################################

#data_list_roots
my_edits <- rbind(c("A","G"),
                  c("T","C"))


posGR <- pos.list.ECT2[["roots"]]

data_list <- data_list_roots[names(design_vector_roots_single)]
data_list <- lapply(data_list,function(x) x[names(posGR),])

edit_matrix_roots <- makeEditPCA(data_list=data_list,
                           editTypes=my_edits,
                           refGR=locsGR,
                           design_vector=design_vector_roots_single,refBased = F,posGR = posGR,give.mat = T)


posGR <- pos.list.ECT2[["shoots"]]

data_list <- data_list_shoots[names(design_vector_shoots_single)]
data_list <- lapply(data_list,function(x) x[names(posGR),])

edit_matrix_shoots <- makeEditPCA(data_list=data_list,
                                 editTypes=my_edits,
                                 refGR=locsGR,
                                 design_vector=design_vector_shoots_single,refBased = F,posGR = posGR,give.mat = T)

p1 <- grid.arrange(edit_matrix_roots[[1]],edit_matrix_shoots[[1]],ncol=2)
p1
ggsave(plot=p1,file="PCA_significant_edit_sites_hyperTRIBE.pdf",width=9,height=4)


#################################################################################
#################################################################################
#Correlation between editing proportions and ADAR levels in roots (shoots in EXT.2B)
#################################################################################
#################################################################################

edit_m <- edit_matrix_roots[[2]]

ADAR <- tpm.mat[1,colnames(edit_m)]
ADARt <- ADAR[grep("e",colnames(edit_m))]
edit_m <- edit_m[,grep("e",colnames(edit_m))]
head(edit_m)

edit.cor <- apply(edit_m,1,function(x) cor(x,ADARt))
edit.cor.bg <- apply(edit_m[,sample(1:4,4)],1,function(x) cor(x,ADARt[sample(1:4,4)]))
  
par(mfrow=c(1,2))
hist(edit.cor,breaks=50,col="blue")
hist(edit.cor.bg,col="grey",add=T,breaks=50)

hist(edit.cor.ctrl,breaks=50,col="blue")
hist(edit.cor.bg.ctrl,col="grey",add=T,breaks=50)


library(easyGgplot2)
df <- melt(data.frame("edit.cor"=edit.cor,"edit.cor.bg"=edit.cor.bg))
df$type <- "mut"
df.ctrl <-  melt(data.frame("edit.cor"=edit.cor.ctrl,"edit.cor.bg"=edit.cor.bg.ctrl))
df.ctrl$type <- "wt"
df <- rbind(df,df.ctrl)
head(df)

p1 <- ggplot(df, aes(x=value,color=variable))+
  #geom_histogram(color="darkblue", fill=c("lightblue"))
  geom_histogram(fill=c("white"), bins=30,alpha=0.5, position="identity") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +facet_wrap(~type) +theme_classic()

p1
ggsave(p1,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/correlation_of_editing_with_ADAR_SHOOTS.pdf",width=4.5,height=2.5)



################################################################################################
###Box plots with hits/p-value/editing proportions vs. expression (no expression-bias for p-values! )
################################################################################################

pl <- pos.list.ECT2[["roots"]]
expressed_genes_roots <- genes.list$genes.roots.expressed

mycuts <- cut2(log2(c(expr.list_ECT2$roots.expr,expr.list_ECT2$shoots.expr)+1),g=9,onlycuts = T)

g10 <- as.numeric(cut2(log(expr.list_ECT2$roots.expr+1),cuts = mycuts))
names(g10) <- names(expr.list_ECT2$roots.expr)
table(g10)

min.p <- tapply(pl$padj,pl$gene,function(x) -log10(min(x,na.rm=T)))
n.hits <- tapply(pl$prop,pl$gene,function(x) length(x))
max.prop <- tapply(pl$prop,pl$gene,function(x) mean(x))


mp <- rep(0,length(expressed_genes_roots))
ge <- rep(0,length(expressed_genes_roots))
nh <- rep(0,length(expressed_genes_roots))
minp <- rep(0,length(expressed_genes_roots))
names(mp) <- expressed_genes_roots
names(ge) <- expressed_genes_roots
names(nh) <- expressed_genes_roots
names(minp) <- expressed_genes_roots


mp[Reduce(intersect,list(names(expr.gene),names(max.prop),names(mp)))] <- max.prop[Reduce(intersect,list(names(expr.gene),names(max.prop),names(mp)))]
ge[intersect(names(ge),names(expr.gene))] <- expr.gene[intersect(names(ge),names(expr.gene))]
minp[Reduce(intersect,list(names(expr.gene),names(min.p),names(minp)))] <- min.p[Reduce(intersect,list(names(expr.gene),names(min.p),names(minp)))]
nh[Reduce(intersect,list(names(expr.gene),names(nh),names(n.hits)))] <- n.hits[Reduce(intersect,list(names(expr.gene),names(nh),names(n.hits)))]

mp
length(ge)
mpgdat <- data.frame("expr"=log2(ge+1),"max.prop"=mp,"min.p"=minp,"n.hits"=nh,"g10"=g10[row.names(mpgdat)])
head(mpgdat)
mpgdat <- mpgdat[(apply(mpgdat,1,function(x) sum(is.infinite(x))==0)),]

table(mpgdat$g10)
boxplot(mpgdat$min.p~mpgdat$g10)
boxplot(mpgdat$n.hits~mpgdat$g10)


#mpgdat$expr <- (cut(log(mpgdat$expr+1),breaks = c(-1,0,1,3,5,10,100000)))
#mpgdat$expr <- mpgdat$expr
#mpgdat$max.prop <- (cut(as.vector(mpgdat$max.prop),breaks = c(-1,0,0.01,0.05,0.1,0.25,0.5,1)))
mpgdat$max.prop <- (cut(as.vector(mpgdat$max.prop),20))

mpgdat$min.p <- (cut(mpgdat$min.p,breaks = c(-1,0,2,3,5,7,10,50000)))
mpgdat$n.hits <- (cut(mpgdat$n.hits,breaks = c(-1,0,1,2,5,10,100000)))

p1  <- ggplot(mpgdat,aes(y=log(expr+1), x=(n.hits))) + geom_boxplot(fill="#E69F00",outlier.size = 0.5) + theme_bw() + xlab("Number of hits in gene")
p1
p2  <- ggplot(mpgdat,aes(y=log(expr+1), x=(max.prop))) + geom_boxplot(fill="#E69F00",outlier.size = 0.5) + theme_bw() + xlab("Max editing proportion in gene")
p2
p3  <- ggplot(mpgdat,aes(y=log(expr+1), x=(min.p))) + geom_boxplot(fill="#E69F00",outlier.size = 0.5) + theme_classic() + xlab("Minimum p-value in gene")
p3

p4 <- grid.arrange(p1,p2,p3,ncol=3)

ggsave(p4,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/ECT2_hyperTRIBE_roots_expr_vs_editing_hits_pval.pdf",width=6.4,height=3.33)




#################################################################################
#################################################################################
#are roots specific genes more expressed in roots?
#################################################################################
#################################################################################

expressed.both <- intersect(genes.list$genes.shoots.expressed,genes.list$genes.roots.expressed)
genes.shoots.expressed <- genes.list$genes.shoots[genes.list$genes.shoots %in% expressed.both]
genes.roots.expressed <- genes.list$genes.roots[genes.list$genes.roots %in% expressed.both]

genes.shoots.expressed.specific <- genes.shoots.expressed[!(genes.shoots.expressed %in% genes.roots.expressed)]
shoots.spec <- tpm.mat.collapsed[genes.shoots.expressed.specific,]

genes.roots.expressed.specific <- genes.roots.expressed[!(genes.roots.expressed %in% genes.shoots.expressed)]
roots.spec <- tpm.mat.collapsed[genes.roots.expressed.specific,]

genes.shared <- intersect(genes.roots.expressed,genes.shoots.expressed)
shared.tpm <- tpm.mat.collapsed[genes.shared,]

genes.non <- intersect(genes.list$genes.nontargets.roots,genes.list$genes.nontargets)
genes.non <- genes.non[(genes.non %in% genes.list$genes.roots.expressed) & (genes.non %in% genes.list$genes.shoots.expressed)]
none.tpm <- tpm.mat.collapsed[genes.non,]

fc <- log(apply(shoots.spec[,grep("Sc",colnames(shoots.spec))],1,mean)/apply(shoots.spec[,grep("Rc",colnames(shoots.spec))],1,mean))
plot(density(fc))
abline(v=0,lty=2)

lgss.shoots <- apply(shoots.spec[,grep("Sc",colnames(shoots.spec))],1,mean)
lgss.roots <- apply(shoots.spec[,grep("Rc",colnames(shoots.spec))],1,mean)
lgrs.shoots <- apply(roots.spec[,grep("Sc",colnames(roots.spec))],1,mean)
lgrs.roots <- apply(roots.spec[,grep("Rc",colnames(roots.spec))],1,mean)
lgbt.shoots <- apply(shared.tpm[,grep("Sc",colnames(shared.tpm))],1,mean)
lgbt.roots <- apply(shared.tpm[,grep("Rc",colnames(shared.tpm))],1,mean)
lgnn.shoots <- apply(none.tpm[,grep("Sc",colnames(none.tpm))],1,mean)
lgnn.roots <- apply(none.tpm[,grep("Rc",colnames(none.tpm))],1,mean)

sspec_shoots <- data.frame("name"="Shoots_spec","type"="Shoots","expr"=lgss.shoots)
sspec_roots <- data.frame("name"="Shoots_spec","type"="Roots","expr"=lgss.roots)
rspec_shoots <- data.frame("name"="Roots_spec","type"="Shoots","expr"=lgrs.shoots)
rspec_roots <- data.frame("name"="Roots_spec","type"="Roots","expr"=lgrs.roots)

type <- c(rep("shoots_spec",length(lgss.shoots)),rep("roots_spec",length(lgrs.roots)))
scatdf <- data.frame( "type"=type,
  "shoots"=c(lgss.shoots,lgrs.shoots),
  "roots"=c(lgss.roots,lgrs.roots))
scatdf$shoots_log <- log(scatdf$shoots+1)
scatdf$roots_log <- log(scatdf$roots+1)
scatdf$fc <- 0
scatdf[scatdf$type=="roots_spec",]$fc <- log2((scatdf$shoots/scatdf$roots)[scatdf$type=="roots_spec"])
scatdf[scatdf$type=="shoots_spec",]$fc <- log2((scatdf$roots/scatdf$shoots)[scatdf$type=="shoots_spec"])

myfg <- rownames(scatdf[which(scatdf$fc>1.5),])

rsids <- row.names(scatdf[which(scatdf$fc>1.5 & scatdf$type=="roots_spec"),])
gdrsids <- gene.desc[rsids,]
gdrsids$fc <- round(as.vector(scatdf[which(scatdf$fc>1.5 & scatdf$type=="roots_spec"),]$fc),3)
gdrsids$iclip_supp <- gdrsids$V1 %in% iCLIP$E2_110kGR$genes[,1]
gdrsids$nano_supp <- gdrsids$V1 %in% nanoGR$genes[,1]
gdrsids$micseq <- gdrsids$V1 %in% micGR$genes[,1]
gdrsids$desc <- as.vector(gdrsids$V3)
gdrsids <- gdrsids[,-c(3:4)]

write.table(gdrsids,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/roots_specific_photosynthesis_genes.txt",quote=F,col.names = T,sep="\t",row.names = F)

scatdf[myfg,]
mybg <- genes.list$genes.expressed

log2(rowMeans(prot_expr[myfg,-c(1,3,5,7,9,11)])/rowMeans(prot_expr[myfg,c(1,3,5,7,9,11)]))
rowMeans(tpm.mat.collapsed[myfg,grep("Rc",colnames(tpm.mat.collapsed))])
rowMeans(tpm.mat.collapsed[myfg,grep("Re",colnames(tpm.mat.collapsed))])

rowMeans(tpm.mat.collapsed[myfg,grep("Sc",colnames(tpm.mat.collapsed))])
rowMeans(tpm.mat.collapsed[myfg,grep("Se",colnames(tpm.mat.collapsed))])

target.set.list$roots[myfg,]
target.set.list$union[myfg,]
par(mfrow=c(2,2))
par(mar=c(4,4,4,4))
plot(de.roots[myfg,]$new.p.levels,de.shoots[myfg,]$new.p.levels,ylim=c(0,8),pch=16,xlim=c(0,8))
abline(0,1)
plot(de.roots[myfg,]$old.p,de.shoots[myfg,]$old.p,ylim=c(0,8),pch=16,xlim=c(0,8))
abline(0,1)
plot(de.roots[myfg,]$new.p.delay,de.shoots[myfg,]$new.p.delay,ylim=c(0,8),pch=16,xlim=c(0,8))
abline(0,1)
plot(de.roots[myfg,]$new.p,de.shoots[myfg,]$new.p,ylim=c(0,8),pch=16,xlim=c(0,8))
abline(0,1)

myfg <- names(rowSums(gdrsids[,6:7])>0)

library(gProfileR)
mysc <- gprofiler(myfg, organism = "athaliana", ordered_query = F,
                  region_query = F, max_p_value = 0.05, min_set_size = 0,
                  min_isect_size = 0, correction_method = "analytical",custom_bg = mybg,
                  hier_filtering = "none", domain_size = "annotated",
                  numeric_ns = "", include_graph = F, src_filter = NULL)

mysc[order(mysc$p.value),]

mydf <- rbind(sspec_shoots,sspec_roots,rspec_shoots,rspec_roots)
head(mydf)


mymin <- min(c(lgss.shoots,lgss.roots,lgrs.shoots,lgrs.roots))
mymax <- max(c(lgss.shoots,lgss.roots,lgrs.shoots,lgrs.roots,lgbt.roots,lgbt.shoots))

par(mfrow=c(1,1))
plot(lgss.roots,lgss.shoots,ylim=c(0,mymax),xlim=c(0,mymax),cex=0.2,xlab="Expression in roots (log(tpm+1))",ylab="Expression in shoots (log(tpm+1))")
points(lgrs.roots,lgrs.shoots,col="red",cex=0.2)
#points(lgss.roots,lgss.shoots,col="blue",cex=0.2)
abline(0,1,lty=2)

plot(lgnn.roots,lgnn.shoots,ylim=c(0,mymax),xlim=c(0,mymax),cex=0.2,xlab="Expression in roots (log(tpm+1))",ylab="Expression in shoots (log(tpm+1))")
abline(0,1,lty=2)

#df <- data.frameshoots"=lgss.shoots)
scatterPlot <- ggplot(scatdf,aes(x=shoots,y=roots,fill=type)) + 
  geom_bin2d(bins=100) +
  scale_color_manual(values = c('#999999','#E69F00')) + 
  geom_abline(lty=2) + xlab("Expression in shoots") + ylab("Expression in roots") +
  theme_classic()
scatterPlot
ggsave(scatterPlot,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/scatter_cell_specific_target_expression.pdf",width=7,height=6)

plot(lgbt.roots,lgbt.shoots,ylim=c(0,mymax),xlim=c(0,mymax),cex=0.2,xlab="Expression in roots (log(tpm+1))",ylab="Expression in shoots (log(tpm+1))")
abline(0,1,lty=2)

points(lgbt.roots,lgbt.shoots,col="red",cex=0.2)

apply(shoots.spec[,grep("Sc",colnames(shoots.spec))],1,sd)

sd(rowMeans(shoots.spec[,grep("Sc",colnames(shoots.spec))]))

mean(rowMeans(shoots.spec[,grep("Rc",colnames(shoots.spec))]))
sd(rowMeans(shoots.spec[,grep("Rc",colnames(shoots.spec))]))

mean(rowMeans(roots.spec[,grep("Rc",colnames(shoots.spec))]))
mean(rowMeans(roots.spec[,grep("Sc",colnames(shoots.spec))]))

mean(rowMeans(shared.tpm[,grep("Rc",colnames(shoots.spec))]))
mean(rowMeans(shared.tpm[,grep("Sc",colnames(shoots.spec))]))

mydf.list <- list()

gene <- "AT1G61520"
for(gene in myfg)
{
  my.m <- tapply((tpm.mat.collapsed[gene,]),rep(1:6,each=5),mean)
  my.sd <- tapply((tpm.mat.collapsed[gene,]),rep(1:6,each=5),sd)
  my.m.prot <- tapply((prot_expr[gene,]),rep(c(1,2),6),mean)
  my.sd.prot <- tapply((prot_expr[gene,]),rep(c(1,2),6),sd)
  my.sd <- c(my.sd,my.sd.prot)
  my.m <- c(my.m,my.m.prot)
  
  name <- c("root_c","root_s","root_t","shoot_c","shoot_s","shoot_t","prot_wt","prot_mut")
  type <- c(rep("root",3),rep("shoot",3),rep("prot",2))
  mydf.list[[gene]] <- data.frame("gene"=gene,"type"=type,"name"=name,"mean"=my.m,"sd"=my.sd)
  
}

mydf <- do.call(rbind,mydf.list)
head(mydf)
mydf <- mydf[-grep("_t",mydf$name),]

p<- ggplot(mydf, aes(x=name, y=mean, fill=type)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) + theme_classic() +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) + facet_wrap(~gene, scales="free_y")

p
ggsave(p,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/pipeline/plots/genes_specific_expression_barplots.pdf",width=8,height=8)


##############PROTOPLAST SPECIFICITY AND HYPERTIRBE editing proportion
heg <- pos.list[[1]][which(pos.list[[1]]$meta$prop>0.25)]$genes[,1]
leg <- pos.list[[1]][which(pos.list[[1]]$meta$prop<0.01)]$genes[,1]

cgenes <- Reduce(intersect,list(heg,row.names(protm1),row.names(tpm.mat.collapsed)))

prote <- rowMeans(protm1[cgenes,])
hte <- rowMeans(tpm.mat.collapsed[cgenes,grep("Rc",colnames(tpm.mat.collapsed))])

cgenes <- Reduce(intersect,list(leg,row.names(protm1),row.names(tpm.mat.collapsed)))

proteL <- rowMeans(protm1[cgenes,])
hteL <- rowMeans(tpm.mat.collapsed[cgenes,grep("Rc",colnames(tpm.mat.collapsed))])

hist(log(proteL/hteL),pch=16,col="blue")
abline(v=0,col="black",lwd=3)
hist(log(prote/hte),col="red",pch=16,add=T)
abline(0,1)

cgenes <- Reduce(intersect,list(row.names(protm1),row.names(tpm.mat.collapsed)))

lgfc <- log(rowMeans(protm1[cgenes,])/rowMeans(tpm.mat.collapsed[cgenes,grep("Rc",colnames(tpm.mat.collapsed))]))
mprop <- tapply(pos.list[[1]]$meta$prop,pos.list[[1]]$genes[,1],max)
lgfcm <- lgfc[intersect(names(mprop),names(lgfc))]
mpropm <- mprop[intersect(names(mprop),names(lgfc))]
plot(lgfcm,mpropm)
abline(v=0,col="red",lty=2,lwd=3)
summary(glm(factor(mpropm[!is.infinite(lgfcm)]>0.25)~lgfcm[!is.infinite(lgfcm)],family=binomial(link="logit")))

#no eveidence that protoplast specific genes are more highly edited

mean(mpropm[lgfcm>2])
mean(mpropm[lgfcm<(-2)])



heg <- unique(pos.list[[1]][which(pos.list[[1]]$meta$prop>0.25)]$genes[,1])
heg
leg <- unique(pos.list[[1]][which(pos.list[[1]]$meta$prop<0.01)]$genes[,1])
bg <- unique(pos.list[[1]]$genes[,1])

#fg <- names(which(lgfcm>1))
#bg <- names(lgfcm)

library(gprofiler2)
mysc <- gprofiler(heg, organism = "athaliana", ordered_query = F,
                  region_query = F, max_p_value = 0.05, min_set_size = 0,
                  min_isect_size = 0, correction_method = "analytical",custom_bg = bg,
                  hier_filtering = "none", domain_size = "annotated",
                  numeric_ns = "", include_graph = F, src_filter = NULL)

mysc[order(mysc$p.value),-c(14)]

#edit_matrix[names(pos.list[[2]][which(pos.list[[2]]$meta$prop>0.25)]),]


##correlaton with ADAR
###

# Change line type

colMeans(edit_matrix[,6:10])
apply(edit_matrix[,6:10],2,sd)


#######highly edited gene
common.pos <- intersect(names(pos.list[[1]]),names(pos.list[[2]]))
prop.roots <- pos.list[[1]][common.pos]$meta$prop
prop.shoots <- pos.list[[2]][common.pos]$meta$prop
lfc.p <- (log(prop.roots/prop.shoots))

expr.roots <- as.numeric(unlist(lapply(strsplit(pos.list[[1]][common.pos]$genes[,4],","),function(x) x[1])))
expr.shoots <- as.numeric(unlist(lapply(strsplit(pos.list[[2]][common.pos]$genes[,4],","),function(x) x[1])))
lfc.e <- log(expr.roots/expr.shoots)

plot(lfc.p,lfc.e,col=densCols(lfc.p,lfc.e))
abline(h=0)
abline(v=0)


#does it tend to be the same genes?
sort(table(pos.list[[1]][pos.list[[1]]$meta$prop>0.25]$genes[,1]))
k <- k+1
pv <- pos.list[[1]][pos.list[[1]]$genes[,1]==hge[k],]$meta$prop
plot(100*pv,type="h")

ca <- queryHits(findOverlaps(pos.list[[1]][pos.list[[1]]$genes[,1]==hge[k]],coldGR))
points(100*pos.list[[1]][pos.list[[1]]$genes[,1]==hge[k],]$meta$prop[ca]~ca,type="h",col="red")

#ca <- queryHits(findOverlaps(pos.list[[1]][pos.list[[1]]$genes[,1]==hge[k]],hotGR))

summary(pos.list[[1]][]$meta$prop)

#summary(pos.list[[1]][queryHits(findOverlaps(pos.list[[1]][pos.list[[1]]$genes[,1]==hge[k]],hotGR))]$meta$prop)

summary(pos.list[[1]][queryHits(findOverlaps(pos.list[[1]]+1000,hotGR))]$meta$prop)



m.prop <- tapply(pos.list[[2]]$meta$prop,pos.list[[2]]$genes[,1],max)

hge <- names(m.prop[m.prop>0.5])
lge <- names(m.prop[m.prop<0.01])
hge.st <- structure[grep(paste(hge,collapse="|"),structure[,1]),]
summary(hge.st)

lge.st <- structure[grep(paste(lge,collapse="|"),structure[,1]),]
summary(lge.st)

t.test(hge.st[,2],lge.st[,2])

reac_shoot <- read.table("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/structure/shoot_control_reactivity",header=T)
summary(reac_shoot[grep(paste(hge,collapse="|"),reac_shoot[,1]),])
summary(reac_shoot[grep(paste(lge,collapse="|"),reac_shoot[,1]),])

htGR <- pos.list[[2]]
comGR <- htGR
pre <- follow(htGR,comGR)
dtn <- abs(start(htGR[!is.na(pre)])-start(comGR[pre[!is.na(pre)]]))
prop <-  htGR$meta$pro[!is.na(pre)]


prop <- prop[dtn<25000]
dtn <- dtn[dtn<25000]
par(mfrow=c(1,1))
boxplot(sqrt(prop)~cut2(log2(dtn+1),g=100),cex=0.5,pch=16,outline=F)
#ss <- smooth.spline(tapply(sqrt(prop),cut2(log2(dtn+1),g=100),median)~1:100)$y
#points(ss,type="l",lwd=5,col="red")


par(mfrow=c(1,2))
htGR <- iCLIP[[3]]
comGR <- nanoGR
dtn <- as.data.frame(distanceToNearest(htGR,comGR))[,3]
hist(log(dtn+1),breaks=50)

htGR <- iCLIP[[4]]
comGR <- nanoGR
dtn <- as.data.frame(distanceToNearest(htGR,comGR))[,3]
hist(log(dtn+1),breaks=50)


htGR <- iCLIP[[3]]
nc <- (as.data.frame(distanceToNearest(htGR,coldGR))[,3])
nh <- (as.data.frame(distanceToNearest(htGR,hotGR))[,3])

hist(log(nc))
hist(log(nh),add=T,col="red")

df <- melt(data.frame("nc"=nc,"nh"=nh))

p1 <- ggplot(df, aes(x=log(value+1),color=variable))+
  #geom_histogram(color="darkblue", fill=c("lightblue"))
  geom_histogram(fill=c("pink"), bins=30,alpha=0.5, position="identity") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))  +theme_classic()
p1













############################### single cell proportions ###################


scm <- read.csv2("/binf-isilon/alab/projects/ECT2_TC/protoplast_markers/denyer_markers/GSE123818_Root_single_cell_wt_datamatrix.csv",sep=",",header = T,row.names = 1)
scm <- as.matrix(scm)


pl <- pos.list.ECT2[["roots"]]
min.p <- tapply(pl$padj,pl$gene,function(x) -log10(min(x,na.rm=T)))
n.hits <- tapply(pl$prop,pl$gene,function(x) length(x))
max.prop <- tapply(pl$prop,pl$gene,function(x) mean(x))

head(scm)
ect2 <- "AT3G13460"
ect3 <- "AT5G61020"
ect2_prop <- scm[ect2,]/colSums(scm)

ect2 <- sample(row.names(scm),1)

ect2_plus <- as.vector(scm[ect2,])>0
ect2_minus <- as.vector(scm[ect2,])==0
  
psrat <- rep(0,nrow(scm))
names(psrat) <- row.names(scm)
for(i in row.names(scm))
{
  ps <- sum(scm[i,ect2_plus])
  ms <- sum(scm[i,ect2_minus])
  psrat[i] <- ps/(ps+ms)
}
mep <- tapply(pos.list.ECT2[[1]]$prop,pos.list.ECT2[[1]]$gene,mean)
  
ci <- (psrat[names(mep)])
co <- (psrat[!(names(psrat) %in% names(mep))])


etab <- rbind(table(c(scm[agene,ect2_minus]>0)),table(c(scm[agene,ect2_plus]>0)))
fisher.test(etab)$estimate
summary(scm[agene,ect2_minus]>0)

mean(pos.list.ECT2[[1]][grep(agene,pos.list.ECT2[[1]]$gene)]$prop)
mean(as.numeric())

scm2 <- scm[row.names(scm) %in% unique(pl$gene),]

cos <- rep(0,nrow(scm))
mi <- rep(0,nrow(scm))
mcor <- rep(0,nrow(scm))
mi_rand <- rep(0,nrow(scm))
for(i in 1:nrow(scm))
{
  #cos[i] <- sum(scm[i,] * scm[ect2,] > 0)/sum(scm[i,]>0)
  #mi[i] <- mutinformation(as.numeric(scm[i,]>0),as.numeric(scm[ect2,]>0), method="emp")
  mi_rand[i] <- mutinformation(as.numeric(scm[i,]>0),as.numeric(scm[sample(1:nrow(scm),1),]>0), method="emp")
  #mcor[i] <- cor(log2(scm[i,]+1),log2(scm[ect2,]+1))
}
names(cos) <- row.names(scm)
names(mi) <- row.names(scm)
names(mi_rand) <- row.names(scm)
names(mcor) <- row.names(scm)

plot(cos,mi)
points(cos[bound],mi[bound],col="red")
points(cos[unbound],mi[unbound],col="blue")


pdf("coexpression_ect2_and_mean_editing_proportion.pdf",width=5,height=5)
boxplot(sqrt(max.prop)~cut2(mi[names(max.prop)],g=100),xlab="% of cells co-expressing ECT2",outline=F,ylab="SQRT mean editing proportion in gene")
boxplot(log2(mi)~cut2(expr.list_ECT2[["roots.expr"]][names(mi)],g=100),xlab="% of cells co-expressing ECT2",outline=F,ylab="SQRT mean editing proportion in gene")
boxplot(expr.list_ECT2[["roots.expr"]][names(max.prop)]~cut2(mi[names(max.prop)],g=100),xlab="% of cells co-expressing ECT2",outline=F,ylab="SQRT mean editing proportion in gene")
dev.off()

unbound <- nanoGR$gene[!(nanoGR$gene %in% genes.list$genes.roots)]
bound <- nanoGR$gene[(nanoGR$gene %in% genes.list$genes.roots)]
other <- row.names(scm)[!(row.names(scm) %in% nanoGR$gene)]

mge <- c(rep("unbound",length(unbound)),rep("bound",length(bound)),rep("other",length(other)))
names(mge) <- c(unbound,bound,other)
library(RColorBrewer)

pdf("mutual_information_with_ECT2_bound_nano_vs_unbound.pdf",width=3.5,height=5)
par(mfrow=c(1,2))
boxplot(mi[names(mge)]~factor(mge,levels=c("bound","unbound","other")),ylim=c(0,0.07),outline=F,ylab="mutual information with ECT2",xlab="",col=brewer.pal(3,"Set1"),main="Denyer et al MI with ECT2")
boxplot(mi_rand[names(mge)]~factor(mge,levels=c("bound","unbound","other")),ylim=c(0,0.07),outline=F,ylab="mutual information with ECT2",xlab="",col=brewer.pal(3,"Set1"),main="Denyer et al MI with ECT2")
#mivr <- mi[names(mge)]/mi_rand[names(mge)]
#boxplot(mivr~factor(mge,levels=c("bound","unbound","other")),outline=F,ylab="mutual information with ECT2",xlab="",col=brewer.pal(3,"Set1"),main="Denyer et al MI with ECT2")

#boxplot(mi[names(mge)]~factor(mge,levels=c("bound","unbound","other")),outline=F,ylab="mutual information with ECT2",xlab="",col=brewer.pal(3,"Set1"),)
dev.off()

summary(cos[unbound[unbound %in% names(cos)]],na.rm=T)
summary(cos[bound[bound %in% names(cos)]],na.rm=T)

summary(mi[unbound[unbound %in% names(cos)]],na.rm=T)
summary(mi[bound[bound %in% names(cos)]],na.rm=T)


library(infotheo)
1-mutinformation(scm[ect3,],scm[ect2,], method="emp")

hist(1-cos,breaks=100)

