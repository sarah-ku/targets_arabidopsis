source("/home/sarah/github/targets_arabidopsis/targets_paper_open.R")

####################################################################################################
####################################################################################################
####################################################################################################
posGRF <- pos.list.ECT2
posGRF <- pos.list.ECT3

my.df.roots <- positionsToGenes(posGR=posGRF[["roots"]],include.ids = T,gene.ids=ids,use.meta = T)
my.df.shoots <- positionsToGenes(posGR=posGRF[["shoots"]],include.ids = T,gene.ids=ids,use.meta = T)
my.df.iclip <- positionsToGenes(posGR=iCLIP[[3]],include.ids = T,gene.ids=ids,use.meta = F)
my.df.m6a <- positionsToGenes(posGR=nanoGR,include.ids = T,gene.ids=ids,use.meta = F)

#my.df.list <- list("my.df.roots"=my.df.roots,"my.df.shoots"=my.df.shoots,"my.df.iclip"=my.df.iclip,my.df.m6a="my.df.m6a")
#######################################################
#######################################################
#######################################################

my.df.list <- lapply(posGRF[1:2],function(x) positionsToGenes(posGR=x,include.ids = T,gene.ids=ids,use.meta = T))
my.df.list <- lapply(my.df.list,function(x) x[!(x$gene %in% grep("_",x$gene,value=T)),])
my.df.pos <- lapply(posGRF[1:2],function(x) positionsToDf(posGR = x,include.ids = T,gene.ids=ids))

for(i in 1:2)
{
  my.df <- my.df.list[[i]]
  my.pos <- my.df.pos[[i]]
  project_id <- names(posGRF)[i]
  save_dir <- "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/"
  
  write.xlsx2(my.df, file=paste(save_dir,"hyperTRIBE_ECT2_hits_",project_id,"_all_genes.xlsx",sep=""), sheetName = "genes_ht_all",
              col.names = TRUE, row.names = FALSE, append = FALSE)
  
  
  write.xlsx(my.pos, file=paste(save_dir,"hyperTRIBE_ECT2_hits_",project_id,"_all_positions.xlsx",sep=""), sheetName = "positions_ht_all",
             col.names = TRUE, row.names = FALSE, append = FALSE)
}


posGRF <- pos.list.ECT3
my.df.list <- lapply(posGRF[1:2],function(x) positionsToGenes(posGR=x,include.ids = T,gene.ids=ids,use.meta = T))
my.df.list <- lapply(my.df.list,function(x) x[!(x$gene %in% grep("_",x$gene,value=T)),])
my.df.pos <- lapply(posGRF[1:2],function(x) positionsToDf(posGR = x,include.ids = T,gene.ids=ids))

for(i in 1:2)
{
  my.df <- my.df.list[[i]]
  my.pos <- my.df.pos[[i]]
  project_id <- names(posGRF)[i]
  save_dir <- "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/"
  
  write.xlsx2(my.df, file=paste(save_dir,"hyperTRIBE_ECT3_hits_",project_id,"_all_genes.xlsx",sep=""), sheetName = "genes_ht_all",
              col.names = TRUE, row.names = FALSE, append = FALSE)
  
  
  write.xlsx(my.pos, file=paste(save_dir,"hyperTRIBE_ECT3_hits_",project_id,"_all_positions.xlsx",sep=""), sheetName = "positions_ht_all",
             col.names = TRUE, row.names = FALSE, append = FALSE)
}

posGRF <- iCLIP
for(i in 1:4)
{
  #posGRF[[i]]$padj <- 1/posGRF[[i]]$score
  posGRF[[i]]$padj <- 1:length(posGRF[[i]])
  posGRF[[i]]$prop <- 1
  posGRF[[i]]$fold_change <- 1
  posGRF[[i]]$ref <- NA
  posGRF[[i]]$targ <- NA
  seqnames(Athaliana) <- c("1","2","3","4","5","Mt","Pt")
  posGRF[[i]]$ref[which(!(as.vector(strand(posGRF[[i]]))=="*"))] <- as.vector(getSeq(Athaliana,posGRF[[i]][-which(as.vector(strand(posGRF[[i]]))=="*")]))
  posGRF[[i]]$targ[which(!(as.vector(strand(posGRF[[i]]))=="*"))] <- as.vector(getSeq(Athaliana,posGRF[[i]][-which(as.vector(strand(posGRF[[i]]))=="*")]))
  posGRF[[i]]$transcript.types <- gsub("exon","CDS",posGRF[[i]]$transcript.types)
}
lapply(posGRF,length)

my.df.list <- lapply(posGRF,function(x) positionsToGenes(posGR=x,include.ids = T,gene.ids=ids,use.meta = T))
my.df.list <- lapply(my.df.list,function(x) x[!(x$gene %in% grep("_",x$gene,value=T)),])

#colnames(my.df.list[[1]])
for(i in 1:4)
{
  my.df.list[[i]]$padj <- unlist(lapply(base::strsplit(as.vector(my.df.list[[i]]$padj),","),function(x) paste(round(1/as.numeric(x),3),collapse=",")))
  my.df.list[[i]] <- my.df.list[[i]][,c("gene","name","strand","transcripts","transcripts.tpm","hits","ids")]
  #colnames(my.df.list[[i]])[ncol(my.df.list[[i]])] <- "score"
}


my.df.pos <- lapply(posGRF,function(x) positionsToDf(posGR = x,include.ids = T,gene.ids=ids))
head(my.df.pos[[1]])
for(i in 1:4)
{
  my.df.pos[[i]] <- my.df.pos[[i]][,c("id","gene","name","score","gtf_strand","transcripts","transcript.tpm","transcript.types","feature_prop","ref")]
  my.df.pos[[i]]$binding_site <- "no"
  my.df.pos[[i]]$binding_site[row.names(my.df.pos[[i]]) %in% names(iCLIPBS[[i]])] <- "yes"
}

for(i in 1:4)
{
  my.df <- my.df.list[[i]]
  my.pos <- my.df.pos[[i]]
  project_id <- names(posGRF)[i]
  save_dir <- "/binf-isilon/alab/projects/ECT2_TC/iCLIP/"
  
  write.xlsx2(my.df, file=paste(save_dir,"iCLIP_hits_",project_id,"_all_genes.xlsx",sep=""), sheetName = paste0(project_id,"_genes_all"),
              col.names = TRUE, row.names = FALSE, append = FALSE)
  
  
  write.xlsx(my.pos, file=paste(save_dir,"iCLIP_hits_",project_id,"_all_positions.xlsx",sep=""), sheetName = paste0(project_id,"_positions_all"),
             col.names = TRUE, row.names = FALSE, append = FALSE)
}



#quant.vec <- rowMeans(tpm.mat[,grep("c",colnames(tpm.mat))])

#######################################################
#######################################################
#######################################################


#"other" hits - all - genes
posGR.filtered <- posGR[!(paste(posGR$meta$ref,posGR$meta$targ,sep=":") %in% c("A:G","T:C"))]
posGR.filtered <- posGR.filtered[posGR.filtered$meta$padj<0.01]
posGR.filtered <- posGR.filtered[posGR.filtered$meta$prop<1]
my.df <- positionsToGenes(posGR.filtered,T,gene.ids=ids)
dim(my.df)

write.xlsx(my.df, file=paste(save_dir,"hyperTRIBE_hits_",project_id,"_other_all_genes.xlsx",sep=""), sheetName = "genes_other_all",
           col.names = TRUE, row.names = FALSE, append = FALSE)

#"other" hits - all - positions
gnames <- my.df$name
names(gnames) <- my.df$gene
my.df <- positionsToDf(posGR.filtered,T,gene.ids=gnames)
dim(my.df)

write.xlsx(my.df, file=paste(save_dir,"hyperTRIBE_hits_",project_id,"_other_all_positions.xlsx",sep=""), sheetName = "positions_other_all",
           col.names = TRUE, row.names = FALSE, append = FALSE)


posGR.filtered <- posGR[posGR$meta$prop==1]
my.df <- positionsToGenes(posGR.filtered,T,gene.ids=ids)
gnames <- my.df$name
names(gnames) <- my.df$gene
my.df <- positionsToDf(posGR.filtered,T,gene.ids=gnames)
write.xlsx(my.df, file=paste(save_dir,"hyperTRIBE_hits_",project_id,"_potential_SNPS.xlsx",sep=""), sheetName = "positions_snps",
           col.names = TRUE, row.names = FALSE, append = FALSE)

#######################################################
#######################################################
#######################################################


#make my.df.expr.roots and my.df.expr.shoots
set.list <- list("HT.roots"=my.df.roots,"HT.shoots"=my.df.shoots,"iCLIP"=my.df.iclip,"m6A"=my.df.m6a)
#set.list$iCLIP[set.list$iCLIP %in% set.list$my.df.roots]

shoots.expr <- rowMeans(tpm.mat.collapsed[,grep("Sc",colnames(tpm.mat.collapsed))])
roots.expr <- rowMeans(tpm.mat.collapsed[,grep("Rc",colnames(tpm.mat.collapsed))])
expr.list_ECT2 <- list("shoots.expr"=shoots.expr,"roots.expr"=roots.expr)


save(expr.list_ECT2,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/expr.list_ECT2.Rdat")

expr.list <- expr.list_ECT3
####################################################################################################
####################################################################################################
####################################################################################################

genes.shoots <- unique(row.names(my.df.shoots))
genes.roots <- unique(row.names(my.df.roots))
genes.single <- unique(c(row.names(my.df.roots),row.names(my.df.shoots)))
genes.union <- unique(c(row.names(my.df.roots),row.names(my.df.shoots),row.names(my.df.iclip)))
genes.intersect <- unique(intersect(genes.single,row.names(my.df.iclip)))
genes.roots.intersect <- unique(intersect(genes.roots,row.names(my.df.iclip)))
genes.shoots.expressed <- names(which(apply(tpm.mat.collapsed[,grep("Sc",colnames(tpm.mat.collapsed))],1,function(x) length(x[x>1]))>1))
genes.roots.expressed <- names(which(apply(tpm.mat.collapsed[,grep("Rc",colnames(tpm.mat.collapsed))],1,function(x) length(x[x>1]))>1))
genes.expressed <- unique(union(genes.shoots.expressed,genes.roots.expressed))
genes.nontargets <- genes.expressed[!(genes.expressed %in% genes.union)]
genes.nontargets.roots <- genes.expressed[!(genes.expressed %in% genes.roots)]

genes.shoots.intersect <- unique(intersect(genes.shoots,row.names(my.df.iclip)))
genes.nontargets.shoots <-  genes.expressed[!(genes.expressed %in% genes.shoots)]

genes.list <- list("genes.shoots"=genes.shoots,
     "genes.roots"=genes.roots,
     "genes.single"=genes.single,
     "genes.union"=genes.union,
     "genes.intersect"=genes.intersect,
     "genes.roots.intersect"=genes.roots.intersect,
     "genes.shoots.intersect"=genes.shoots.intersect,
     "genes.shoots.expressed"=genes.shoots.expressed,
     "genes.roots.expressed"=genes.roots.expressed,
     "genes.expressed"=genes.expressed,
     "genes.nontargets"=genes.nontargets,
     "genes.nontargets.roots"=genes.nontargets.roots,
     "genes.nontargets.shoots"=genes.nontargets.shoots)

unlist(lapply(genes.list,length))

genes.list <- lapply(genes.list,function(x) x[!(x %in% grep("_",x,value=T))])

save(genes.list,file = "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/genes.list_ECT3.Rdat")

length(genes.shoots.expressed)
length(genes.roots.expressed)
length(genes.roots)
length(genes.shoots)
length(genes.expressed)
length(genes.union)  
length(genes.intersect)
length(row.names(my.df.iclip))

df.expressed <- getSupportDf(genes=genes.expressed,set.list = set.list,expr.list=expr.list,gene.ids=ids)
df.union <- getSupportDf(genes=genes.union,set.list = set.list,expr.list=expr.list,gene.ids=ids)
df.roots.nontargets <- getSupportDf(genes=genes.roots.expressed[!(genes.roots.expressed %in% genes.union)],set.list = set.list,expr.list=expr.list,gene.ids=ids)
df.shoots.nontargets <- getSupportDf(genes=genes.shoots.expressed[!(genes.shoots.expressed %in% genes.union)],set.list = set.list,expr.list=expr.list,gene.ids=ids)
df.nontargets <- getSupportDf(genes=genes.expressed[!(genes.expressed %in% genes.union)],set.list = set.list,expr.list=expr.list,gene.ids=ids)

library(rJava)
library(xlsx)
options(java.parameters=c("-Xmx100g"))
save.dir <- "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/"

#ALL TISSUES (for comparison with m6A and FA-CLIP datasets)
#"STRINGENT TARGET SET ALL-TISSUES": Intersection between E2T-single[roots+shoots] and iCLIP110kDa
df.intersect_all <- getSupportDf(genes=genes.intersect,set.list = set.list,expr.list=expr.list,gene.ids=ids)

#"LARGE TARGET SET ALL-TISSUES": Union between E2T-single[roots+shoots] and iCLIP110kDa
df.union_all <- getSupportDf(genes=genes.union,set.list = set.list,expr.list=expr.list,gene.ids=ids)

#"NON-TARGETS ALL-TISSUES": Genes expressed in E2T-single[roots+shoots] that are NOT in the "LARGE TARGET SET ALL TISSUES".
df.nontargets_all <- getSupportDf(genes=genes.nontargets,set.list = set.list,expr.list=expr.list,gene.ids=ids)

write.xlsx2(df.intersect_all, file=paste(save.dir,"ALL_TISSUES_TARGETS_stringent.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

write.xlsx2(df.union_all, file=paste(save.dir,"ALL_TISSUES_TARGETS_permissive.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

write.xlsx2(df.nontargets_all, file=paste(save.dir,"ALL_TISSUES_TARGETS_non_targets.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

#ONLY ROOTS (for protoplast-DEA and Austria-analyses)
#"STRINGENT TARGET SET ONLY-ROOTS": Intersection between E2T-single[roots] and iCLIP110kDa
df.intersect_roots <- getSupportDf(genes=genes.roots.intersect,set.list = set.list,expr.list=expr.list,gene.ids=ids)
dim(df.intersect_roots)

write.xlsx2(df.intersect_roots, file=paste(save.dir,"ROOTS_ONLY_TARGETS_stringent.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

#"LARGE TARGET SET ONLY-ROOTS": E2T-single[roots]
df_roots <- getSupportDf(genes=genes.roots,set.list = set.list,expr.list=expr.list,gene.ids=ids)
dim(df_roots)

write.xlsx2(df_roots, file=paste(save.dir,"ROOTS_ONLY_TARGETS_permissive.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

#"NON-TARGETS ONLY-ROOTS": Genes expressed in E2T-single[roots] that are NOT in the "LARGE TARGET SET ONLY-ROOTS".
df.non_targets_roots <- getSupportDf(genes=genes.nontargets.roots,set.list = set.list,expr.list=expr.list,gene.ids=ids)
dim(df.non_targets_roots)

write.xlsx2(df.non_targets_roots, file=paste(save.dir,"ROOTS_ONLY_TARGETS_non_targets.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

target.set.list <- list("union"=df.union_all,"intersect"=df.intersect_all,"non_targets_all"=df.nontargets_all,"roots"=df_roots,"roots.intersect"=df.intersect_roots,"non_targets_roots"=df.non_targets_roots)
save.dir <- "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/"
save(target.set.list,file=paste0(save.dir,"target.set.list_new.Rdat"))


getSupportDf <- function(genes,set.list,expr.list,gene.ids=ids)
{
  #remove strays
  if(length(grep("_[0-9]+",genes))>0)
  {
    genes <- genes[-grep("_[0-9]+",genes)]
  }
  
  supp <- matrix(0,ncol=length(set.list),nrow=length(genes))
  row.names(supp) <- genes
  colnames(supp) <- names(set.list)
  
  for(k in 1:length(set.list))
  {
    supp[intersect(genes,row.names(set.list[[k]])),names(set.list)[k]] <- as.vector(set.list[[k]][intersect(genes,row.names(set.list[[k]])),]$hits)
  }
  
  has.padj <- which(unlist(lapply(set.list,function(x) length(which(colnames(x)=="padj"))))>0)
  pval.ht <- matrix(99,nrow=length(genes),ncol=length(has.padj))
  row.names(pval.ht) <- genes 
  colnames(pval.ht) <- names(has.padj)
  
  for(k in has.padj)
  {
    pval.ht[intersect(genes,row.names(set.list[[k]])),names(has.padj)[k]] <- as.numeric(unlist(lapply(strsplit(as.vector(set.list[[k]][intersect(genes,row.names(set.list[[k]])),]$padj),","),function(x) x[1])))
  }
  
  pval.ht.min <- apply(pval.ht,1,min)
  
  length(pval.ht.min)
  dim(supp)
  
  #get gene expr
  expr.mat <- matrix(NA,nrow=length(genes),ncol=length(expr.list))
  row.names(expr.mat) <- genes
  colnames(expr.mat) <- names(expr.list)
  for(k in 1:length(expr.list))
  {
    expr.vec <- expr.list[[k]]
    expr.mat[intersect(genes,names(expr.vec)),names(expr.list)[k]] <- expr.vec[intersect(genes,names(expr.vec))]
  }
  
  expr.mat <- round(expr.mat,2)
  
  gene.names <- genes
  names(gene.names) <- genes
  gene.names[intersect(genes,names(gene.ids))] <- gene.ids[intersect(genes,names(gene.ids))]
  #put together
  res.df <- data.frame(genes,gene.names,supp,pval.ht.min,expr.mat)
  res.df <- res.df[order(res.df$pval.ht.min),]
  return(res.df)
}


####################################################################################################
####################################################################################################
####################################################################################################
library(Hmisc)



getMatchedSet <- function(expr.vec,genes_targets,genes_nontargets,groups=10,keep.targets=F)
{
  mycut <- cut2(expr.vec,g=groups) 
  table(mycut)
  levels(mycut) <- c(1:length(levels(mycut)))
  expr.df <- data.frame("expr"=expr.vec,"bin"=mycut)
  row.names(expr.df) <- names(expr.vec)
  
  expr.df$is.target <- FALSE
  expr.df[genes_targets,]$is.target <- TRUE
  
  expr.df$is.nontarget <- FALSE
  expr.df[genes_nontargets,]$is.nontarget <- TRUE
  
  
  #head(expr.df)
  print(table(expr.df$bin[expr.df$is.target]))
  print(table(expr.df$bin[expr.df$is.nontarget]))
  
  pos.vec <- c()
  neg.vec <- c()
  expr.cat <- c()
  for(i in 1:length(levels(mycut)))
  {
    is.not.ta <- which(expr.df$bin==i & expr.df$is.nontarget==TRUE)
    is.ta <- which(expr.df$bin==i & expr.df$is.target==TRUE)
    
    if(!keep.targets)
    {
      to.samp <- min(c(length(is.not.ta),length(is.ta)))
      samp.ta <- sample(is.ta,to.samp,replace=FALSE)
      samp.not.ta <- sample(is.not.ta,to.samp,replace=FALSE)
    }else{
      to.samp <- length(is.ta)
      samp.ta <- is.ta
      samp.not.ta <- sample(is.not.ta,to.samp,replace=TRUE)
    }
    pos.vec <- c(pos.vec,samp.ta)
    neg.vec <- c(neg.vec,samp.not.ta)
    expr.cat <- c(expr.cat,rep(i,(to.samp)))
  }
  pos.genes <- row.names(expr.df)[pos.vec]
  neg.genes <- row.names(expr.df)[neg.vec]
  res <- data.frame("targets"=pos.genes,"non_targets"=neg.genes,"bin"=expr.cat)
  #head(res)
  return(res)
}

expr.list <- expr.list_ECT2
genes.roots <- genes.list$genes.roots
expr.vec <- expr.list[["roots.expr"]]
expr.vec <- expr.vec[genes.list$genes.roots.expressed]

matched.roots <- getMatchedSet(expr.vec = expr.vec,genes_targets = genes.roots,genes_nontargets = genes.list$genes.nontargets.roots,groups=10,keep.targets = T)
dim(matched.roots)
head(matched.roots)

genes.roots.intersect <- genes.list$genes.roots.intersect
matched.roots.intersect <- getMatchedSet(expr.vec = expr.vec,genes_targets=genes.roots.intersect,genes_nontargets = genes.list$genes.nontargets.roots,groups=10,keep.targets=T)
dim(matched.roots.intersect)

save.dir <- "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/matched_sets_paper/"

write.xlsx2(matched.roots, file=paste(save.dir,"MATCHED_SET_targets_non_targets_ROOTS_permissive.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

write.xlsx2(matched.roots.intersect, file=paste(save.dir,"MATCHED_SET_targets_non_targets_ROOTS_stringent.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)


expr.vec <- (expr.list[["roots.expr"]]+expr.list[["shoots.expr"]])/2
expr.vec <- expr.vec[genes.list$genes.expressed]
genes.union <- genes.list$genes.union

matched.union <- getMatchedSet(expr.vec = expr.vec,genes_targets=genes.union,genes_nontargets = genes.list$genes.nontargets,groups=10,keep.targets = T)
dim(matched.union)

genes.intersect <- genes.list$genes.intersect
matched.intersect <- getMatchedSet(expr.vec = expr.vec,genes_targets=genes.intersect,genes_nontargets = genes.list$genes.nontargets,groups=10,keep.targets=T)
dim(matched.intersect)

write.xlsx2(matched.union, file=paste(save.dir,"MATCHED_SET_targets_non_targets_ALL_permissive.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

write.xlsx2(matched.intersect, file=paste(save.dir,"MATCHED_SET_targets_non_targets_ALL_stringent.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)


matched.set.list <- list("union"=matched.union,"intersect"=matched.intersect,"roots"=matched.roots,"roots.intersect"=matched.roots.intersect)
#save.dir <- "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/"
save(matched.set.list,file=paste0(save.dir,"matched.set.list.Rdat"))







############# ECT3+ECT2 based for protoplast analysis

posGRF <- pos.list.ECT2

my.df.roots.2 <- positionsToGenes(posGR=pos.list.ECT2[["roots"]],include.ids = T,gene.ids=ids,use.meta = T)
my.df.roots.3 <- positionsToGenes(posGR=pos.list.ECT3[["roots"]],include.ids = T,gene.ids=ids,use.meta = T)
my.df.iclip <- positionsToGenes(posGR=iCLIP[[3]],include.ids = T,gene.ids=ids,use.meta = F)

#make my.df.expr.roots and my.df.expr.shoots
set.list <- list("HT.roots.2"=my.df.roots.2,"HT.roots.3"=my.df.roots.3,"iCLIP"=my.df.iclip)

#shoots.expr <- rowMeans(tpm.mat.collapsed[,grep("Sc",colnames(tpm.mat.collapsed))])
roots.expr <- rowMeans(tpm.mat.collapsed[,grep("Rc",colnames(tpm.mat.collapsed))])
expr.list <- list("roots.expr"=roots.expr)

genes.roots <- union(row.names(my.df.roots.2),row.names(my.df.roots.3))
#genes.single <- unique(c(row.names(my.df.roots),row.names(my.df.shoots)))
genes.union <- unique(c(genes.roots,row.names(my.df.iclip)))
genes.stringent <- unique(intersect(genes.roots,row.names(my.df.iclip)))
genes.super.stringent <- Reduce(intersect,list(row.names(my.df.roots.2),row.names(my.df.roots.3),row.names(my.df.iclip)))
#genes.roots.intersect <- unique(intersect(genes.roots,row.names(my.df.iclip)))
#genes.shoots.expressed <- names(which(apply(tpm.mat.collapsed[,grep("Sc",colnames(tpm.mat.collapsed))],1,function(x) length(x[x>1]))>1))
genes.roots.expressed <- names(which(apply(tpm.mat.collapsed[,grep("Rc",colnames(tpm.mat.collapsed))],1,function(x) length(x[x>1]))>1))
#genes.expressed <- unique(union(genes.shoots.expressed,genes.roots.expressed))
genes.nontargets.roots <- genes.expressed[!(genes.expressed %in% genes.union)]

genes.list.proto <- list("genes.roots"=genes.roots,
                   "genes.union"=genes.union,
                   "genes.stringent"=genes.stringent,
                   "genes.super.stringent"=genes.super.stringent,
                   "genes.roots.expressed"=genes.roots.expressed,
                   "genes.nontargets.roots"=genes.nontargets.roots)


genes.list.proto <- lapply(genes.list.proto,function(x) x[!(x %in% grep("_",x,value=T))])
unlist(lapply(genes.list.proto,length))

save(genes.list.proto,file = "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/genes.list_proto.Rdat")


##################################################
########### MATCHED SETS #########################
##################################################

expr.list <- expr.list_ECT2
genes.union <- genes.list.proto$genes.union
expr.vec <- expr.list[["roots.expr"]]
expr.vec <- expr.vec[genes.list.proto$genes.roots.expressed]

matched.roots <- getMatchedSet(expr.vec = expr.vec,genes_targets=genes.union,genes_nontargets = genes.list.proto$genes.nontargets.roots,groups=10,keep.targets = T)
dim(matched.roots)
head(matched.roots)

genes.roots.intersect <- genes.list.proto$genes.stringent
length(genes.roots.intersect)
matched.roots.intersect <- getMatchedSet(expr.vec = expr.vec,genes_targets=genes.roots.intersect,genes_nontargets = genes.list.proto$genes.nontargets.roots,groups=10,keep.targets=T)
dim(matched.roots.intersect)

save.dir <- "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/matched_sets_paper/"

write.xlsx2(matched.roots, file=paste(save.dir,"MATCHED_SET_targets_non_targets_PROTO_ROOTS_permissive.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

write.xlsx2(matched.roots.intersect, file=paste(save.dir,"MATCHED_SET_targets_non_targets_PROTO_ROOTS_stringent.xlsx",sep=""), sheetName = "genes",
            col.names = TRUE, row.names = FALSE, append = FALSE)

matched.set.list.proto <- list("union"=matched.roots,"intersect"=matched.roots.intersect)
save(matched.set.list.proto,file = "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/matched_sets_paper/matched.set.list.proto.Rdat")
lapply(matched.set.list.proto,nrow)
lapply(genes.list.proto,length)


################################################################################
###################Length based matched sets ###################################
################################################################################


library(openxlsx)
genes_info_lengths <- read.xlsx("/binf-isilon/alab/projects/ECT2_TC/transcript_lengths/full_list_arab_genes.xlsx", sheet=2, rowNames = FALSE)
head(genes_info_lengths)
tlen <- genes_info_lengths$transcript_length
names(tlen) <- genes_info_lengths$gene_id
tlen <- tlen[unique(genes_info_lengths$gene_id)]
expr.vec <- tlen
length(tlen)

library(Hmisc)
matched_length_roots_stringent <- getMatchedSet(expr.vec = tlen,genes_targets =  genes.list$genes.roots.intersect,genes_nontargets = genes.list$genes.nontargets.roots,groups=10,keep.targets = T)
head(matched_length_roots_stringent)

summary(genes.list$genes.roots %in% names(tlen))
matched_length_roots_large <- getMatchedSet(expr.vec = tlen,genes_targets =  genes.list$genes.roots[genes.list$genes.roots %in% names(tlen)],genes_nontargets = genes.list$genes.nontargets.roots,groups=10,keep.targets = T)
head(matched_length_roots_large)

summary(genes.list$genes.union %in% names(tlen))
matched_length_union <- getMatchedSet(expr.vec = tlen,genes_targets =  genes.list$genes.union[genes.list$genes.union %in% names(tlen)],genes_nontargets = genes.list$genes.nontargets,groups=10,keep.targets = T)
head(matched_length_roots_large)

summary(genes.list$genes.intersect %in% names(tlen))
matched_length_intersect <- getMatchedSet(expr.vec = tlen,genes_targets =  genes.list$genes.intersect[genes.list$genes.intersect %in% names(tlen)],genes_nontargets = genes.list$genes.nontargets,groups=10,keep.targets = T)
head(matched_length_intersect)

matched.set.list.length <- list("union"=matched_length_union,"intersect"=matched_length_intersect,"roots"=matched_length_roots_large,"roots.intersect"=matched_length_roots_stringent)

save(matched.set.list.length,file = "/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/matched_sets_paper/matched.set.list.length.Rdat")


##########
load("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/matched_sets_paper/matched.set.list.Rdat")
load("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/genes.list_ECT2.Rdat")

gene <- "AT1G16780"
grep(gene,matched.set.list$roots$non_targets,value=F)
grep(gene,matched.set.list$roots$targets,value=F)
grep(gene,genes.list$genes.nontargets.roots)

gene <- "AT5G03110"
grep(gene,matched.set.list$roots$non_targets,value=F)
grep(gene,matched.set.list$roots$targets,value=F)
grep(gene,genes.list$genes.nontargets.roots)

length(intersect(matched.set.list$roots$non_targets,genes.list$genes.nontargets.roots))
length(intersect(matched.set.list$roots$targets,genes.list$genes.nontargets.roots))
