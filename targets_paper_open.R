library(GenomicRanges)
library(reshape2)
library(rtracklayer)
library(systemPipeR)
library(eulerr)
library(RNAeditR)

library(GeneOverlap)
library(gridExtra)
library(ggplot2)
library(Biostrings)
library(BSgenome.Athaliana.TAIR.TAIR9)
require(ggseqlogo)
library(Hmisc)

library(RColorBrewer)
#library(viridis)

load("/binf-isilon/alab/projects/ECT2_TC/binding sites unfiltered/iCLIPBS_with_meta.Rdat")


setwd("/binf-isilon/alab/projects/ECT2_TC/")

load(file="./hyperTRIBE/pos_list_ECT2_final.Rdat")
load(file="./hyperTRIBE/pos_list_ECT3_final.Rdat")
#summary(unique(pos.list.ECT3$roots$gene) %in% unique(pos.list.ECT2$roots$gene))

load(file="./hyperTRIBE/iCLIP/ECT2_iCLIP_181109/BED/called_peaks_broad/iCLIP_with_meta.Rdat")
#summary(intersect(iCLIP[[3]]$gene,pos.list.ECT3$roots$gene) %in% intersect(iCLIP[[3]]$gene,pos.list.ECT2$roots$gene))
load(file="./iCLIP/newSites/iclipBS.Rdat")

load(file="./parker19/parker19_nanoGR.Rdat")
load(file="./parker19/parker19_micGR.Rdat")
#names(nanoGR) <- paste0(as.vector(seqnames(nanoGR)),end(nanoGR))
names(micGR) <- paste0(as.vector(seqnames(micGR)),end(micGR))
start(nanoGR) <- end(nanoGR)
names(nanoGR) <- paste0(as.vector(seqnames(nanoGR)),"_",start(nanoGR))




load(file="/binf-isilon/alab/projects/ECT2_TC/FA_CLIP/weiGRC_singlebp.Rdat")
#weiGRC

load(file="./hyperTRIBE/m6A/m6AGR_centered.Rdat")


load("/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/tpm.mat.Rdat")
tpm.mat3 <- tpm.mat

load("/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/tpm.mat.collapsed.Rdat")
tpm.mat.collapsed3 <- tpm.mat.collapsed

load("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/quantification/transcripts_matrix_araport11.Rdat")
#head(tpm.mat)
#colnames(tpm.mat)[colnames(tpm.mat)=="E2T_Se5"] <- "E2T_St5"

load(file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/quantification/tpm.mat.collapsed.Rdat")
#colnames(tpm.mat.collapsed)[colnames(tpm.mat.collapsed)=="E2T_Se5"] <- "E2T_St5"
#head(tpm.mat.collapsed)

gtf <- rtracklayer::import("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/iCLIP/Araport11_GFF3_genes_transposons.201606.gtf")
seqlevels(gtf) <- c("1","2","3","4","5","Pt","Mt")
gtfGR <- gtf


load(paste0("./hyperTRIBE/gene_ids.Rdat"))
head(ids)

gene.desc <- read.csv("./annotations/gene_descriptions_arabidopsis.tsv",sep="\t",header=F)
row.names(gene.desc) <- gene.desc$V1

gtfE <- gtf
colset <- brewer.pal(8,"Set1")
mycols <- brewer.pal(8,"Set1")

load(file="./hyperTRIBE/targets/matched.set.list.Rdat")
load(file="./hyperTRIBE/targets/target.set.list_new.Rdat")
load(file="./hyperTRIBE/targets/genes.list_ECT2.Rdat")
load(file="./hyperTRIBE/targets/genes.list_proto.Rdat")
load(file="./hyperTRIBE/targets/expr.list_ECT2.Rdat")

#names(expr.list_ECT2)
#expr.list_ECT3 <- list("shoots.expr"=rowMeans(tpm.mat.collapsed3[,grep("Sc",colnames(tpm.mat.collapsed3))]),"roots.expr"=rowMeans(tpm.mat.collapsed3[,grep("Rc",colnames(tpm.mat.collapsed3))]))
#save(expr.list_ECT3,file="/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/targets/expr.list_ECT3.Rdat")
#tpm.mat.collapsed3

load(file="./hyperTRIBE/targets/expr.list_ECT3.Rdat")
load(file="/binf-isilon/alab/projects/ECT2_TC/anderson_et_al_andGR.Rdat")

htl <- unique(unlist(GRangesList(pos.list.ECT2[1:2])))


load(paste0("/binf-isilon/PBgrp/qbp693/ect2_protoplast/output_rnaseq_smartseq_analyses/pf2_de_output_smartseq.Rdata"))

load("/binf-isilon/alab/projects/ECT2_TC/binding sites unfiltered/iCLIPBS_with_meta.Rdat")



makeEulerDiagram <- function(sets_list, color_palette, plot_name, plot_title){
  
  plot_id <- as.symbol("plot_name")
  venset <- overLapper(sets_list, type="vennsets")
  plot_id <- plot(
    euler(intersectmatrix(venset), regionError = T),
    fills=list(fill=color_palette, alpha=0.8),
    edges=list(lwd=0.3),
    labels = list(fontfamily = "sans",
                  main.fontfamily = "sans",
                  cat.fontfamily = "sans"),
    legend = TRUE,
    main=plot_title,
    quantities=TRUE,
    percentages = TRUE)
  
  plot_id
}
