library(GenomicRanges)
library(reshape2)
library(rtracklayer)
library(systemPipeR)
library(eulerr)
library(RNAeditR)

setwd("/binf-isilon/alab/projects/ECT2_TC/")

load(file="./hyperTRIBE/pos_list_ECT2_final.Rdat")
load(file="./hyperTRIBE/pos_list_ECT3_final.Rdat")
#summary(unique(pos.list.ECT3$roots$gene) %in% unique(pos.list.ECT2$roots$gene))

load(file="./hyperTRIBE/iCLIP/ECT2_iCLIP_181109/BED/called_peaks_broad/iCLIP_with_meta.Rdat")
#summary(intersect(iCLIP[[3]]$gene,pos.list.ECT3$roots$gene) %in% intersect(iCLIP[[3]]$gene,pos.list.ECT2$roots$gene))
load(file="./iCLIP/newSites/iclipBS.Rdat")

load(file="./parker19/parker19_nanoGR.Rdat")
load(file="./parker19/parker19_micGR.Rdat")
names(nanoGR) <- paste0(as.vector(seqnames(nanoGR)),end(nanoGR))
names(micGR) <- paste0(as.vector(seqnames(micGR)),end(micGR))


load(file="./hyperTRIBE/m6A/m6AGR_centered.Rdat")
load("./hyperTRIBE/quantification/transcripts_matrix_araport11.Rdat")

load(file="./hyperTRIBE/quantification/tpm.mat.collapsed.Rdat")
colnames(tpm.mat.collapsed)[colnames(tpm.mat.collapsed)=="E2T_Se5"] <- "E2T_St5"

gtf <- rtracklayer::import("./hyperTRIBE/iCLIP/Araport11_GFF3_genes_transposons.201606.gtf")
seqlevels(gtf) <- c("1","2","3","4","5","Pt","Mt")
gtfGR <- gtf


load(paste0("./hyperTRIBE/gene_ids.Rdat"))
head(ids)

gene.desc <- read.csv("./annotations/gene_descriptions_arabidopsis.tsv",sep="\t",header=F)
row.names(gene.desc) <- gene.desc$V1


load(file="./hyperTRIBE/targets/matched.set.list.Rdat")
load(file="./hyperTRIBE/targets/target.set.list_new.Rdat")
load(file="./hyperTRIBE/targets/genes.list_ECT2.Rdat")
load(file="./hyperTRIBE/targets/genes.list_proto.Rdat")
load(file="./hyperTRIBE/targets/expr.list_ECT2.Rdat")




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
