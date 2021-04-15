#ect3 script

#FROM SORTED BAM FILES START HERE: start with samtools mpileup
#mpileup + compile counts for roots and shoots separately:
#samtools mpileup --max-depth 50000 -Q 30 --skip-indels -f /binf-isilon/alab/projects/ECT2_TC/annotations/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa /binf-isilon/PBgrp/qbp693/hyperTRIBE_ect3/output/*R*sort.bam | perl /home/sarah/R_code/ECT2/hyperTRIBE_mpileup2bases.pl> baseCounts_roots_hyperTRIBE_ECT3.txt &
#samtools mpileup --max-depth 50000 -Q 30 --skip-indels -f /binf-isilon/alab/projects/ECT2_TC/annotations/Arabidopsis_thaliana/Ensembl/TAIR10/Sequence/WholeGenomeFasta/genome.fa /binf-isilon/PBgrp/qbp693/hyperTRIBE_ect3/output/*S*sort.bam | perl /home/sarah/R_code/ECT2/hyperTRIBE_mpileup2bases.pl> baseCounts_shoots_hyperTRIBE_ECT3.txt &

################################ needed later



library(tximport)

 my.path <- paste0("/binf-isilon/PBgrp/qbp693/hyperTRIBE_ect3/quantifications_araport/quants/")
 filenames <- list.files(path = my.path,full.names = T)
# #filenames <- filenames.list[[i]]
 filenames <- paste0(filenames,"/quant.sf")
# 
 transcripts <- as.vector(read.table(filenames[1])[,1])
 genes <- gsub("\\.[0-9]+","",transcripts)
 gene.df <- data.frame("transcript"=transcripts,"gene"=genes)
# 
 txi <- tximport(files = filenames, tx2gene=gene.df,type = "salmon")
# tpm.mat <- (txi$abundance)
# colnames(tpm.mat) <- gsub(".*/(.*)_quant/quant.sf","\\1",filenames)
# save(tpm.mat.collapsed,file="/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/tpm.mat.collapsed.Rdat")
# 
# gene.df2 <- data.frame("transcript"=transcripts,"gene"=transcripts)
# txi <- tximport(files = filenames, tx2gene=gene.df2,type = "salmon")
# tpm.mat <- (txi$abundance)
# colnames(tpm.mat) <- gsub(".*/(.*)_quant/quant.sf","\\1",filenames)
# save(tpm.mat,file="/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/tpm.mat.Rdat")

#open quantifications
#load("/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/tpm.mat.collapsed.Rdat")
load("/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/tpm.mat.Rdat")

head(tpm.mat)
barplot(tpm.mat[1,],las=2)
#open annotations
gtf <- rtracklayer::import("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/iCLIP/Araport11_GFF3_genes_transposons.201606.gtf")
seqlevels(gtf) <- c("1","2","3","4","5","Pt","Mt")

load("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/gene_ids.Rdat")

################################### for modelling

dat_shoots <- read.table("/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/baseCounts_shoots_hyperTRIBE_ECT3.txt",header=F)
dat_roots <- read.table("/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/baseCounts_roots_hyperTRIBE_ECT3.txt",header=F)
dim(dat_shoots)
dim(dat_roots)

locsGR_shoots <- GRanges(Rle(dat_shoots$V1),IRanges(dat_shoots$V2,width=1),ref=dat_shoots$V3,names=paste(dat_shoots$V1,dat_shoots$V2,sep="_"))
locsGR_roots <- GRanges(Rle(dat_roots$V1),IRanges(dat_roots$V2,width=1),ref=dat_roots$V3,names=paste(dat_roots$V1,dat_roots$V2,sep="_"))

filenames <- list.files(path="/binf-isilon/PBgrp/qbp693/hyperTRIBE_ect3/output/",pattern = "sort.bam$")

samp.names_roots <- grep("E3T_R",filenames,value=T)
samp.names_roots <- gsub(".sort.bam","",samp.names_roots)

samp.names_shoots <- grep("E3T_S",filenames,value=T)
samp.names_shoots <- gsub(".sort.bam","",samp.names_shoots)

data_list_roots <- extractCountData(dat_roots,samp.names_roots,strand=F)
data_list_shoots <- extractCountData(dat_shoots,samp.names_shoots,strand=F)

#check data is populated for all samples
lapply(data_list_roots,nrow)
lapply(data_list_shoots,nrow)



#for the roots we remove Re8 since it's got very low coverage (indicating a problem with the sample)
#data_list_roots <- data_list_roots[-which(names(data_list_roots)=="E2T_Re8")]

#now produce one design vector per experiment (total of 4 for the roots and the shoots, single and triple mutants combinations)
design_vector_roots_single <- c(E3T_Rc11 = "control", E3T_Rc12 = "control", E3T_Rc13 = "control", 
                                E3T_Rc14 = "control", E3T_Rc15 = "control", E3T_Re6 = "treat", 
                                E3T_Re7 = "treat",E3T_Re8 = "treat", E3T_Re9 = "treat", E3T_Re10 = "treat")
table(design_vector_roots_single)

design_vector_roots_triple <- c(E3T_Rc11 = "control", E3T_Rc12 = "control", E3T_Rc13 = "control", 
                                E3T_Rc14 = "control", E3T_Rc15 = "control", E3T_Rt1 = "treat", 
                                E3T_Rt2 = "treat", E3T_Rt3 = "treat", E3T_Rt4 = "treat", E3T_Rt5 = "treat")
table(design_vector_roots_triple)


design_vector_shoots_single <- c(E3T_Sc11 = "control", E3T_Sc12 = "control", E3T_Sc13 = "control", 
                                 E3T_Sc14 = "control", E3T_Sc15 = "control", E3T_Se6 = "treat", 
                                 E3T_Se7 = "treat", E3T_Se8 = "treat", E3T_Se9 = "treat",E3T_Se10 = "treat")
table(design_vector_shoots_single)


design_vector_shoots_triple <- c(E3T_Sc11 = "control", E3T_Sc12 = "control", E3T_Sc13 = "control", 
                                 E3T_Sc14 = "control", E3T_Sc15 = "control", E3T_St1 = "treat", E3T_St2 = "treat", E3T_St3 = "treat",
                                 E3T_St4 = "treat", E3T_St5 = "treat")
table(design_vector_shoots_triple)

design_vector_roots_single_triple <- c(E3T_Re6 = "control", E3T_Re7 = "control", E3T_Re9 = "control", E3T_Re10 = "control", 
                                      E3T_Rt1 = "treat", E3T_Rt2 = "treat", E3T_Rt3 = "treat", E3T_Rt4 = "treat", E3T_Rt5 = "treat")
table(design_vector_roots_single_triple)

design_vector_shoots_single_triple <- c( E3T_Se6 = "control", E3T_Se7 = "control", E3T_Se8 = "control", E3T_Se9 = "control",E3T_Se10 = "control", 
                                         E3T_St1 = "treat", E3T_St2 = "treat", E3T_St3 = "treat", E3T_St4 = "treat", E3T_St5 = "treat")

table(design_vector_shoots_single_triple)

my_edits <- rbind(c("A","G"),
                  c("G","A"),
                  c("T","C"),
                  c("C","T"),
                  c("A","T"),
                  c("T","A"),
                  c("G","T"),
                  c("T","G"),
                  c("G","C"),
                  c("C","G"),
                  c("A","C"),
                  c("C","A"))
