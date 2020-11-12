setwd("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE")
library(RNAeditR)

################################### for modelling
dat_shoots <- read.table("./pipeline/output/baseCounts_shoots_hyperTRIBE.txt",header=F)
dat_roots <- read.table("./pipeline/output/baseCounts_roots_hyperTRIBE.txt",header=F)
dim(dat_shoots)
dim(dat_roots)

locsGR_shoots <- GRanges(Rle(dat_shoots$V1),IRanges(dat_shoots$V2,width=1),ref=dat_shoots$V3,names=paste(dat_shoots$V1,dat_shoots$V2,sep="_"))
locsGR_roots <- GRanges(Rle(dat_roots$V1),IRanges(dat_roots$V2,width=1),ref=dat_roots$V3,names=paste(dat_roots$V1,dat_roots$V2,sep="_"))


samp.names_roots <- c("E2T_Rc11","E2T_Rc12","E2T_Rc13","E2T_Rc14","E2T_Rc15",
                      "E2T_Re6","E2T_Re7","E2T_Re8","E2T_Re9","E2T_Re10",
                      "E2T_Rt1", "E2T_Rt2", "E2T_Rt3" ,"E2T_Rt4" ,"E2T_Rt5")

samp.names_shoots <- c("E2T_Sc11","E2T_Sc12","E2T_Sc13","E2T_Sc14","E2T_Sc15",
                       "E2T_Se6", "E2T_Se7", "E2T_Se8" , "E2T_Se9",  "E2T_Se10",
                       "E2T_St1", "E2T_St2", "E2T_St3" ,"E2T_St4" ,"E2T_Se5")



data_list_roots <- extractCountData(dat_roots,samp.names_roots,strand=F)
data_list_shoots <- extractCountData(dat_shoots,samp.names_shoots,strand=F)

#check data is populated for all samples
lapply(data_list_roots,nrow)
lapply(data_list_shoots,nrow)


#for the roots we remove Re8 since it's got very low coverage (indicating a problem with the sample)
data_list_roots <- data_list_roots[-which(names(data_list_roots)=="E2T_Re8")]

#now produce one design vector per experiment (total of 4 for the roots and the shoots, single and triple mutants combinations)
design_vector_roots_single <- c(E2T_Rc11 = "control", E2T_Rc12 = "control", E2T_Rc13 = "control", 
                                E2T_Rc14 = "control", E2T_Rc15 = "control", E2T_Re6 = "treat", 
                                E2T_Re7 = "treat", E2T_Re9 = "treat", E2T_Re10 = "treat")
table(design_vector_roots_single)

design_vector_roots_triple <- c(E2T_Rc11 = "control", E2T_Rc12 = "control", E2T_Rc13 = "control", 
                                E2T_Rc14 = "control", E2T_Rc15 = "control", E2T_Rt1 = "treat", 
                                E2T_Rt2 = "treat", E2T_Rt3 = "treat", E2T_Rt4 = "treat", E2T_Rt5 = "treat")
table(design_vector_roots_triple)


design_vector_shoots_single <- c(E2T_Sc11 = "control", E2T_Sc12 = "control", E2T_Sc13 = "control", 
                                 E2T_Sc14 = "control", E2T_Sc15 = "control", E2T_Se6 = "treat", 
                                 E2T_Se7 = "treat", E2T_Se8 = "treat", E2T_Se9 = "treat",E2T_Se10 = "treat")
table(design_vector_shoots_single)


design_vector_shoots_triple <- c(E2T_Sc11 = "control", E2T_Sc12 = "control", E2T_Sc13 = "control", 
                                 E2T_Sc14 = "control", E2T_Sc15 = "control", E2T_Se5 = "treat", E2T_St1 = "treat", E2T_St2 = "treat",
                                 E2T_St4 = "treat", E2T_St3 = "treat")
table(design_vector_shoots_triple)


design_vector_roots_single_triple <- c(E2T_Re6 = "control", E2T_Re7 = "control", E2T_Re9 = "control", E2T_Re10 = "control",
                                       E2T_Rt1 = "treat", E2T_Rt2 = "treat", E2T_Rt3 = "treat", E2T_Rt4 = "treat", E2T_Rt5 = "treat")
table(design_vector_roots_single_triple)

design_vector_shoots_single_triple <- c( E2T_Se6 = "control", E2T_Se7 = "control", E2T_Se8 = "control", E2T_Se9 = "control",E2T_Se10 = "control", 
                                         E2T_Se5 = "treat",E2T_St1 = "treat", E2T_St2 = "treat", E2T_St4 = "treat", E2T_St3 = "treat")

table(design_vector_shoots_single_triple)
