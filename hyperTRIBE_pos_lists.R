setwd("/binf-isilon/alab/projects/ECT2_TC/hyperTRIBE/")

load(file="posGRL_list.Rdat")
#posGRL[[1]]
posGR <- list()
for(i in 1:length(posGRL))
{
  posGR[[i]] <- posGRL[[i]][posGRL[[i]]$padj<0.01 & posGRL[[i]]$fold_change>1 & posGRL[[i]]$prop<1]
  posGR[[i]] <- posGR[[i]][paste(posGR[[i]]$ref,posGR[[i]]$targ) %in% c("A G", "T C")]
  posGR[[i]] <- posGR[[i]][posGR[[i]]$tags_treat>10]
}
names(posGR) <- names(posGRL)
pos.list.ECT2 <- posGR
save(pos.list.ECT2,file="pos_list_ECT2_final.Rdat")

setwd("/binf-isilon/alab/projects/ECT2_TC/ECT3_hyperTRIBE/")

load(file="posGRL_list.Rdat")
#posGRL[[1]]
posGR <- list()
for(i in 1:length(posGRL))
{
  posGR[[i]] <- posGRL[[i]][posGRL[[i]]$padj<0.01 & posGRL[[i]]$fold_change>1 & posGRL[[i]]$prop<1]
  posGR[[i]] <- posGR[[i]][paste(posGR[[i]]$ref,posGR[[i]]$targ) %in% c("A G", "T C")]
  posGR[[i]] <- posGR[[i]][posGR[[i]]$tags_treat>10]
}
names(posGR) <- names(posGRL)
pos.list.ECT3 <- posGR
save(pos.list.ECT3,file="pos_list_ECT3_final.Rdat")
