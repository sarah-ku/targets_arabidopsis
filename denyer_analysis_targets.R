setwd("/binf-isilon/alab/projects/ECT2_TC/protoplast_markers/denyer_markers/")

#open denyer markers
denyer <- read.table("/binf-isilon/alab/projects/ECT2_TC/protoplast_markers/denyer_markers/denyer_markers_large.txt",sep="\t",header=T)
denyer <- denyer[denyer$logFC>=0.5,]

dnames <- read.table("/binf-isilon/alab/projects/ECT2_TC/protoplast_markers/denyer_markers/cluster_names")
dnames <- as.vector(dnames[,6])
names(dnames) <- paste(0:14)

denyer_sub <- read.table("/binf-isilon/alab/projects/ECT2_TC/protoplast_markers/denyer_markers/denyer_markers_large_subclusters.txt",sep="\t",header=T)

#open denyer single cell root matrix
scm <- read.csv2("/binf-isilon/alab/projects/ECT2_TC/protoplast_markers/denyer_markers/GSE123818_Root_single_cell_wt_datamatrix.csv",sep=",",header = T,row.names = 1)
scm <- as.matrix(scm)

#generate cell assignments
denyer$target <- TRUE
denyer$cluster.name <- dnames[paste(as.vector(denyer$cluster))]

######################
m6agenes <- unique(as.vector(nanoGR$gene))
denyer$cluster.name <- dnames[paste(as.vector(denyer$cluster))]

target_genes <- genes.list.proto$genes.roots
denyer$target <- as.vector(denyer$id) %in% target_genes

denyer$target.strict <- as.vector(denyer$id) %in% genes.list.proto$genes.stringent
denyer$m6a <- as.vector(denyer$id) %in% m6agenes
denyer.m6a <- denyer[denyer$m6a,]
mean.target <- mean(as.vector(unique(denyer$id)) %in% target_genes)
mean.target.strict <- mean(as.vector(unique(denyer$id)) %in% genes.list.proto$genes.stringent)
mean.target.m6a <- mean(as.vector(unique(denyer.m6a$id)) %in% genes.list.proto$genes.roots)
mean.target.strict.m6a <- mean(as.vector(unique(denyer.m6a$id)) %in% genes.list.proto$genes.stringent)
denyer$m6a_only <- denyer$m6a & (!denyer$target)

setwd("/binf-isilon/alab/projects/ECT2_TC/protoplast_markers/denyer_markers/")
pdf("denyer_single_cell_large_marker_proportions_main_targs_m6a_ECT2_AND_ECT3.pdf",width=11,height=4)
par(mfrow=c(1,3))
mycols <- c(brewer.pal(9,"Set1"),brewer.pal(8,"Set2"))[1:length(table(ptypes))]

tab <- rbind(tapply(denyer$target,denyer$cluster,mean),1-tapply(denyer$target,denyer$cluster,mean))
colnames(tab) <- paste(dnames[paste(colnames(tab))],paste(colnames(tab)))
barplot(tab,ylim=c(0,1),col=c("dark grey","light grey"),las=2,ylab="Proportion of markers which are targets (roots all)",xlab="",main="prop. targets")
mbarp <- barplot(tab[1,],col=mycols,add=T,axes=F,names.arg = "")
abline(h=mean.target,lty=2)
text(x=mbarp, y = .9, label = round(tab[1,],2), pos = 3,las=2, cex = 0.8, col = "red")


tab <- rbind(tapply(denyer$target.strict,denyer$cluster,mean),1-tapply(denyer$target.strict,denyer$cluster,mean))
colnames(tab) <- paste(dnames[paste(colnames(tab))],paste(colnames(tab)))
barplot(tab,ylim=c(0,1),col=c("dark grey","light grey"),las=2,ylab="Proportion of markers which are targets (roots strict)",xlab="",main="prop. strict targets")
mbarp <- barplot(tab[1,],col=mycols,add=T,axes=F,names.arg = "")
abline(h=mean.target.strict,lty=2)
text(x=mbarp, y = .9, label = round(tab[1,],2), pos = 3,las=2, cex = 0.8, col = "red")

tab <- rbind(tapply(denyer$m6a,denyer$cluster,mean),1-tapply(denyer$m6a,denyer$cluster,mean))
colnames(tab) <- paste(dnames[paste(colnames(tab))],paste(colnames(tab)))
barplot(tab,col=c("dark grey","light grey"),ylim=c(0,1),las=2,ylab="Proportion of markers which are m6A (roots all)",xlab="",main="prop. m6A")
mbarp <-barplot(tab[1,],col=mycols,add=T,axes=F,names.arg = "")
abline(h=mean.target.m6a,lty=2)
text(x=mbarp, y = .9, label = round(tab[1,],2), pos = 3,las=2, cex = 0.8, col = "red")

dev.off()