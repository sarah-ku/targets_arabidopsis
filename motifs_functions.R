
returnDens <- function(tGR,gtfGR){
  props_list_plus <- list()
  props_list_minus <- list()
  props_list <- list()
  
  names(tGR) <- paste0(seqnames(tGR),":",start(tGR),",",strand(tGR))
  
  for(type in c("5UTR","CDS","3UTR"))
  {
    cdsOL <- findOverlaps(tGR,gtfGR[gtfGR$type==type])
    gtf_hits <- gtfGR[gtfGR$type==type][subjectHits(cdsOL)]
    tgr_hits <- tGR[queryHits(cdsOL)]
    
    ss <- start(gtf_hits)
    sq <- start(tgr_hits)
    es <- end(gtf_hits)
    
    props_plus <- (sq-ss+1)/width(gtf_hits)
    names(props_plus) <- names(tgr_hits)
    props_plus <- props_plus[which(strand(tgr_hits)=="+")]
    
    props_minus <- (es - sq +1)/width(gtf_hits)
    names(props_minus) <- names(tgr_hits)
    props_minus <- props_minus[which(strand(tgr_hits)=="-")]
    
    props_plus <- tapply(props_plus,names(props_plus),function(x) mean(x))
    props_minus <- tapply(props_minus,names(props_minus),function(x) mean(x))
    
    props_list_plus[[type]] <- props_plus
    props_list_minus[[type]] <- props_minus
    props_list[[type]] <- c(props_plus,props_minus)
  }
  lst <- props_list_plus 
  pvec_p <- c(lst[[1]],lst[[2]]+1,lst[[3]]+2)
  lst <- props_list_minus
  pvec_m <- c(lst[[1]],lst[[2]]+1,lst[[3]]+2)
  tort <- c(pvec_p,pvec_m)
  #pvec_p <- c(lst[[4]])
  return(tort)
}

plotNE <- function(bs=0.005,dr_rand,dr,dr_ind,inc_groups=T,inc_int=T,to.add=F,col.no=1,no.samps=10,to.smooth=T,mytitle="")
{
  rtabcut <- as.vector(table(cut(dr_rand,breaks=seq(0,3,by=bs))))
  rtab <- as.vector(table(cut(dr,breaks=seq(0,3,by=bs))))/rtabcut
  nl <- length(dr)/length(dr_rand)
  rtab <- (rtab/nl)
  #plot(rtab,type="l")
  
  if(inc_int)
  {
    rtab_rand <- list()
    for(i in 1:no.samps)
    {
      dr_r <- dr[sample(1:length(dr),length(dr),replace=T)]
      rtab_rand[[i]] <- as.vector(table(cut(dr_r,breaks=seq(0,3,by=bs))))/rtabcut
      rtab_rand[[i]] <- rtab_rand[[i]]/nl
      #nl <- length(dr)/length(dr_rand)
    }
    
    rtab_mat <- do.call(rbind,rtab_rand)
    rtab_mean <- colMeans(rtab_mat)
    rtab_upper <- apply(rtab_mat,2,function(x) quantile(x,.95))
    rtab_lower <- apply(rtab_mat,2,function(x) quantile(x,.05))
    
    mymin <- min(rtab_lower)
    mymax <- max(rtab_upper)
    
  }else{
    mymin <- min(rtab)
    mymax <- max(rtab)
  }
  
  xls <- seq(0,3,by=bs)[-1]
  
  #summary(rtab_lower)
  
  if(to.add)
  {
    if(to.smooth)
    {
      points(y=smooth.spline(rtab)$y,x=xls,col=brewer.pal("Set1",n = 9)[col.no],type="l",pch=16,xaxt="n",
             ylab="Normalised enrichment",xlab="Genomic position",
             ylim=c(mymin,mymax),lty=1,lwd=1)
    }else{
      points(y=rtab,x=xls,col=brewer.pal("Set1",n = 9)[col.no],type="l",pch=16,xaxt="n",
             ylab="Normalised enrichment",xlab="Genomic position",
             ylim=c(mymin,mymax),lty=1,lwd=1)
    }
  }else{
    if(to.smooth)
    {
      plot(y=smooth.spline(rtab)$y,x=xls,col=brewer.pal("Set1",n = 9)[col.no],pch=16,xaxt="n",type="l",
             ylab="Normalised enrichment",xlab="Genomic position",
             ylim=c(mymin,mymax),lty=1,lwd=1)
    }else{
      plot(y=rtab,x=xls,col=brewer.pal("Set1",n = 9)[col.no],type="l",pch=16,xaxt="n",
           ylab="Normalised enrichment",xlab="Genomic position",
           ylim=c(mymin,mymax),lty=1,lwd=1,main=mytitle)
    }
  }
  
  
  tp <- rtab
  x <- c(xls,rev(xls))
  
  if(inc_int)
  {
    if(to.smooth)
    {
      y <- c(smooth.spline(rtab_lower)$y, rev(smooth.spline(rtab_upper)$y))
      
      polygon(x,y, col = brewer.pal("Set1",n = 9)[col.no], border = NA)
      points(y=smooth.spline(rtab_mean)$y,x=xls,col="black",lty=1,type="l")
      #points(y=smooth.spline(rtab_lower)$y,x=xls,type="l",col="grey",lty=2)
      #points(y=smooth.spline(rtab_upper)$y,x=xls,type="l",col="grey",lty=2)
    }else{
      y <- c(rtab_lower, rev(rtab_upper))
      
      polygon(x,y, col = brewer.pal("Set1",n = 9)[col.no], border = NA)
      points(y=rtab,x=xls,col="black",lty=1,type="l")
      #points(y=rtab_lower,x=xls,type="l",col="grey",lty=2)
      #points(y=rtab_upper,x=xls,type="l",col="grey",lty=2)
    }
  }
  #else{
    #points(y=rtab,x=xls,col=brewer.pal("Set1",n = 9)[col.no],lty=1,lwd=2,type="l")
  #}
  
  
  
  abline(h=1,lty=2)
  axis(1,c(0,1,2,3),c(0,1,2,3))
  abline(v=c(0,1,2,3),lty=2)
}

plotDens <- function(bs=0.005,dr_rand,dr,dr_ind,inc_groups=T,inc_int=T,to.add=F,col.no=1,no.samps=10)
{
  rtab <- as.vector(table(cut(dr,breaks=seq(0,3,by=bs*3))))/as.vector(table(cut(dr_rand,breaks=seq(0,3,by=bs*3))))

  
  rtab_rand <- list()
  for(i in 1:no.samps)
  {
    dr_r <- dr_rand[sample(1:length(dr_rand),length(dr_rand),replace=T)]
    rtab_rand[[i]] <- as.vector(table(cut(dr,breaks=seq(0,3,by=bs*3))))/as.vector(table(cut(dr_r,breaks=seq(0,3,by=bs*3))))
  }
  
  rtab_mat <- do.call(rbind,rtab_rand)
  rtab_upper <- apply(rtab_mat,2,function(x) quantile(x,.9))
  rtab_lower <- apply(rtab_mat,2,function(x) quantile(x,.1))
  
  rtab_upper_SMOOTH <- smooth.spline(rtab_upper)$y
  rtab_lower_SMOOTH <- smooth.spline(rtab_lower)$y
  

  rtabv <- list()
  for(i in 1:length(dr_ind))
  {
    rtabv[[i]] <- as.vector(table(cut(dr_ind[[i]],breaks=seq(0,3,by=bs))))/as.vector(table(cut(dr_rand,breaks=seq(0,3,by=bs))))
    nl <- length(dr_ind[[i]])/length(dr_rand)
    rtabv[[i]] <- rtabv[[i]]/nl
  }
  #}
  
  
  #rtab
  xls <- seq(0,3,by=bs*3)[-1]
  
  par(mfrow=c(1,1))
  tp <- rtab
  
  x <- c(xls,rev(xls))
  #y <- c(smooth.spline(rtab_lower)$y, rev(smooth.spline(rtab_upper)$y))
  y <- c(rtab_lower, rev(rtab_upper))
  
  
  rtabvi <- unlist(rtabv)[!is.infinite(unlist(rtabv))]
  
  if(inc_groups)
  {
    mymin <- min(c(rtab_lower,rtabvi))
    mymax <- max(c(rtab_upper,rtabvi))
  }else{
    mymin <- min(c(rtab_lower,rtabvi))
    mymax <- max(c(rtab_upper,rtabvi))+0.2
  }
  
  
  #plot(y=tp,x=xls,type="l",xaxt="n",lwd=1,col="red",ylab="Fold enrichment over CDS",xlab="Genomic position",ylim=c(0,15))
  #plot(y=smooth.spline(tp)$y,x=xls,type="l",xaxt="n",ylab="Fold enrichment over CDS",xlab="Genomic position",ylim=c(0.6,1.8),col="black",lty=2,lwd=2)
  if(to.add)
  {
    points(y=rtab,x=xls,col="grey",type="n",pch=16,xaxt="n",ylab="Normalised enrichment",xlab="Genomic position",ylim=c(min(c(rtab_lower,unlist(rtabv))),max(c(rtab_upper,unlist(rtabv)))),lty=2,lwd=1)
    
  }else{
    
    plot(y=rtab,x=xls,col="grey",type="n",pch=16,xaxt="n",ylab="Normalised enrichment",xlab="Genomic position",ylim=c(mymin,mymax),lty=2,lwd=1)
  }
  
  if(inc_int)
  {
    #polygon(x,y, col = brewer.pal("Set1",n = 6)[col.no], border = NA)
    #points(y=smooth.spline(rtab)$y,x=xls,col="black",lty=1,type="l")
    #points(y=rtab_lower_SMOOTH,x=xls,type="l",col="black",lty=2)
    #points(y=rtab_upper_SMOOTH,x=xls,type="l",col="black",lty=2)
    
    polygon(x,y, col = brewer.pal("Set1",n = 6)[col.no], border = NA)
    points(y=rtab,x=xls,col="black",lty=1,type="l")
    points(y=rtab_lower,x=xls,type="l",col="black",lty=2)
    points(y=rtab_upper,x=xls,type="l",col="black",lty=2)
    
  }
  
  #points(y=)
  abline(h=1,lty=2)
  
  axis(1,c(0,1,2,3),c(0,1,2,3))
  abline(v=c(0,1,2,3),lty=2)
  
  
  if(inc_groups)
  {
    for(i in 1:length(dr_ind))
    {
      tp <- rtabv[[i]]
      points(tp,x=xls,type="l",col=brewer.pal("Set1",n = 6)[i],lwd=1,lty=1)
    }
    legend("topleft",names(dr_ind),col=brewer.pal("Set1",n = 6)[1:length(dr_ind)],lwd=1,lty=1)
  }
  
}




getMatchedSets <- function(tGR,rGR=randGR,mdb=0.1)
{
  dr_list <- list()
  dr_list_break <- list()
  for(k in 1:length(tGR))
  {
    dr_list[[k]] <- returnDens(tGR[[k]],gtfGR)
    dr_list_break[[k]] <- as.numeric(as.vector(seq(0,3,by=mdb)[(as.numeric(cut(dr_list[[k]],seq(0,3,by=mdb))))]))
    names(dr_list_break[[k]]) <- names(dr_list[[k]])
    dr_list_break[[k]] <- dr_list_break[[k]][!(names(dr_list_break[[k]])=="")]
  }
  
  random_locations <- returnDens(rGR,gtfGR)
  #hist(random_locations)
  random_locations <- tapply(random_locations,names(random_locations),function(x) sample(rep(x,2),1))
  #hist(random_locations)
  
  seqvec <- seq(0,3,by=mdb)
  random_round <- as.numeric(as.vector(seqvec[(as.numeric(cut(random_locations,seqvec)))]))
  names(random_round) <- names(random_locations)
  
  names(rGR) <- paste(as.vector(seqnames(rGR)),":",start(rGR),",",as.vector(strand(rGR)),sep="")
  #rGR
  rGR$name <- names(rGR)
  rGR$pos <- NA
  rGR[names(random_round)]$pos <- random_round
  rGR <- rGR[!is.na(rGR$pos)]
  #rGR
  
  match_list <- list()
  for(i in 1:length(dr_list_break))
  {
    match_list[[i]] <- rep(NA,length(dr_list_break[[i]]))
    names(match_list[[i]]) <- names(dr_list_break[[i]])
  }
  
  for(mynum in as.vector(seq(0,3,by=mdb)))
  {
    random_names <- names(which(random_round==mynum))
    random_names <- random_names[!is.na(random_names)]
    for(r in 1:length(match_list))
    {
      inds <- which(dr_list_break[[r]]==mynum)
      match_list[[r]][inds] <- sample(random_names,length(inds),replace=T)
    }
  }
  
  for(r in 1:length(match_list))
  {
    tGR[[r]]$match <- NA
    tGR[[r]][names(match_list[[r]])]$match <- match_list[[r]]
    tGR[[r]]$pos <- NA
    tGR[[r]][names(match_list[[r]])]$pos <- dr_list_break[[r]]
    tGR[[r]] <- tGR[[r]][!is.na(tGR[[r]]$pos)]
  }

  names(match_list) <- names(tGR)
  tGR[["random"]] <- rGR
  
 return(tGR) 
}


getCountDens <- function(matchGR,twrGR,tol=5)
{
  fgbg_list <- list()
  pcom <- rep(0,length(matchGR)-1)
  glmdat_list <- list()
  for(k in 1:(length(matchGR)-1))
  {
    ol1 <- countOverlaps(matchGR[[k]]+tol,twrGR+2)
    ol1r <- countOverlaps(rGR[matchGR[[k]]$match]+tol,twrGR+2)
    ol1[ol1>0] <- 1
    ol1r[ol1r>0] <- 1
    ol1sum <- tapply(ol1,matchGR[[k]]$pos,sum)
    ol1rsum <- tapply(ol1r,matchGR[[k]]$pos,sum)
    fgbg_list[[names(matchGR)[k]]] <- cbind(ol1sum,ol1rsum)
    glmdf <- data.frame(c(ol1,ol1r),c(rep("fg",length(ol1)),rep("bg",length(ol1r))))
    colnames(glmdf) <- c("count","condition")
    myglm <- (glm(as.numeric(condition=="fg")~count,data=glmdf))
    pcom[k] <- (summary(myglm)$coefficients)[2,4]
    glmdat_list[[k]] <- glmdf
  }
  return(list("fgbg"=fgbg_list,"sigs"=pcom,"glmdat"=glmdat_list))
}


plot_ol <- function(regsGR,nanoGRp,myn=100,para=F)
{
  nnano <- list()
  if(para)
  {
    nnano <- foreach(i=c(-myn:myn)) %dopar%
    {
      mshiftGR <- shift(regsGR,i) 
      gene_overlap <- which((countOverlaps(mshiftGR,gtfE))>0)
      mol <- sum(countOverlaps(mshiftGR[gene_overlap],nanoGRp))
      tlen <- length(gene_overlap)
      1000*mol/tlen
    }
  }else{
    for(i in -myn:myn)
    {
      mshiftGR <- shift(regsGR,i) 
      gene_overlap <- which((countOverlaps(mshiftGR,gtfE))>0)
      mol <- sum(countOverlaps(mshiftGR[gene_overlap],nanoGRp))
      tlen <- length(gene_overlap)
      nnano[[paste(i)]] <- 1000*mol/tlen
    }
  }

  return(unlist(nnano))
}

makeRelPlot <- function(setGR1,setGR2,myn=50,to.add=F,my.col="black",my.ylim=c(0,50),my.title="",to.smooth=T,use.para=F,to.norm=F)
{
  regsGRp <- setGR1[strand(setGR1)=="+"]
  nanoPLUSp <- setGR2[strand(setGR2)=="+"]
  
  regsGRm <- setGR1[strand(setGR1)=="-"]
  nanoPLUSm <- setGR2[strand(setGR2)=="-"]
  
  np <- plot_ol(nanoPLUSp,regsGRp,myn,para=use.para)
  nm <- plot_ol(nanoPLUSm,regsGRm,myn,para=use.para)
  
  nn <- np + rev(nm)
  
  
  wins <- myn
  if(!to.add)
  {
    #tp <- 1000*nn/length(setGR2)
    tp <- nn
    if(to.smooth)
    {
      plot(smooth.spline(tp),type="l",ylim=my.ylim,xaxt="n",ylab="Site per 1000",col=my.col,xlab="Distance from m6A",main=my.title)
    }else{
      if(to.norm)
      {
        tp <- tp/sum(tp)
        plot((tp),type="l",ylim=my.ylim,xaxt="n",ylab="Normalised sites per 1000",col=my.col,xlab="Distance from m6A",main=my.title)
        
      }else{
        plot((tp),type="l",ylim=my.ylim,xaxt="n",ylab="Site per 1000",col=my.col,xlab="Distance from m6A",main=my.title)
        
      }
    }
    #abline(v=wins+1+10,lty=2)
    #abline(v=wins+1-10,lty=2)
    abline(v=wins+1,lty=2)
    axis(1,c(1:length(tp))[c(1,wins+1,length(tp))],c(-wins:wins)[c(1,wins+1,length(tp))],las=2)
    
  }else{
    #tp <- 1000*nn/length(setGR2)
    tp <- nn
    if(to.smooth)
    {
      points(smooth.spline(tp),type="l",xaxt="n",ylab="Site per 1000",col=my.col,xlab="Distance from m6A",main="iCLIP overlap relative to m6A")
    }else{
      if(to.norm)
      {
        tp <- tp/sum(tp)
        points(tp,type="l",xaxt="n",ylab="Site per 1000",col=my.col,xlab="Distance from m6A",main="iCLIP overlap relative to m6A")
        
      }else{
        points(tp,type="l",xaxt="n",ylab="Site per 1000",col=my.col,xlab="Distance from m6A",main="iCLIP overlap relative to m6A")
        
      }
    }
  }
  
}


makeDBplot <- function(aSet_list,motGR,colset,inc.bar=F,y.min=0,y.max=50,my.name="WEI",my.set.names,my.size=50,to.smooth=F)
{
  makeRelPlot(motGR,aSet_list[[1]],to.add=F,my.col=colset[1],my.ylim=c(y.min,y.max),my.title=my.name,myn=my.size,to.smooth = to.smooth)
  legend("topright",my.set.names,col=colset[c(1:length(aSet_list))],lwd=3)
  for(k in 2:length(aSet_list))
  {
    makeRelPlot(motGR,aSet_list[[k]],to.add=T,my.col=colset[k],my.ylim=c(y.min,y.max),myn=my.size,to.smooth = to.smooth)
  }
  
  if(inc.bar & (length(aSet_list)>1))
  {
    co_list <- list()
    for(k in 1:length(aSet_list))
    {
      tmp <- (countOverlaps(aSet_list[[k]]+25,motGR))
      co_list[[k]] <- 1000*sum(tmp)/length(aSet_list[[k]])
    }
    tp <- unlist(co_list)
    names(tp) <- my.set.names
    barplot(tp,col=colset[1:length(aSet_list)],ylab="Sites per 1000 m6A",main="m6A +/- 25bp motif count",las=2)
  }
}




makeRelPlotSets <- function(setGR1,setGR2_list,myn=50,my.col="black",my.ylim=NULL,my.title="",use.para=F)
{
  regsGRp <- setGR1[strand(setGR1)=="+"]
  regsGRm <- setGR1[strand(setGR1)=="-"]
  
  tp_list <- list()
  for(k in 1:length(setGR2_list))
  {
    setGR2 <- setGR2_list[[k]]
    nanoPLUSp <- setGR2[strand(setGR2)=="+"]
    nanoPLUSm <- setGR2[strand(setGR2)=="-"]
    
    np <- plot_ol(nanoPLUSp,regsGRp,myn,para=use.para)
    nm <- plot_ol(nanoPLUSm,regsGRm,myn,para=use.para)
    
    nn <- np + rev(nm)
    tp_list[[k]] <- nn
  }

  wins <- myn

  
  if(!is.null(my.ylim))
  {
    my.min <- my.ylim[1]
    my.max <- my.ylim[2]
  }else{
    my.max <- max(unlist(tp_list))
    my.min <- min(unlist(tp_list))
  }
  
  for(i in 1:length(tp_list))
  {
    tp <- tp_list[[i]]
    if(i==1)
    {
      plot((tp),type="l",ylim=c(my.min,my.max),xaxt="n",ylab="Site per 1000",col=my.col[i],xlab="Distance from m6A",main=my.title)
      
    }else{
      points((tp),type="l",col=my.col[i])
    }
  }
  
  abline(v=wins+1,lty=2)
  axis(1,c(1:length(tp))[c(1,wins+1,length(tp))],c(-wins:wins)[c(1,wins+1,length(tp))],las=2)

  
}




########################################################
#####HELPER FUNCTIONS FOR RANDOM FOREST ANALYSIS########
########################################################

createMerged <- function(xGR)
{
  ngrf <- xGR[as.vector(strand(xGR)) %in% c("+","-")]
  ngred <- reduce(ngrf+25)
  newlocs <- round(start(ngred)+(width(ngred)-1)/2)
  start(ngred) <- newlocs
  end(ngred) <- newlocs
  ngrf_merged <- ngrf[nearest(ngred,ngrf)]
  return(ngrf_merged)
}

makeGRF <- function(X,resp)
{
  X[X>10] <- 10
  dim(X)
  X <- apply(X,2,function(x) scale(x)[,1])
  
  my.sets <- sample(1:5,nrow(X),replace=T)
  X_test <- X[which((my.sets==5)),]
  resp_test <- resp[which((my.sets==5))]
  dim(X_test)
  summary(as.factor(resp_test))
  
  X_train <- X[which(!(my.sets==5)),]
  resp_train <- resp[which(!(my.sets==5))]
  dim(X_train)
  summary(as.factor(resp_train))
  
  #library(gbm)
  boost <- gbm(as.factor(resp_train)~.,
               data=as.data.frame(X_train),
               distribution="gaussian",
               shrinkage=0.05,interaction.depth=6,
               cv.folds = 5,
               n.trees=c(2000))
  return(boost)
}

getRFtest <- function(resp,X)
{
  X[X>10] <- 10
  dim(X)
  X <- apply(X,2,function(x) scale(x)[,1])
  
  my.sets <- sample(1:5,nrow(X),replace=T)
  X_test <- X[which((my.sets==5)),]
  resp_test <- resp[which((my.sets==5))]
  return(list("resp_test"=resp_test,"X_test"=X_test))
}


makeSummaryPlot <- function(boost,X_test)
{
  bres <- summary(boost,plot=F)
  preds <- predict.gbm(boost,as.data.frame(X_test),2000)
  #hist(preds-1)
  
  #require(pROC)
  rf.roc<-roc(as.factor(resp_test),preds-1)
  
  #pdf("random_forest_motifs_results_gradient_boost.pdf",width=10,height=5)
  #par(mfrow=c(1,2))
  
  plot(rf.roc)
  text(0.5,0.5,paste("AUC=",round(auc(rf.roc),4)))
  
  
  tp <- bres$rel.inf[1:30]
  names(tp) <- bres$var[1:30]
  mycols <- brewer.pal(3,"Set1")
  colvec <- rep(mycols[1],length(tp))
  colvec[grep("up",names(tp))] <- mycols[2]
  colvec[grep("down",names(tp))] <- mycols[3]
  
  barplot(tp,las=2,cex.names = 0.75,col = colvec)
  #dev.off()
}



makeRFdata <- function(add_size=25,center_size=10,nanoGRme,nanoMatchGR)
{
  shift_size <- add_size+center_size
  
  nanoGRme_US <- nanoGRme
  nanoGRme_US[strand(nanoGRme_US)=="+"] <- shift(nanoGRme_US[strand(nanoGRme_US)=="+"],-shift_size) + add_size
  nanoGRme_US[strand(nanoGRme_US)=="-"] <- shift(nanoGRme_US[strand(nanoGRme_US)=="-"],+shift_size) + add_size
  
  nanoGRme_DS <- nanoGRme
  nanoGRme_DS[strand(nanoGRme_DS)=="+"] <- shift(nanoGRme_DS[strand(nanoGRme_DS)=="+"],+shift_size) + add_size
  nanoGRme_DS[strand(nanoGRme_DS)=="-"] <- shift(nanoGRme_DS[strand(nanoGRme_DS)=="-"],-shift_size) + add_size
  
  
  nanoGRmeM_US <- nanoMatchGR
  nanoGRmeM_US[strand(nanoGRmeM_US)=="+"] <- shift(nanoGRmeM_US[strand(nanoGRmeM_US)=="+"],-shift_size) + add_size
  nanoGRmeM_US[strand(nanoGRmeM_US)=="-"] <- shift(nanoGRmeM_US[strand(nanoGRmeM_US)=="-"],+shift_size) + add_size
  
  nanoGRmeM_DS <- nanoMatchGR
  nanoGRmeM_DS[strand(nanoGRmeM_DS)=="+"] <- shift(nanoGRmeM_DS[strand(nanoGRmeM_DS)=="+"],+shift_size) + add_size
  nanoGRmeM_DS[strand(nanoGRmeM_DS)=="-"] <- shift(nanoGRmeM_DS[strand(nanoGRmeM_DS)=="-"],-shift_size) + add_size
  
  nanoGRme_WI <- nanoGRme+center_size
  nanoGRmeM_WI <- nanoMatchGR+center_size
  
  nol <- list()
  nol_bg <- list()
  
  nol_us <- list()
  nol_ds <- list()
  nol_wi <- list()
  
  nol_bg_us <- list()
  nol_bg_ds <- list()
  nol_bg_wi <- list()
  
  for(i in 1:length(motGR))
  {
    nol[[i]] <- countOverlaps(nanoGRme+50,motGR[[i]])
    nol_bg[[i]] <- countOverlaps(nanoMatchGR+50,motGR[[i]])
    
    nol_us[[i]] <- countOverlaps(nanoGRme_US,motGR[[i]])
    nol_ds[[i]] <- countOverlaps(nanoGRme_DS,motGR[[i]])
    nol_wi[[i]] <- countOverlaps(nanoGRme_WI,motGR[[i]])
    
    nol_bg_us[[i]] <- countOverlaps(nanoGRmeM_US,motGR[[i]])
    nol_bg_ds[[i]] <- countOverlaps(nanoGRmeM_DS,motGR[[i]])
    nol_bg_wi[[i]] <- countOverlaps(nanoGRmeM_WI,motGR[[i]])
  }
  
  nol <- do.call(cbind,nol)
  nol_bg <- do.call(cbind,nol_bg)
  
  colnames(nol) <- names(motGR)
  colnames(nol_bg) <- names(motGR)
  
  nol_us <- do.call(cbind,nol_us)
  nol_ds <- do.call(cbind,nol_ds)
  nol_wi <- do.call(cbind,nol_wi)
  
  nol_bg_us <- do.call(cbind,nol_bg_us)
  nol_bg_ds <- do.call(cbind,nol_bg_ds)
  nol_bg_wi <- do.call(cbind,nol_bg_wi)
  
  colnames(nol_us) <- paste0(names(motGR),"_up")
  colnames(nol_ds) <- paste0(names(motGR),"_down")
  colnames(nol_wi) <- paste0(names(motGR),"_at")
  
  colnames(nol_bg_us) <- paste0(names(motGR),"_up")
  colnames(nol_bg_ds) <- paste0(names(motGR),"_down")
  colnames(nol_bg_wi) <- paste0(names(motGR),"_at")
  
  resp <- c(rep(1,nrow(nol)),rep(0,nrow(nol_bg)))
  
  #X <- rbind(nol,nol_bg)
  X <- rbind(cbind(nol_us,nol_ds,nol_wi),cbind(nol_bg_us,nol_bg_ds,nol_bg_wi))
  #colnames(X)
  
  return(list("resp"=resp,"X"=X))
}




