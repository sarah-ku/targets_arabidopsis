
returnDens <- function(tGR,gtfGR){
  props_list_plus <- list()
  props_list_minus <- list()
  props_list <- list()
  for(type in c("5UTR","CDS","3UTR"))
  {
    cdsOL <- findOverlaps(tGR,gtfGR[gtfGR$type==type])
    ss <- start(gtfGR[gtfGR$type==type][subjectHits(cdsOL)])
    sq <- start(tGR[queryHits(cdsOL)])
    es <- end(gtfGR[gtfGR$type==type][subjectHits(cdsOL)])
    props_plus <- (sq-ss+1)/width(gtfGR[gtfGR$type==type][subjectHits(cdsOL)])
    names(props_plus) <- names(tGR[queryHits(cdsOL)])
    props_plus <- props_plus[which(strand(tGR[queryHits(cdsOL)])=="+")]
    
    props_minus <- (es - sq +1)/width(gtfGR[gtfGR$type==type][subjectHits(cdsOL)])
    props_minus <- props_minus[which(strand(tGR[queryHits(cdsOL)])=="-")]
    
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
    points(y=rtab,x=xls,col="grey",type="n",pch=16,xaxt="n",
           ylab="Normalised enrichment",xlab="Genomic position",
           ylim=c(mymin,mymax),lty=2,lwd=1)
  }else{
    plot(y=rtab,x=xls,col="grey",type="n",pch=16,xaxt="n",
         ylab="Normalised enrichment",xlab="Genomic position",
         ylim=c(mymin,mymax),lty=2,lwd=1,main=mytitle)
  }
  
  
  tp <- rtab
  x <- c(xls,rev(xls))
  
  if(inc_int)
  {
    if(to.smooth)
    {
      y <- c(smooth.spline(rtab_lower)$y, rev(smooth.spline(rtab_upper)$y))
      
      polygon(x,y, col = brewer.pal("Set1",n = 6)[col.no], border = NA)
      points(y=smooth.spline(rtab_mean)$y,x=xls,col="black",lty=1,type="l")
      #points(y=smooth.spline(rtab_lower)$y,x=xls,type="l",col="grey",lty=2)
      #points(y=smooth.spline(rtab_upper)$y,x=xls,type="l",col="grey",lty=2)
    }else{
      y <- c(rtab_lower, rev(rtab_upper))
      
      polygon(x,y, col = brewer.pal("Set1",n = 6)[col.no], border = NA)
      points(y=rtab,x=xls,col="black",lty=1,type="l")
      #points(y=rtab_lower,x=xls,type="l",col="grey",lty=2)
      #points(y=rtab_upper,x=xls,type="l",col="grey",lty=2)
    }
  }else{
    points(y=rtab,x=xls,col=brewer.pal("Set1",n = 6)[col.no],lty=1,lwd=2,type="l")
  }
  
  
  
  abline(h=1,lty=2)
  axis(1,c(0,1,2,3),c(0,1,2,3))
  abline(v=c(0,1,2,3),lty=2)
}

plotDens <- function(bs=0.005,dr_rand,dr,dr_ind,inc_groups=T,inc_int=T,to.add=F,col.no=1,no.samps=10)
{
  #bs <- 0.005
  #dr_r <- dr_rand[sample(1:length(dr_rand),length(dr))]
  rtab <- as.vector(table(cut(dr,breaks=seq(0,3,by=bs*3))))/as.vector(table(cut(dr_rand,breaks=seq(0,3,by=bs*3))))
  #  nl <- length(dr)/length(dr_rand)
  #  rtab <- (rtab/nl)
  
  #plot(density(dr),ylim=c(0,0.7))
  #points(density(dr_rand)$x,density(dr_rand)$y,type="l")
  #plot(rtab,type="l")
  #rtab <- rtab/mean(rtab[xls>1 & xls<2])
  #rtab <- smooth.spline(rtab)$y
  #plot(rtab,pch=16,col="blue",type="l")
  #points(rtab,type="l")
  
  rtab_rand <- list()
  for(i in 1:no.samps)
  {
    dr_r <- dr_rand[sample(1:length(dr_rand),length(dr_rand),replace=T)]
    
    #dr_r <- dr_rand[sample(1:length(dr_rand),length(dr),replace=T)]
    rtab_rand[[i]] <- as.vector(table(cut(dr,breaks=seq(0,3,by=bs*3))))/as.vector(table(cut(dr_r,breaks=seq(0,3,by=bs*3))))
    #nl <- length(dr)/length(dr_rand)
  }
  
  rtab_mat <- do.call(rbind,rtab_rand)
  rtab_upper <- apply(rtab_mat,2,function(x) quantile(x,.9))
  rtab_lower <- apply(rtab_mat,2,function(x) quantile(x,.1))
  #summary(rtab_lower)
  
  rtab_upper_SMOOTH <- smooth.spline(rtab_upper)$y
  rtab_lower_SMOOTH <- smooth.spline(rtab_lower)$y
  
  #points(rtab_upper,type="l")
  
  #points(rtab_lower,type="l")
  
  #if(inc_groups)
  #{
  rtabv <- list()
  for(i in 1:length(dr_ind))
  {
    rtabv[[i]] <- as.vector(table(cut(dr_ind[[i]],breaks=seq(0,3,by=bs))))/as.vector(table(cut(dr_rand,breaks=seq(0,3,by=bs))))
    nl <- length(dr_ind[[i]])/length(dr_rand)
    rtabv[[i]] <- rtabv[[i]]/nl
    #rtabv[[i]] <- rtabv[[i]]/mean(rtabv[[i]][xls>1 & xls<2])
    #rtabv[[i]] <- smooth.spline(rtabv[[i]])$y
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
    #mymin <- min(c(rtab_lower_SMOOTH))
    #mymax <- max(c(rtab_upper_SMOOTH))
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

