#' Compute the trimmed mean, effective sample size, and squared standard error.
#'  The default amount of trimming is tr=.2.
#
trimparts <- function(x,tr=.2){
  tm<-mean(x,tr)
  h1<-length(x)-2*floor(tr*length(x))
  sqse<-(length(x)-1)*winvar(x,tr)/(h1*(h1-1))
  trimparts<-c(tm,sqse)
  trimparts
}

#'  Compute the Winsorized variance for the data in the vector x.
#'  tr is the amount of Winsorization which defaults to .2.
#'
winvar <- function(x,tr=.2,na.rm=FALSE,STAND=NULL){
  remx=x
  x<-x[!is.na(x)]
  y<-sort(x)
  n<-length(x)
  ibot<-floor(tr*n)+1
  itop<-length(x)-ibot+1
  xbot<-y[ibot]
  xtop<-y[itop]
  y<-ifelse(y<=xbot,xbot,y)
  y<-ifelse(y>=xtop,xtop,y)
  wv<-var(y)
  if(!na.rm)if(sum(is.na(remx)>0))wv=NA
  wv
}

  #'  Estimate the standard error of the trimmed mean
  #'  The default amount of trimming is tr=.2.
  #'
trimse<-function(x,tr=.2,na.rm=FALSE){
  if(na.rm)x<-x[!is.na(x)]
  trimse<-sqrt(winvar(x,tr))/((1-2*tr)*sqrt(length(x)))
  trimse
}