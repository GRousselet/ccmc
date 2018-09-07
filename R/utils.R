#' Remove any rows of data having missing values
#'
elimna<-function(m){
  DONE=FALSE
  if(is.list(m) && is.matrix(m)){
    z=pool.a.list(m)
    m=matrix(z,ncol=ncol(m))
    DONE=TRUE
  }
  if(!DONE){
    if(is.list(m) && is.matrix(m[[1]])){
      for(j in 1:length(m))m[[j]]=na.omit(m[[j]])
      e=m
      DONE=TRUE
    }}
  if(!DONE){
    if(is.list(m) && is.null(dim(m))){ #!is.matrix(m))
      for(j in 1:length(m))m[[j]]=as.vector(na.omit(m[[j]]))
      e=m
      DONE=TRUE
    }}
  if(!DONE){
    #if(!is.list(m)){
    #if(is.null(dim(m)))
    m<-as.matrix(m)
    ikeep<-c(1:nrow(m))
    for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
    e<-m[ikeep[ikeep>=1],]
    #}
  }
  e
}

listm<-function(x){
  #
  # Store the data in a matrix or data frame in a new
  # R variable having list mode.
  # Col 1 will be stored in y[[1]], col 2 in y[[2]], and so on.
  #
  if(is.null(dim(x)))stop("The argument x must be a matrix or data frame")
  y<-list()
  for(j in 1:ncol(x))y[[j]]<-x[,j]
  y
}