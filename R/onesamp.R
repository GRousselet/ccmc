#' One sample percentile bootstrap test with cluster correction for multiple comparisons
#' 
#' 
#'
# onesampb.ccmc <- function(x, nboot=599){
#   
# }

#' One sample t-test on trimmed means with cluster correction for multiple comparisons
#   
#
#   grp is used to specify some subset of the groups, if desired.
#   By default, all J groups are used.
#   g=NULL, x is assumed to be a matrix or have list mode
#
#   if g is specifed, it is assumed that column g of x is
#   a factor variable and that the dependent variable of interest is in column
#   dp of x, which can be a matrix or data frame.
#
#   The default number of bootstrap samples is nboot=599
#' @param x A list or a matrix. In the first case
#'   x[[1]] contains the data for the first group, x[[2]] the data
#'   for the second group, etc. Length(x) = the number of groups = J.
#'   If stored in a matrix, columns correspond to groups.
#' @param nullval The null value against which to compare the trimmed mean - default = 0.
#' @param tr The amount of trimming - default 0.2. Set tr=0 for a t-test on means. 
#' @param alpha Alpha level - default 0.05.
#' @param bt Set to TRUE to compute critical t values, p values and confidence intervals based on boostrap-t distributions, otherwise use standard calculations - default = TRUE
#' @param nboot Number of bootstrap samples - default = 599
#' @return A list of univariate t values, p values and confidence intervals, as well as cluster-based statistics.
#' @examples 
#' # no effect: 
#' set.seed(21)
#' x <- matrix(rnorm(100), ncol = 5)
#' # cluster of length 3:
#' set.seed(21)
#' x <- matrix(rnorm(100), ncol = 5)
#' x[,3:5] <- x[,3:5] + 1
#' @section References
#' @seealso \code{}
#' @include
trimbt.ccmc <- function(x,nullval=0,tr=.2,alpha=.05,bt=TRUE,nboot=599){
  if(is.data.frame(x)){
    x=as.matrix(x)
  }
  if(is.matrix(x)){
    x <- listm(x)
  }
  if(!is.list(x)){
    stop("Data must be stored in a matrix or in list mode.")
  }
  J <- length(x)
  x <- lapply(x,elimna) # remove missing values
  m <- unlist(lapply(x,mean,tr=tr))
  se <- unlist(lapply(x,trimse,tr=tr))
  tval <- (m - nullval) / se # t-test on original data
  bootsamp <- array(0,c(J,2,nboot)) # declare matrix of bootstrap samples
  hval <- vector("numeric",J)
  print("Taking bootstrap samples. Please wait.")
  for(j in 1:J){
    hval[j] <- length(x[[j]])-2*floor(tr*length(x[[j]]))
    # hval is the number of observations in the jth group after trimming.
    print(paste("Working on group",j))
    xcen <- x[[j]] - mean(x[[j]],tr) # centre distribution
    data <- matrix(sample(xcen,size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
    bootsamp[j,,] <- apply(data,1,trimparts,tr) # A 2 by nboot matrix. The first row
    #                     contains the bootstrap trimmed means, the second row
    #                     contains the bootstrap squared standard errors.
  }
  # bootsamp[,1,] = J by nboot matrix containing the bootstrap trimmed means
  # bootsamp[,2,] = J by nboot matrix containing the bootstrap sq standard errors
  boot.tval <- (bootsamp[,1,]-nullval) / sqrt(bootsamp[,2,]) # bootstrap t values
  boot.tval <- apply(abs(boot.tval), 2, sort)
  icrit <- round((1-alpha)*nboot) # 95th quantile
  if(bt){ # use bootstrap-t thresholds
    # confidence intervals
    ci <- matrix(0, nrow = J, ncol = 2)
    tcrit <- boot.tval[,icrit]
    ci[,1] <- m - tcrit * se
    ci[,2] <- m + tcrit * se
    # p values  
    pval <- vector(mode = "numeric", length = J)
    for(j in 1:J){
      pval[j] <- ( sum(abs(tval[j]) <= abs(boot.tval[j,])) ) / nboot
    }
    # cluster test ====================
    cmap <- cluster.make(pval <= alpha) # form clusters in original data
    boot.cmap <- t(apply(boot.tval > matrix(rep(tcrit,nboot),nrow = J), 1, cluster.make))
    boot.max.sums <- vector(mode = "numeric", length = J)
    for(B in 1:nboot){ # max cluster sum for each bootstrap sample
      boot.max.sums[B] <- max(cluster.sum(values = boot.tval[,B]^2, cmap = boot.cmap[,B]))
    }
    # cluster sum threshold
    boot.th <- sort(boot.max.sums)[icrit]
    # cluster significance
    cluster.sig <- cluster.test(values = tval^2, cmap = cmap, boot.th)
  } else { # use standard calculations
    # confidence intervals and p values
    df <- hval-1
    ci <- matrix(0, nrow = J, ncol = 2)
    tcrit <- qt(1-alpha/2,df)
    ci[,1] <- m - tcrit * se
    ci[,2] <- m + tcrit * se
    pval <- 2*(1-pt(abs(tval),df))
    # cluster test ====================
    cmap <- cluster.make(pval <= alpha) # form clusters in original data
    boot.cmap <- t(apply(boot.tval > matrix(rep(tcrit,nboot),nrow = J), 1, cluster.make))
    boot.max.sums <- vector(mode = "numeric", length = J)
    for(B in 1:nboot){ # max cluster sum for each bootstrap sample
      boot.max.sums[B] <- max(cluster.sum(values = boot.tval[,B]^2, cmap = boot.cmap[,B]))
    }
    # cluster sum threshold
    boot.th <- sort(boot.max.sums)[icrit]
    # cluster significance
    cluster.sig <- cluster.test(values = tval^2, cmap = cmap, boot.th)
  }
  # outputs
  list(estimate = m,
       ci = ci,
       tval = tval,
       pval = pval,
       cluster.map = cmap,
       cluster.sig = cluster.sig)
}