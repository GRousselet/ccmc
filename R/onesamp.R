#' One sample percentile bootstrap test with cluster correction for multiple comparisons
#' 
#' 
#'
onesampb.ccmc <- function(x, nboot=599){
   
}

#' One sample t-test on trimmed means with cluster correction for multiple 
#' comparisons
#' 
#' The t-test is done using Yuen's extension of the t-test to trimmed means. Set
#' the argument tr=0 to use means instead of trimmed means. The cluster 
#' correction is achieved using a bootstrap-t approach: each marginal 
#' distribution is centred, so that the null hypothesis is true, then samples
#' with replacement are taken. For each of these bootstrap samples, t-tests are computed, 
#' providing data driven t distributions under the null hypothesis.
#' Clusters are formed by proximity, so that if 2 contiguous t-tests are p < alpha,
#' their values are included in the cluster calculation. The cluster statistics in the 
#' sum of the squared t-values. 
#' 
#' @param x A list or a matrix. In the first case x[[1]] contains the data for 
#'   the first group, x[[2]] the data for the second group, etc. Length(x) = the
#'   number of groups = J. If stored in a matrix, columns correspond to groups.
#' @param nullval The null value against which to compare the trimmed mean - 
#'   default = 0.
#' @param tr The amount of trimming - default 0.2. Set tr=0 for a t-test on 
#'   means.
#' @param alpha Alpha level - default 0.05.
#' @param bt Set to TRUE to compute critical t values, p values and confidence 
#'   intervals based on boostrap-t distributions, otherwise use standard 
#'   calculations - default = TRUE
#' @param nboot Number of bootstrap samples - default = 599
#' @return A list of univariate and cluster-based statistics:
#'   
#'   estimate = mean or trimmed mean
#'   
#'   ci = confidence interval
#'   
#'   tval = t values
#'   
#'   pval = p values
#'   
#'   cluster.th = cluster threshold
#'   
#'   cluster.map = vector of cluster IDs
#'   
#'   cluster.sig = statistical significance based on cluster test
#' @examples 
#' # no effect: 
#' set.seed(21)
#' x <- matrix(rnorm(100), ncol = 5)
#' trimbt.ccmc(x,nullval=0,tr=0,alpha=.05,bt=FALSE,nboot=599)
#' 
#' # cluster of length 3:
#' set.seed(21)
#' x <- matrix(rnorm(100), ncol = 5)
#' x[,3:5] <- x[,3:5] + 1
#' trimbt.ccmc(x,nullval=0,tr=0,alpha=.05,bt=FALSE,nboot=599)
#' 
#' # use bootstrap-t thresholds
#' trimbt.ccmc(x,nullval=0,tr=0,alpha=.05,bt=TRUE,nboot=599)
#' 
#' # get cluster statistics for cluster 1:
#' out <- trimbt.ccmc(x,nullval=0,tr=0,alpha=.05,bt=TRUE,nboot=599)
#' c2sum <- sum(out$tval[out$cluster.map==1]^2)
#' @section References
#' @seealso \code{trimbt.tfce} 
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
  # print("Taking bootstrap samples. Please wait.")
  for(j in 1:J){
    hval[j] <- length(x[[j]])-2*floor(tr*length(x[[j]]))
    # hval is the number of observations in the jth group after trimming.
    # print(paste("Working on group",j))
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
       cluster.th = boot.th,
       cluster.map = cmap,
       cluster.sig = cluster.sig)
}

#' One sample t-test on trimmed means with TFCE correction for multiple 
#' comparisons
#'  
#' The t-test is done using Yuen's extension of the t-test to trimmed means. Set
#' the argument tr=0 to use means instead of trimmed means. The cluster 
#' correction is achieved using a bootstrap-t approach: each marginal 
#' distribution is centred, so that the null hypothesis is true, then samples
#' with replacement are taken. For each of these bootstrap samples, t-tests are computed, 
#' providing data driven t distributions under the null hypothesis.
#' Clusters are formed by proximity, but without using p values as thresholds as in 
#' \code{trimbt.ccmc}. 
#' Instead, multiple thresholds are used, exploring the range of the data, and cluster
#' statistics are formed by integrating values across multiple thresholds, following
#' the \code{tfce} equation. To correct for multiple comparisons, the observed t values
#' are squared (F values) and then TFCE transformed. These values are then compared to the
#' (1-alpha) percentile of the distribution of max bootstrap TFCE scores. 
#' 
#' @param x A list or a matrix. In the first case x[[1]] contains the data for 
#'   the first group, x[[2]] the data for the second group, etc. Length(x) = the
#'   number of groups = J. If stored in a matrix, columns correspond to groups.
#' @param nullval The null value against which to compare the trimmed mean - 
#'   default = 0.
#' @param tr The amount of trimming - default 0.2. Set tr=0 for a t-test on 
#'   means.
#' @param alpha Alpha level - default 0.05.
#' @param nboot Number of bootstrap samples - default = 599
#' @return A list of TFCE results:
#'
#'   tfce.score = TFCE transformation of the observed squared t values    
#'   
#'   tfce.th = TFCE threshold
#'   
#'   tfce.sig = statistical significance based on the TFCE test
#'   
#' @examples 
#' # no effect: 
#' set.seed(21)
#' x <- matrix(rnorm(100), ncol = 5)
#' trimbt.ccmc(x,nullval=0,tr=0,alpha=.05,bt=FALSE,nboot=599)
#' 
#' # cluster of length 3:
#' set.seed(21)
#' x <- matrix(rnorm(100), ncol = 5)
#' x[,3:5] <- x[,3:5] + 1
#' trimbt.ccmc(x,nullval=0,tr=0,alpha=.05,bt=FALSE,nboot=599)
#' 
#' # use bootstrap-t thresholds
#' trimbt.ccmc(x,nullval=0,tr=0,alpha=.05,bt=TRUE,nboot=599)
#' 
#' # get cluster statistics for cluster 1:
#' out <- trimbt.ccmc(x,nullval=0,tr=0,alpha=.05,bt=TRUE,nboot=599)
#' c2sum <- sum(out$tval[out$cluster.map==1]^2)
#' @section References
#' @seealso \code{tfce} \code{trimbt.ccmc}
#' @include
trimbt.tfce <- function(x,nullval=0,tr=.2,alpha=.05,nboot=599,E=2, H=1, dh=0.1){
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
  tfce.score <- tfce(tval^2,E=E, H=H, dh=dh)
  bootsamp <- array(0,c(J,2,nboot)) # declare matrix of bootstrap samples
  # hval <- vector("numeric",J)
  # print("Taking bootstrap samples. Please wait.")
  for(j in 1:J){
    # hval[j] <- length(x[[j]])-2*floor(tr*length(x[[j]]))
    # hval is the number of observations in the jth group after trimming.
    # print(paste("Working on group",j))
    xcen <- x[[j]] - mean(x[[j]],tr) # centre distribution
    data <- matrix(sample(xcen,size=length(x[[j]])*nboot,replace=TRUE),nrow=nboot)
    bootsamp[j,,] <- apply(data,1,trimparts,tr) # A 2 by nboot matrix. The first row
    #                     contains the bootstrap trimmed means, the second row
    #                     contains the bootstrap squared standard errors.
  }
  # bootsamp[,1,] = J by nboot matrix containing the bootstrap trimmed means
  # bootsamp[,2,] = J by nboot matrix containing the bootstrap sq standard errors
  boot.tval <- (bootsamp[,1,]-nullval) / sqrt(bootsamp[,2,]) # bootstrap t values
  # TFCE test ====================
  tfce.boot <- tfce(boot.tval^2,E=E, H=H, dh=dh)
  tfce.max <- sort(apply(tfce.boot, 2, max))
  icrit <- round((1-alpha)*nboot) # 95th quantile
  tfce.sig <- tfce.score >= tfce.max[icrit]
  # outputs
  list(tfce.score = tfce.score,
       tfce.th = tfce.max[icrit],
       tfce.sig = tfce.sig)
}