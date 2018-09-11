#' Make clusters based on binary maps
cluster.make <- function(x){
  y <- rle(x)
  cmap <- vector(mode = "numeric", length = 0)
  nC <- length(y$values) # number of clusters
  indx <- 0 # cluster counter
  for(CL in 1:nC){
    if(y$values[CL] == 0){
      val <- 0
    } else {
      indx <- indx + 1
      val <- indx
    }
    cmap <- c(cmap, rep(val, y$lengths[CL]))
  }
  cmap
}

#' Save sum for each cluster
cluster.sum <- function(values, cmap){
  csum <- vector(mode = "numeric", length = max(cmap))
  if(max(cmap)>0){
    for(CL in 1:max(cmap)){
      csum[CL] <- sum(values[cmap==CL])
    }
  } else {
    csum <- 0
  }
  csum
}

cluster.test <- function(values, cmap, boot.th){
  csig <- vector(mode = "logical", length = length(cmap))
  if(max(cmap)>0){
    for(CL in 1:max(cmap)){
      csig[cmap==CL] <- sum(values[cmap==CL]) > boot.th
    }
  } else {
    csig <- FALSE
  }
  csig
}

#' Threshold Free Cluster Enhancement
#'  
#' Implementation of the Threshold-free cluster enhancement method
#' developped for fMRI by Smith & Nichols, NeuroImage 44(2009), 83-98
#' tfce = sum(extent(h)^E*height^H*dh)
#' R adaptation of Matlab code from the LIMO EEG toolbox:
#' https://github.com/LIMO-EEG-Toolbox/limo_eeg/blob/master/limo_tfce.m
#' @param data A vector or matrix of positive scores or statistics. 
#' If matrix, the columns are bootstrap samples.
#' @param E Value of extend parameter - default = 2
#' @param H Value of height parameter - default = 1
#' @param dh Value of integration step parameter - default = 0.1
#' @return A map (vector or matrix) of TFCE scores.
#' @examples
#' # Vector
#' # Matrix
#' @section References
#' Pernet, C., Latinus, M., Nichols, T.E., & Rousselet, G.A. (2015)
#' Cluster-based computational methods for mass univariate analyses
#' of event-related brain potentials/fields: a simulation study
#' Journal Of Neuroscience Method 250, Pages 85–93
#' <10.1016/j.jneumeth.2014.08.003>
#'
#' Pernet, Cyril; Rousselet, Guillaume (2014): Type 1 error rate using TFCE for ERP.
#' figshare. http://dx.doi.org/10.6084/m9.figshare.1008325
#' @export
tfce <- function(data, E = 2, H = 1, dh = 0.1){
  
  if(min(data) < 0){
    error("tfce() only accepts data with positive values")
  }
  
  if(is.vector(data)){ # vector
    
    min.data <- min(data)
    max.data <- max(data)
    length.data <- length(data)
    
    # define increment size forced by dh
    data_range <- max(data) - min(data)
    if(data_range > 1){
      precision <- round(data_range / dh)
      # arbitrary decision to limit precision to 200th of the data range - 
      # needed as sometime under H0 one value can be very wrong
      if(precision > 200){ 
        increment <- data_range / 200
      } else {
        increment <- data_range / precision
      }
    } else {
      increment <- data_range * dh
    }
    
    # select a height, obtain cluster map, obtain extent map
    # = cluster map but with extent of cluster rather than number of the cluster
    # then tfce score for that height
    index <- 1
    nsteps <- length(seq(from=min.data, to=max.data, by=increment))
    tfce <- matrix(NA, nrow = length.data, ncol = nsteps)
    for(h in seq(from=min.data, to=max.data, by=increment)){
      clustered_map <- cluster.make(data > h)
      if(max(clustered_map) == 0){
        tfce[,index] <- 0
        index <- index + 1
      } else {
        extent_map <- vector(mode = "numeric", length = length.data) # same as cluster map but contains extent value instead
        for(i in 1:max(clustered_map)){
          idx <- clustered_map == i
          extent_map[idx] <- sum(idx)
        }
        tfce[,index] <- (extent_map^E)*h^H*increment
        index <- index + 1
      }
    }
    # compute final score
    tfce_score <- apply(tfce, 1, sum, na.rm = TRUE)
  } else if(is.matrix(data)){ # matrix of bootstrap samples
    out <- dim(data)
    length.data <- out[[1]]
    nboot <- out[[2]]
    tfce_score <- matrix(NA, nrow = length.data, nboot)
    
    # select a height, obtain cluster map, obtain extent map
    # then tfce score for that height
    
    for(boot in 1:nboot){
      tmp.data <- data[,boot]
      min.data <- min(tmp.data)
      max.data <- max(tmp.data)
      # define increment size forced by dh
      data_range <- max(tmp.data) - min(tmp.data)
      if(data_range > 1){
        precision <- round(data_range / dh)
        if(precision > 200){
          increment <- data_range / 200
        } else {
          increment <- data_range / precision
        }
      } else {
        increment <- data_range * dh
      }
      
      index <- 1
      nsteps <- length(seq(from=min.data, to=max.data, by=increment))
      tfce <- matrix(NA, nrow = length.data, ncol = nsteps)
      
      for(h in seq(from=min.data, to=max.data, by=increment)){
        clustered_map <- cluster.make(tmp.data > h)
        if(max(clustered_map) == 0){
          tfce[,index] <- 0
          index <- index + 1
        } else {
        extent_map <- vector(mode = "numeric", length = length.data)
        for(i in 1:max(clustered_map)){
          idx <- clustered_map == i
          extent_map[idx] <- sum(idx)
        }
        tfce[,index] <- (extent_map^E)*h^H*increment
        index <- index + 1
        }
      }
      tfce_score[, boot] <- apply(tfce, 1, sum, na.rm = TRUE)
    }
  } else {
    error("data must be a vector or a matrix")
  }
  tfce_score
}
  
  