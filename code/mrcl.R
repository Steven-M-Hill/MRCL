mrcl = function(X, Gobs, nbins = 30, binsRange = c(-3,3), outOfRange = "ignore",
                q = 100, kernelType = "rbfdot", neighboursCount = 50, 
                batchSize = 100, lambda1 = 0.001, lambda2 = 0.001) {
  # Manifold Regularized Causal Learning (MRCL) as described in Hill, Oates,
  # Blythe & Mukherjee (2019) arXiv preprint arXiv:1612.05678.
  #
  # Args: 
  #   X: n x p data matrix of samples x variables (nodes); columns will be scaled 
  #      internally to have zero mean, unit variance.
  #   Gobs: p x p matrix of edge labels; each entry is in {NA, y1, y2} where 
  #         y1<y2 are real numbers and NA = unknown, y1 = non-causal, y2 = causal.
  #   nbins: (featurization) number of bins to use for each axis of 2d-histogram.
  #          Default is 30 bins, resulting in 30^2=900 features.
  #   binsRange: (featurization) vector giving axis limits for 2d-histogram. 
  #              Default is c(-3,3). Note that data has been internally scaled
  #              to have unit variance.
  #   outOfRange: (featurization) action to take for data points not within binsRange.
  #               "ignore" (default) excludes the data points from the histogram;
  #               "include" includes the data points in the first or last bin.
  #   q: (featurization) number of principal components to be used in PCA applied
  #      to the nbins^2 features. Default is 100.
  #   kernelType: string giving type of (ambient) kernel to use within Laplacian
  #               Regularized Least Squares (see Belkin et al., 2006); also serves
  #               as similarity matrix (matrix of weights) for the Laplacian. 
  #               Default is "rbfdot".
  #   neighboursCount: (kernel width parameter calculation) average distance is 
  #                    calculated between a randomly selected label and its 
  #                    neighboursCount nearest neighbour labels. This is repeated
  #                    batchSize times and results averaged. Default is 50 neighbours.
  #   batchSize: (kernel width parameter calculation) see neighboursCount.
  #              Default is 100.
  #   lambda1: vector of intrinsic regularization parameter(s). Default is 0.001.
  #            Paired elementwise with lambda2 (i.e. lambda1[i] is paired with
  #            lambda2[i]); length must therefore be the same as lambda2.
  #            A value of 0 means the Laplacian is not used and standard 
  #            Regularized Least Squares is performed.
  #   lambda2: vector of ambient regularization parameter(s). Default is 0.001.
  #            Paired elementwise with lambda1 (i.e. lambda1[i] is paired with
  #            lambda2[i]); length must therefore be the same as lambda1.
  #
  # Returns:
  #     Ghat: list of same length as lambda1 (and lambda2). Element i gives the
  #           p x p matrix of predicted edge label scores from applying MRCL with
  #           tuning parameters lambda1[i] and lambda2[i].
  #     feat: a p^2 x q matrix giving the featurization for each of the p^2 labels.
  #
  # Requires:
  #   install.packages('kernlab')
  #
  # References:
  # Hill SM, Oates C, Blythe DA, Mukherjee S (2019) Causal Learning via Manifold 
  #   Regularization. arXiv preprint arXiv:1612.05678.
  # Belkin M, Niyogi P, Sindhwani V (2006) Manifold regularization: A geometric
  #   framework for learning from labeled and unlabeled examples. JMLR 7:2399-2434.
  
  library(kernlab)
  
  stopifnot(ncol(X) == ncol(Gobs), length(lambda1) == length(lambda2))
  
  # convert observed label matrix Gobs to a vector, row-wise
  Y <- c(t(Gobs))
  
  # compute 2d-histogram counts featurization for each pair of nodes
  feat <- computeFeaturization(X, Gobs, nbins = nbins, binsRange = binsRange,
                                outOfRange = outOfRange, q = q)
  
  # compute width parameter for the kernel
  width <- computeKernelWidth(feat$postPCA, neighboursCount = neighboursCount,
                               batchSize = batchSize)
  
  # compute kernel matrix; also serves as similarity matrix (matrix of weights)
  # for Laplacian
  W <- computeKernel(feat$postPCA, type = kernelType, sigma = 1 / (width^2))
  
  # compute normalized Laplacian matrix
  d <- rowSums(W)
  dSqrtInv <- 1 / sqrt(d)
  D <- diag(d) # degree matrix
  L <- D - W # Laplacian matrix
  Lnorm <- t(dSqrtInv*t(dSqrtInv*L)) # equivalent to, but faster than, 
  #                                      (diag(dSqrtInv) %*% L) %*% diag(dSqrtInv)

  idxU <- is.na(Y) # unobserved label index
  u <- sum(idxU) # unobserved count
  l <- length(idxU) - u # observed count

  Y[idxU] <- 0 # set unobserved entries in Y to 0
  
  I <- diag(length(idxU)) # identity matrix

  # precompute some matrix multiplications
  JW <- W
  JW[idxU,] <- 0 # equivalent to, but faster than, J%*%W where J <- diag(!idxU)
  
  if(all(lambda1==0)) {
    LW <- 0 # no Laplacian regularization
  } else {
    LW <- Lnorm%*%W
  }
  
  # compute Ghat for each tuning parameter combination
  Ghats <- mapply(function(lambda1,lambda2) {
    
    # Equation (8) in Belkin et al. (2006)
    toBeInv <- JW + lambda2 * l * I + ((lambda1 * l) / (u + l)^2) * LW
    alphaStar <- solve(toBeInv, Y)
    
    # classification function scores
    fhat <- W %*% alphaStar
    Ghat <- matrix(fhat, nrow(Gobs), ncol(Gobs), byrow = TRUE)
    
    return(list(lambda1=lambda1,lambda2=lambda2,Ghat=Ghat))
    
  },lambda1,lambda2,SIMPLIFY=FALSE)
  
  return( list(Ghat=Ghats,feat=feat$postPCA) )

}


computeFeaturization <- function(X, Gobs, nbins = 30, binsRange = NULL,
                                 outOfRange = NULL, q = 100) {
  # Compute 2d-histogram featurization for each pair of nodes (variables)
  #
  # Args: 
  #   X: n x p data matrix of samples x variables (nodes); columns will be scaled 
  #      internally to have zero mean, unit variance.
  #   Gobs: p x p matrix of edge labels; each entry is in {NA, y1, y2} where 
  #         y1<y2 are real numbers and NA = unknown, y1 = non-causal, y2 = causal.
  #   nbins: number of bins to use for each axis of 2d-histogram.
  #          Default is 30 bins, resulting in 30^2=900 features.
  #   binsRange: vector giving axis limits for 2d-histogram. If NULL, the default,
  #              range(X) is used for all histograms.
  #   outOfRange: action to take for data points not within binsRange. 
  #               Not applicable if binsRange is NULL, otherwise must be set to 
  #               one of the following:
  #               "ignore" excludes the data points from the histogram;
  #               "include" includes the data points in the first or last bin.
  #   q: number of principal components to be used in PCA applied
  #      to the nbins^2 features. Default is 100.
  #
  # Returns:
  #     prePCA: a p^2 x nbins^2 matrix giving the featurization for each of the
  #             p^2 labels, prior to PCA.
  #     postPCA: a p^2 x q matrix giving the featurization for each of the
  #             p^2 labels, after applying PCA.

  if(!is.null(binsRange) && is.null(outOfRange)) {
    stop("Argument outOfRange needs to be specified: options are 'ignore' and 'include'")
  }
  
  stopifnot(!anyNA(X)) # check for missing values
  
  # standardize each column of X to have zero mean, unit variance
  X <- scale(X)
  
  # index mapping for labels
  getIndexMapping <- function(n, p) {
    K <- cbind("ij" = 1:(n*p), "i"  = rep(1:n, each = p), "j"  = rep(1:p, times = n))
    return(K)
  }
  p <- ncol(Gobs)
  K <- getIndexMapping(p, p)
  
  if(is.null(binsRange)) {
    binsRange <- range(X)
  } else {
    if(outOfRange=="include") {
      X[X<binsRange[1]] <- binsRange[1]
      X[X>binsRange[2]] <- binsRange[2]
    } else if(outOfRange!="ignore") {
      stop("Argument outOfRange should be set to 'ignore' or 'include'")
    }
  }
  
  computeFeatSingleLabel <- function(ij, K, X, nbins, binsRange, outOfRange) {
    # Computes featurization for a single label (pair of nodes)
    
    i <- K[ij, "i"]
    j <- K[ij, "j"]
    
    Xij <- X[, c(i, j)] # data for pair of nodes indexed by ij
    
    H <- hist2d(Xij[, 1], Xij[, 2], nbins = nbins, binsRange = binsRange,
                prob = TRUE)$z
    
    return(c(H)) # return histogram counts as a vector
  }
  
  # compute featurization for each label (pair of nodes)
  featPrePCA <- sapply(K[,"ij"], function(ij) 
    computeFeatSingleLabel(ij, K, X, nbins, rep(binsRange,2), outOfRange))
  featPrePCA <- t(featPrePCA)
  
  stopifnot(nrow(featPrePCA) == nrow(K))
  
  # reduce dimensionality (number of features) using PCA to q features
  pcaRes <- prcomp(featPrePCA, center = TRUE, scale. = FALSE, rank. = q)
  
  featPostPCA <- pcaRes$x
  
  return(list(postPCA = featPostPCA, prePCA = featPrePCA))
  
}

computeKernelWidth = function(feat, neighboursCount = 50, batchSize = 100) {
  # Compute the kernel width using an empirical nearest-neighbour approach
  #
  # Args: 
  #   feat: featurization matrix (labels x features).
  #   neighboursCount: average distance is calculated between a randomly selected
  #                    label and its neighboursCount nearest neighbour labels. 
  #                    This is repeated batchSize times and results averaged. 
  #                    Default is 50 neighbours.
  #   batchSize:  see neighboursCount. Default is 100.
  #
  # Returns:
  #     width: kernel width parameter
  
  perm <- sample.int(nrow(feat))
  num_keep <- min(nrow(feat), batchSize)
  perm <- perm[1:num_keep]
  
  distance <- rep(NA, length(perm))
  
  for (i in 1:length(perm)) {
    tmp <- sweep(feat, 2, feat[perm[i], ], "-")^2
    d <- sqrt(rowSums(tmp))
    stopifnot(nrow(d) == nrow(feat))
    nn <- order(d)[1:neighboursCount]
    distance[i] <- mean(d[nn])
  }
  
  width <- mean(distance)
  return(width)
}


computeKernel <- function(feat, type, ...) {
  # Compute a kernel matrix from the featurization
  #
  # Args: 
  #   feat: featurization matrix (labels x features).
  #   type: string giving type of kernel. Should match one of the kernel 
  #         functions provided in the kernlab package.
  #   ...: arguments for the kernlab package kernel function specified by type.
  #
  # Returns:
  #     kernMat: kernel matrix. Number of rows/columns is nrow(feat).
  #
  # Requires:
  #   install.packages('kernlab')
  
  fcn <- get(type)
  kernel <- fcn(...)
  kernMat <- kernelMatrix(kernel, feat)
  
  return(kernMat)
}


hist2d <- function(x, y, nbins = 30, binsRange = c(range(x), range(y)),
                   prob = FALSE) {
  # Compute 2d-histogram for a pair of variables
  #
  # Args: 
  #   x: vector
  #   y: vector, same length as x
  #   nbins: number of bins to use for each axis of 2d-histogram.
  #          Default is 30 bins, resulting in 30^2=900 features.
  #   binsRange: vector of length 4 giving axis limits for 2d-histogram. First 2
  #              elements are limits for x-axis, last 2 elements are for y-axis.
  #              Default is c(range(x),range(y)).
  #   prob: logical; if TRUE, normalizes histogram counts so they sum to one.
  #
  # Returns:
  #   z: table of histogram counts, normalized if prob is TRUE
  #   xBins: x-axis bin boundaries
  #   yBins: y-axis bin boundaries
  
  if (missing(x)) stop('Missing X input vector.')
  if (missing(y)) stop('Missing Y input vector.')
  if (length(x) != length(y)) stop('X and Y have different lengths.')
  
  # determine x, y limits of the 2D histogram
  xrange <- binsRange[1:2]
  yrange <- binsRange[3:4]
  
  # determine bin width of the 2D histogram
  dx <- (xrange[2] - xrange[1]) / nbins
  dy <- (yrange[2] - yrange[1]) / nbins
  
  # determine the bin boundries
  xBins <- seq(xrange[1], xrange[2], by = dx)
  yBins <- seq(yrange[1], yrange[2], by = dy)
  
  # for each data point (x,y) find which bin it falls into
  x[x==xBins[1]] <- xBins[1] + dx/2 # deal with lower boundary case, so counted in first bin
  y[y==yBins[1]] <- yBins[1] + dy/2
  xBinned <- ceiling((x-xBins[1]) / dx)
  yBinned <- ceiling((y-yBins[1]) / dy)
  
  # compute table, showing counts in each 2d bin
  # Note: we need the explicit factor conversion here, otherwise
  # rows/columns with zero counts will be dropped from the table.
  freq2D <- table(factor(xBinned, levels = 1:nbins), 
                   factor(yBinned, levels = 1:nbins))
  
  # if prob is TRUE then normalise so that the sum of the counts is 1.
  if (prob == TRUE)  {
    norm <- sum(freq2D)
    freq2D <- freq2D / norm
  }
  
  # output the boundaries x and y ordinates, and the 2D histogram (z
  return(list(z = freq2D, xBins = xBins, yBins = yBins))
  
}