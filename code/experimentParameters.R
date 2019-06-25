getParameters <- function(experiment) {
  
  if(experiment == "TCPA") {
    # TCPA ---------------------------------------------------------------------
    
    rhoSet <- c(0.2, 0.4, 0.6, 0.8) # proportion of edge labels to treat as observed labels for training
    nSet <- c(200, 500, "max") # training data sample sizes
    
    dataPar <- list(
      # dataset parameters
      experiment = experiment,
      subtype = "BRCA"
      )
    
    labelSplitPar <- list(
      # label train/test split parameters
      experiment = experiment,
      propnEdgesInclDiag = FALSE, # exclude self-loops for stratified sampling of edges and non-edges?
      M = 25, # number of train/test splits to generate
      method = "random", # sample at random or sample whole rows at a time
      propnEdgesRowLB = NULL # if method="rows", stratifies sampling by rows with and without at least propnEdgesRowLB edges
    )
    
    dataSampPar <- list(
      # data subsampling parameters
      experiment = experiment,
      M = labelSplitPar$M, # number of subsampled datasets to generate
      scaled = TRUE # standardise each gene after subsampling
      )
    
    methPar <- list(
      # method parameters
      experiment = experiment,
      doMRCL = TRUE,
      doCor = TRUE,
      doLasso = TRUE,
      doIDA = TRUE,
      doRFCI = TRUE,
      doKNN = TRUE,
      doGIES = FALSE,
      M_MRCL = 1, # number of iterations of MRCL
      lambda1_MRCL = 0.001, # MRCL intrinsic (Laplacian) regularisation tuning parameter
      lambda2_MRCL = 0.001, # MRCL ambient (RKHS) regularisation tuning parameter
      M_Lasso = 1, # number of iterations of Lasso
      sigLevelIDA = 0.01, # IDA significance level for PC algorithm CI tests
      sigLevelRFCI = 0.01, # RFCI significance level for PC algorithm CI tests
      nPrCmpKNN = 2, # number of principal components to keep for KNN
      Ktune = data.frame(k = seq(3, 49, by = 2)), # tuning parameters for KNN (# neighbours)
      M = max(labelSplitPar$M, dataSampPar$M) # number of iterations
    )
    
    
    
  } else if(experiment %in% c("yeast_random", "yeast_row-wise", "yeast_row-wise_GIES")) {
    # yeast --------------------------------------------------------------------
    
    rhoSet <- c(0.2, 0.4, 0.6, 0.8) # proportion of edge labels to treat as observed labels for training
    if(experiment == "yeast_row-wise_GIES") {
      nSet <- NA # for this experiment, sample size is a function of rho: n=nrow(yeastObsData)/2 + floor(rho*dataPar$nGeneSubset)
    } else {
      nSet <- c(200, 500, 1000) # training data sample sizes
    }
    
    dataPar <- list(
      experiment = experiment,
      tau = 5, # threshold for calling edges
      nGeneSubset = 50, # number of genes in network
      M = 25, # number of datasets to generate, each with a different gene subset
      propnEdgesGlobalLB = 0.025, # lower bound on proportion of labels called as edges across whole network
      propnEdgesInclDiag = FALSE # exclude self-loops in the calculation of proportion of edges?
    )
    if(experiment == "yeast_random") {
      dataPar$propnEdgesRowLB <- NULL # lower bound on proportion of labels called as edges for a single row of adjacency matrix
      dataPar$propnEdgesRowPropn <- NULL # proportion of rows of adjacency matrix that need to satisfy the propnEdgesRowLB constraint
      dataPar$trainInterventions <- "non-nodes" # set to "nodes"/"non-nodes" to use interventions that do/do not target network nodes as training samples
    } else if(experiment == "yeast_row-wise") {
      dataPar$propnEdgesRowLB <- 1/50 
      dataPar$propnEdgesRowPropn <- 0.5
      dataPar$trainInterventions <- "non-nodes" 
    } else if(experiment == "yeast_row-wise_GIES") {
      dataPar$propnEdgesRowLB <- 1/50 
      dataPar$propnEdgesRowPropn <- 0.5
      dataPar$trainInterventions <- "nodes" 
    }
    
    labelSplitPar <- list(
      experiment = experiment,
      propnEdgesInclDiag = dataPar$propnEdgesInclDiag, # exclude self-loops for stratified sampling of edges and non-edges?
      M = NULL, # if NULL, a single train/test split of labels is generated for each of the dataPar$M datasets
      propnEdgesRowLB = dataPar$propnEdgesRowLB # if method="rows", stratifies sampling by rows with and without at least propnEdgesRowLB edges
    )
    if(experiment == "yeast_random") {
      labelSplitPar$method = "random" # sample at random
    } else if(experiment %in% c("yeast_row-wise", "yeast_row-wise_GIES")) {
      labelSplitPar$method = "rows" # sample whole rows at a time
    }
    
    dataSampPar <- list(
      experiment = experiment,
      M = NULL, # if NULL, a single subsampling is generated for each of the dataPar$M datasets
      scaled = TRUE, # standardise each gene after subsampling
      fixed = "obs", # if set to "obs", always include all observational training data
      exclude = NULL, # if set to "obs", exclude observational training data from training data
      trainInterventions = dataPar$trainInterventions # set to "nodes"/"non-nodes" to use interventions that do/do not target network nodes as training samples
    )
    if(dataSampPar$trainInterventions == "nodes" && labelSplitPar$method == "random") {
      stop("invalid parameters: interventions on nodes in training data must be
            used with row-wise subsampling of labels")
    }
    
    methPar <- list(
      experiment = experiment,
      doMRCL = TRUE,
      doCor = TRUE,
      doLasso = TRUE,
      doIDA = TRUE,
      doRFCI = TRUE,
      doKNN = TRUE,
      doGIES = FALSE,
      M_MRCL = 1, # number of iterations of MRCL
      lambda1_MRCL = 0.001, # MRCL intrinsic (Laplacian) regularisation tuning parameter
      lambda2_MRCL = 0.001, # MRCL ambient (RKHS) regularisation tuning parameter
      M_Lasso = 1, # number of iterations of Lasso
      sigLevelIDA = 0.01, # IDA significance level for PC algorithm CI tests
      sigLevelRFCI = 0.01, # RFCI significance level for PC algorithm CI tests
      nPrCmpKNN = 2, # number of principal components to keep for KNN
      Ktune = data.frame(k = seq(3, 49, by = 2)), # tuning parameters for KNN (# neighbours)
      M = max(dataPar$M, labelSplitPar$M, dataSampPar$M) # number of iterations
    )
    
    if(experiment == "yeast_row-wise_GIES") {
      methPar$maxDegreeGIES <- 0 # tuning parameter for GIES. Set to 0 for default of no maxDegree restriction.
      methPar$doGIES <- TRUE
    }
    
    
  } else if(experiment=="cellLine") {
    # cellLine -----------------------------------------------------------------
    
    rhoSet <- NA
    nSet <- NA
    
    dataPar <- list(
      experiment = experiment,
      propStimsThresh = 0.5,
      scaled = TRUE
    )
    
    labelSplitPar <- list()
    dataSampPar <- list()
    
    methPar <- list(
      experiment = experiment,
      doMRCL = TRUE,
      doCor = TRUE,
      doLasso = TRUE,
      doIDA = TRUE,
      doRFCI = TRUE,
      doKNN = TRUE,
      doGIES = FALSE,
      M_MRCL = 10, # number of iterations of MRCL
      lambda1_MRCL = 0.001, # MRCL intrinsic (Laplacian) regularisation tuning parameter
      lambda2_MRCL = 0.001, # MRCL ambient (RKHS) regularisation tuning parameter
      M_Lasso = 10, # number of iterations of Lasso
      sigLevelIDA = 0.01, # IDA significance level for PC algorithm CI tests
      sigLevelRFCI = 0.01, # RFCI significance level for PC algorithm CI tests
      nPrCmpKNN = 2, # number of principal components to keep for KNN
      Ktune = data.frame(k = seq(3, 49, by = 2)), # tuning parameters for KNN (# neighbours)
      M = 4 # 4 cell lines
    )
    
  }
  
  return(list(dataPar = dataPar, labelSplitPar = labelSplitPar,
              dataSampPar = dataSampPar, methPar = methPar, rhoSet = rhoSet, nSet = nSet))
  
}