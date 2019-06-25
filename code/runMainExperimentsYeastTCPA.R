library(dplyr)
library(doParallel)
library(doRNG)

experiment <- "yeast_row-wise" # one of "TCPA", "yeast_random", "yeast_row-wise", "yeast_row-wise_GIES"
generateDatasets <- TRUE # set to FALSE if want to use training data / label data generated in a previous run

# set seed
if(experiment == "TCPA") {
  seed <- 49254205
} else if(experiment == "yeast_random") {
  seed <- 4634629
} else if(experiment == "yeast_row-wise") {
  seed <- 172301
} else if(experiment == "yeast_row-wise_GIES") {
  seed <- 3834740
}
set.seed(seed)

cl <- makeCluster(3) 
registerDoParallel(cl)

source("code/experimentFunctions.R")
source("code/experimentParameters.R")
source("code/mrcl.R")

###### parameters --------------------------------------------------------------

out <- getParameters(experiment)
list2env(out,envir=globalenv())
rm(out)

if(experiment == "TCPA") {
  
  # load TCPA data
  load("data/TCGA-PANCAN19-L4-BRCA-35phosphoSubset.RData")
  
} else {
  
  # load yeast data
  fname <- paste0("data/yeastData.RData")
  if(!file.exists(fname)) {
    importYeastData()
  }
  load(fname)
}

###### generate M sets of training data and label data -------------------------
# each set is for a different subset of nGeneSet genes and has its own training/test split for the observational data

if(!dir.exists("output")) {
  dir.create("output")
}

if(!dir.exists(file.path("output", experiment))) {
  dir.create(file.path("output", experiment))
  dir.create(file.path("output", experiment, "data"))
  dir.create(file.path("output", experiment, "results"))
  dir.create(file.path("output", experiment, "plots"))
}

fname <- file.path("output", experiment, "data", "data.RData")
if(!generateDatasets) {
  load(fname)
} else {
  
  if(experiment == "TCPA") {
    out <- getTCPAtrainingDataLabelData(tcpaData)  
  } else {
    out <- getYeastTrainingDataLabelData(yeastObsData, yeastIntData, yeastIntIndices, dataPar$nGeneSubset, dataPar$tau, dataPar$M, dataPar$propnEdgesGlobalLB, dataPar$propnEdgesRowLB, dataPar$propnEdgesRowPropn, dataPar$propnEdgesInclDiag, dataPar$trainInterventions)  
  }
  
  list2env(out, envir = globalenv())
  rm(out)
  
  params <- list(data = dataPar)
  postDataSeed <- .Random.seed
  save(dataTrain, dataLabels, params, postDataSeed, file = fname)
}

generateSubsampDatasets <- generateDatasets

for(rho in rhoSet) {

  # training labels and test labels
  fname <- file.path("output", experiment, "data", paste0("labelSplit_rho=", rho, ".RData"))
  if(!generateDatasets) {
    load(fname)
  } else {
    out <- generateTrainTestLabelData(dataLabels, rho, M = labelSplitPar$M, 
                                      propnEdgesInclDiag = labelSplitPar$propnEdgesInclDiag,
                                      method = labelSplitPar$method, 
                                      propnEdgesRowLB = labelSplitPar$propnEdgesRowLB)
    list2env(out,envir = globalenv())
    
    params <- list(data = dataPar, labelSplit = labelSplitPar)
    params$labelSplit$rho <- rho
    save(Gobs, Gtest, params, file = fname)
  } 
  
  ###### run methods -------------------------------------------------------------
  
  # set n if training data interventions are on network nodes
  if(experiment != "TCPA" && dataSampPar$trainInterventions == "nodes") {
    nSet <- floor(nrow(yeastObsData) / 2) + floor(rho * dataPar$nGeneSubset)
  } 
  out <- foreach(n = nSet, .packages=c("pROC", "glmnet", "pcalg", "kernlab", "caret"), 
                 .final = function(x) setNames(x, paste0("n=", nSet))) %dorng% { 
    
    # generate subsampled training data
    fname <- file.path("output", experiment, "data", paste0("dataSubsample_n=", n, ".RData"))
    if(!generateSubsampDatasets) {
      load(fname)
    } else {
      fixed <- NULL
      exclude <- NULL
      
      if(!is.null(dataSampPar$fixed) && dataSampPar$fixed == "obs") {
        # include observational samples in training data
        fixed <- lapply(dataTrain, function(x) grep("obs", rownames(x)))
      } else if(!is.null(dataSampPar$exclude) && dataSampPar$exclude == "obs") {
        # exclude observational samples from training data
        exclude <- lapply(dataTrain, function(x) grep("obs", rownames(x)))
      }
      # if training data has interventions on network nodes, exclude those samples that are used as test label data in Gtest
      if(experiment != "TCPA" && dataPar$trainInterventions == "nodes") {
        excludeTestInterventions <- mapply(function(x, y) 
          which(rownames(x) %in% rownames(y)[rowSums(!is.na(y)) > 0]), 
          dataTrain, Gtest, SIMPLIFY = FALSE)
        if(!is.null(exclude)) {
          exclude <- mapply(function(x, y) c(x, y), 
                            exclude, excludeTestInterventions, SIMPLIFY = FALSE)
        } else {
          exclude <- excludeTestInterventions
        }
      }
      
      if(n == "max") {
        if(experiment == "TCPA") {
          nActual <- nrow(dataTrain)
        } else {
          nActual <- nrow(dataTrain[[1]]) - sapply(exclude,length)
        }
      } else {
        nActual <- as.numeric(n)
      }
      dataTrainSubsample <- generateSubsampledData(dataTrain, nActual, M = dataSampPar$M,
                                                   geneIdxInt = geneIdxInt, fixed = fixed,
                                                   exclude = exclude, scaled = dataSampPar$scaled)
      
      params <- list(data = dataPar, dataSamp = dataSampPar)
      params$dataSamp$n <- nActual
      save(dataTrainSubsample, params, file = fname)
    } 
    
    # store params
    params$labelSplit <- labelSplitPar
    params$meth <- methPar
    params$labelSplit$rho <- rho
    
    # initialise list to store results, loading previous results if applicable
    fname <- file.path("output", experiment, "results", paste0("n=", n, "_rho=", rho, ".RData"))
    res <- initRes(params)$res
    
    res <- runMethods(dataTrainSubsample, methPar, res, Gobs, Gtest, geneIdxNode)
    
    save(res, params, file = fname)
  }
  if(experiment == "TCPA" || dataSampPar$trainInterventions != "nodes") {
    generateSubsampDatasets <- FALSE # only need to be generated for first iteration of rho
  }  
}
stopCluster(cl)
  
####### collate auc scores across values of n and rho, calculate means and std devs ------------------------------------------------

if(experiment != "TCPA" && dataSampPar$trainInterventions == "nodes") {
  nSet <- floor(nrow(yeastObsData) / 2) + floor(rhoSet * dataPar$nGeneSubset)
  single_n_per_rho <- TRUE
} else {
  single_n_per_rho <- FALSE
}

fname <- file.path("output", experiment, "results", "collatedResults.RData") 
out <- collateAUCs(experiment, rhoSet, nSet, single_n_per_rho = single_n_per_rho)
list2env(out, envir = globalenv())
aucsMeanSd <- aucMeanSd(aucs)

out <- collatePointEstimateScores(experiment, "pointEval", rhoSet, nSet, single_n_per_rho = single_n_per_rho)
list2env(out, envir = globalenv())
scoresMeanSd <- pointEstimateScoresMeanSd(scores)
save(scores, scoresMeanSd, aucs, aucsMeanSd, avgROCs, params, file=fname)
