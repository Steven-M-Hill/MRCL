library(pROC)
library(glmnet)
library(caret)
library(pcalg)
library(dplyr)
library(R.matlab)

experiment <- "cellLine"

seed <- 2865018
set.seed(seed)

source("code/experimentFunctions.R")
source("code/experimentParameters.R")
source("code/mrcl.R")

###### parameters ------------------------------------------------------

out <- getParameters(experiment)
list2env(out, envir = globalenv())
rm(out)

###### load RPPA cell line data and CDMs -----------------------
if(!file.exists("data/cellLineData.RData")) {
  importCellLineData()
}
load("data/cellLineData.RData")

####################

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
if(file.exists(fname)) {
  load(fname)
} else {
  out <- getCellLineTrainingDataLabelData(rppaData, CDMs, 
                                          propStimsThresh = dataPar$propStimsThresh,
                                          scaled = dataPar$scaled)
  list2env(out, envir = globalenv())
  rm(out)
  
  params <- list(data = dataPar)
  save(dataTrain, dataLabels, params, file = fname)
}

# training labels and test labels
fname <- file.path("output", experiment, "data", "labelSplit.RData")
if(file.exists(fname)) {
  load(fname)
} else {
  Gobs <- lapply(dataLabels, function(x) {
    x[testInhibitorTarget, ] <- NA
    return(x)
  })
  Gtest <- lapply(dataLabels, function(x) {
    tmp <- matrix(nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))
    tmp[testInhibitorTarget, ] <- x[testInhibitorTarget, ]
    return(tmp)
  })
  
  params <- list(data = dataPar)
  save(Gobs, Gtest, testInhibitor, testInhibitorTarget, params, file = fname)
} 

params$meth <- methPar

fname <- file.path("output", experiment, "results", "results.RData")
res <- initRes(params)$res
res <- lapply(res, lapply, function(x) {
  if(is.vector(x) && length(x) == length(cellLines)) {
    names(x) <- cellLines
  } else if(is.matrix(x) & ncol(x) == length(cellLines)) {
    colnames(x) <- cellLines
  }
  return(x)
})

res <- runMethods(dataTrain, methPar, res, Gobs, Gtest)

save(res, params, testInhibitor, testInhibitorTarget, file=fname)


#####


fname <- file.path("output", experiment, "results", "collatedResults.RData") 
out <- collateAUCs(experiment)
list2env(out, envir=globalenv())
aucs <- bind_rows(aucs, .id = "Method")

out <- collatePointEstimateScores(experiment, "pointEval")
list2env(out, envir = globalenv())
save(scores, aucs, avgROCs, params, file = fname)
