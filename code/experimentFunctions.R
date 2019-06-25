# functions for the application of the MRCL method

importCellLineData <- function() {
  
  # import RPPA data from https://github.com/Steven-M-Hill/causal-signaling-networks-CellSystems2016 ---------------------------------------------------------
  githubURL <- "http://raw.githubusercontent.com/Steven-M-Hill/causal-signaling-networks-CellSystems2016/master/"
  download.file(paste0(githubURL,"data/DataS1.zip"),"data/DataS1.zip") # download data in CSV format
  cellLines <- c("UACC812","MCF7","BT20","BT549")
  unzip("data/DataS1.zip",exdir = "data",files = paste0("DataS1/core/",cellLines,".csv"))
  unlink("data/DataS1.zip")
  
  rppaData <- vector("list", length(cellLines))
  names(rppaData) <- cellLines
  for(i in 1:length(cellLines)) {
    rppaData[[i]] <- read.csv(paste0("data/DataS1/core/", cellLines[i], ".csv"),stringsAsFactors=FALSE, check.names = FALSE)
    colnames(rppaData[[i]])[1:6] <- rppaData[[i]][2,1:6]
    rppaData[[i]] <- rppaData[[i]][-c(1:2), -c(1,5)]
    rownames(rppaData[[i]]) <- NULL
    rppaData[[i]][1:3] <- lapply(rppaData[[i]][1:3], factor)
    rppaData[[i]][4] <- as.numeric(rppaData[[i]][[4]])
  }
  
  # remove any samples labelled as excluded and remove this metadata column
  for(i in 1:length(cellLines)) {
    rppaData[[i]] <- subset(rppaData[[i]],`Sample excluded from analyses?`==0)
    rppaData[[i]] <- rppaData[[i]][,-4]
    rppaData[[i]][1:3] <- lapply(rppaData[[i]][1:3],factor)
  }
  
  # take logs and average replicates
  rppaData <- lapply(rppaData,function(x) cbind(x[,1:3],log2(x[,-c(1:3)])))
  rppaData <- lapply(rppaData,function(x) {
    aggregate(x[,-c(1:3)],x[,1:3],mean,na.rm=TRUE)
  })
  
  unlink("data/DataS1", recursive = TRUE)
  
  # # import CDM data from https://github.com/Steven-M-Hill/causal-signaling-networks-CellSystems2016 --------------------------------------------------------- ---------------------------------------------
  
  matData <- readMat(url(paste0(githubURL,"results/CDMs/CDMs.mat")))
  
  CDMs <- matData$CDMs
  CDMproteinNames <- unname(unlist(matData$proteinNames))
  CDMinhibitorRegimes <- unname(unlist(matData$inhibitor))
  CDMcontexts <- unname(unlist(matData$contextLabels))
  dimnames(CDMs) <- list(CDMproteinNames,CDMcontexts,CDMinhibitorRegimes)
  
  proteins <- CDMproteinNames
  
  save(rppaData, CDMs, cellLines, proteins, file="data/cellLineData.RData")
  
}


importYeastData <- function() {
  # import yeast data from the Deleteome Collection (http://deleteome.holstegelab.nl/) ---------------------------------------------------------
  # RData file Kemmermen.RData on Downloads/Causal inference subpage of above website
  fileURL <- "http://deleteome.holstegelab.nl/data/downloads/causal_inference/Kemmeren.RData"
  load(url(fileURL))
  
  # extract interventional data
  yeastIntData <- data$int
  yeastIntIndices <- data$intpos
  
  # extract observational data (use only the 161 samples taken from file limma_MARG_wt_pool.txt on above website)
  yeastObsData <- data$obs[1:161, ]
  # check for duplicated samples and remove them
  tmp <- cor(t(yeastObsData))
  idx <- which(tmp*upper.tri(tmp)>0.999, arr.ind = TRUE)
  yeastObsData <- yeastObsData[-(idx[,2]), ]
  
  save(yeastObsData, yeastIntData, yeastIntIndices, file="data/yeastData.RData")
}


getTCPAtrainingDataLabelData <- function(data) {
  
  dataTrain <- data[,-(1:4)] # remove metadata columns
  dataLabels <- getPriorGraph()
  
  return(list(dataTrain = dataTrain, dataLabels = dataLabels, geneIdxInt = NULL))
  
}


getPriorGraph <- function() {
  # prior graph as in Cell Systems paper
  
  priorGraph <- matrix(0, nrow = 35, ncol = 35)
  
  EBP1 <- 1
  ACC <- 2
  AKTpS<-3; AKTpT<-4; AKT <- c(AKTpS, AKTpT)
  AMPK <- 5
  BAD <- 6
  cMET <- 7
  cRAF <- 8
  CHK1 <- 9
  CHK2 <- 10
  EGFRpY1068 <- 11; EGFRpY1173 <- 12; EGFR <- c(EGFRpY1068, EGFRpY1173)
  ERa <- 13
  GSK3 <- 14
  HER2 <- 15
  JNK <- 16
  MAPK <- 17
  MEK <- 18
  mTOR <- 19
  NFkB <- 20
  p27pT157 <- 21; p27pT198 <- 22; p27 <- c(p27pT157, p27pT198)
  p38 <- 23
  p70 <- 24
  p90 <- 25
  PDK1 <- 26
  PKCa <- 27
  PRAS40 <- 28
  RB <- 29
  S6 <- 30
  SRCpY416 <- 31; SRCpY527 <- 32; SRC <- c(SRCpY416, SRCpY527)
  STAT3 <- 33
  YAP <- 34
  YB1 <- 35
  
  priorGraph[c(EGFR, HER2, cMET), PKCa] <- 1
  priorGraph[c(EGFR, HER2, cMET), PDK1] <- 1
  priorGraph[c(EGFR, HER2, cMET), cRAF] <- 1
  priorGraph[c(EGFR, HER2, cMET), SRC] <- 1
  priorGraph[c(EGFR, HER2, cMET), p38] <- 1
  priorGraph[c(EGFR, HER2, cMET), JNK] <- 1
  priorGraph[AKT, PRAS40] <- 1
  priorGraph[AKT, mTOR] <- 1
  priorGraph[AKT, GSK3] <- 1
  priorGraph[AKT, BAD] <- 1
  priorGraph[AKT, YAP] <- 1
  priorGraph[AMPK, ACC] <- 1
  priorGraph[AMPK, mTOR] <- 1
  priorGraph[cRAF, MEK] <- 1
  priorGraph[GSK3, mTOR] <- 1
  priorGraph[MAPK, p90] <- 1
  priorGraph[MAPK, cRAF] <- 1
  priorGraph[MAPK, MEK] <- 1
  priorGraph[MAPK, mTOR] <- 1
  priorGraph[MEK, MAPK] <- 1
  priorGraph[mTOR, p70] <- 1
  priorGraph[mTOR, EBP1] <- 1
  priorGraph[mTOR, AKT] <- 1
  priorGraph[p70, S6] <- 1
  priorGraph[p70, PDK1] <- 1
  priorGraph[p70, mTOR] <- 1
  priorGraph[p90, S6] <- 1
  priorGraph[p90, mTOR] <- 1
  priorGraph[p90, GSK3] <- 1
  priorGraph[p90, BAD] <- 1
  priorGraph[PDK1, AKT] <- 1
  priorGraph[PDK1, PKCa] <- 1
  priorGraph[PDK1, p70] <- 1
  priorGraph[PRAS40, mTOR] <- 1
  priorGraph[SRC, STAT3] <- 1
  
  return(priorGraph)
}


getYeastTrainingDataLabelData <- function(obsData, intData, intIndices, nGeneSubset,
                                          tau, M, propnEdgesGlobalLB = NULL, 
                                          propnEdgesRowLB = NULL, 
                                          propnEdgesRowPropn = NULL, 
                                          propnEdgesInclDiag = TRUE,
                                          trainInterventions = "non-nodes") {
  
  dataTrain <- vector("list", M) # training data
  dataZscores <- vector("list", M) # causal effects
  dataLabels <- vector("list", M) # labels
  geneIdxNode <- vector("list", M) # gene index for network nodes
  geneIdxInt <- vector("list", M) # gene index for interventional targets in training data
  propnEdgesGlobal <- rep(-1, M)
  propnEdgesRow <- matrix(-1, nrow = nGeneSubset, ncol = M)
  if(is.null(propnEdgesGlobalLB)) {
    propnEdgesGlobalLB <- 0
  }
  if(is.null(propnEdgesRowLB)) {
    propnEdgesRowLB <- 0
    propnEdgesRowPropn <- 0
  }
  
  nObsTot <- nrow(obsData)
  nIntTot <- nrow(intData)
  
  for(m in 1:M) {
    
    # split observational samples into train and test, and 
    #     select gene subset C \subset T (set of all targeted genes) to take as network nodes.
    #     (expression levels for genes not in C are not used in this analysis)
    
    done <- propnEdgesGlobal[m] >= propnEdgesGlobalLB
    done <- done & mean(propnEdgesRow[,m] >= propnEdgesRowLB) >= propnEdgesRowPropn 
    
    cnt <- 0
    
    while(!done) {
      
      cnt <- cnt + 1
      
      trainIdx <- sample(nObsTot, floor(nObsTot/2)) # obs training data samples, the rest are test samples
      intIdxNode <- sample(nIntTot, nGeneSubset) # subset C of targeted genes to take as network nodes - interventional data samples index
      geneIdxNode_m <- intIndices[intIdxNode] # subset C of targeted genes to take as network nodes - genes index
      
      obsDataTrain <- obsData[trainIdx, geneIdxNode_m] # obs training data for gene subset C
      obsDataTest <- obsData[-trainIdx, geneIdxNode_m] # obs test data for gene subset C
      
      # interventional data for network nodes, used to obtain "ground truth" edge labels
      # matrix of size nGeneSubset x nGeneSubset. 
      # entry (i,j) is expression of gene j under intervention of gene i, where i and j are in the subset C
      intDataNode <- intData[intIdxNode, geneIdxNode_m] 
      
      # robust z-score for changes under intervention, relative to observational test data
      dataZscores[[m]] <- scale(intDataNode, center = apply(obsDataTest, 2, median),
                                scale = apply(obsDataTest, 2, IQR, type = 2))
      # robust z scores thersholded at tau to obtain binary "edge labels" matrix G, for genes in the subset C under intervention of genes in C
      dataLabels[[m]] <- (abs(dataZscores[[m]]) >= tau)*1 # the *1 converts logical to numeric
      
      if(propnEdgesInclDiag) {
        propnEdgesGlobal[m] <- sum(dataLabels[[m]]) / (nGeneSubset ^ 2)
        propnEdgesRow[, m] <- apply(dataLabels[[m]], 1, sum) / nGeneSubset
      } else {
        tmp <- dataLabels[[m]]
        diag(tmp) <- NA
        propnEdgesGlobal[m] <- sum(tmp, na.rm = TRUE) / (nGeneSubset ^ 2 - nGeneSubset)
        propnEdgesRow[, m] <- apply(tmp, 1, sum, na.rm = TRUE) / (nGeneSubset - 1)
      }
      
      if(trainInterventions == "non-nodes") {
        # full training data matrix X, consists of observational training data and 
        #     interventions for non-network nodes (interventions on network nodes are used as test data/labels)
       
        intDataTrain <- intData[-intIdxNode, geneIdxNode_m] # expression levels for genes in C, under interventions of genes not in C (i.e. not in test data)
        geneIdxInt_m <- intIndices[-intIdxNode]
      } else if(trainInterventions == "nodes") {
        # full training data matrix X, consists of observational training data and interventions for network nodes 
        #     (some of these samples will be removed later and used as test data/labels)
        intDataTrain <- intData[intIdxNode, geneIdxNode_m] # expression levels for genes in C, under interventions of genes not in C (i.e. not in test data)
        geneIdxInt_m <- geneIdxNode_m
      } else {
        stop("trainInterventions should be set to 'nodes' or 'non-nodes'")
      }
      dataTrain[[m]] <- scale(rbind(obsDataTrain, intDataTrain))
      geneIdxInt[[m]] <- c(rep(NA, length(trainIdx)), geneIdxInt_m)
      
      geneIdxNode[[m]] <- geneIdxNode_m
      
      done <- propnEdgesGlobal[m] >= propnEdgesGlobalLB
      done <- done & mean(propnEdgesRow[,m] >= propnEdgesRowLB) >= propnEdgesRowPropn 
      
    }
    print(cnt)
  }
  return(list(dataTrain = dataTrain, dataLabels = dataLabels, dataZscores = dataZscores,
              geneIdxNode = geneIdxNode, geneIdxInt = geneIdxInt, 
              propnEdgesGlobal = propnEdgesGlobal, propnEdgesRow = propnEdgesRow))
}


getCellLineTrainingDataLabelData <- function(data, CDMs, propStimsThresh, 
                                             scaled = FALSE) {
  
  cellLine <- names(data)
  testInhibitor <- "AZD8055"
  holdOut <- c("AZD8055", "BEZ235") # BEZ235 data are not used at all, AZD8055 data are used to obtain the test labels and so are excluded from the training data matrix
  
  ##### form training data matrix ---------------------
  
  dataTrain <- lapply(data, function(x) x[!(x$Inhibitor %in% holdOut), -(1:3)]) # hold out BEZ and AZD
  if(scaled) {
     dataTrain <- lapply(dataTrain, scale)
  }
  
  ##### labels, derived from the causal descendancy matrices (which in turn were
  #     derived from the interventional data) ----------

  dataLabels <- lapply(dataTrain, function(x) {
    x <- matrix(ncol = ncol(x), nrow = ncol(x), dimnames = list(colnames(x), colnames(x)))
  })
    
  # set missing data in CDMs to NA, these values are not valid due to an issue with the experiment
  CDMs[, grepl("BT549", colnames(CDMs)), 5] <- NA
  CDMs[, grepl("BT20,PBS", colnames(CDMs)) | grepl("BT20,NRG1", colnames(CDMs)), 5] <- NA
  
  # labels for the test inhibitor
  getDesc <- function(inhibitor) {
    descCellLineStims <- lapply(cellLine, function(x) CDMs[, grepl(x, colnames(CDMs)), inhibitor]) # CDM data for one cell line (all 8 stimuli) and one inhibitor
    descCellLine <- lapply(descCellLineStims, function(x) (rowMeans(x, na.rm = TRUE) >= propStimsThresh) * 1 )
    return(descCellLine)
  }
  desc <- getDesc(testInhibitor)
  
  inhibTarget <- "mTOR_pS2448"
  setLabels <- function(inhibTarget, desc) {
    mapply(function(x, y) {
      x[inhibTarget, ] <- matrix(rep(y, length(inhibTarget)), 
                                 nrow = length(inhibTarget), byrow = TRUE)
      return(x)
    }, dataLabels, desc, SIMPLIFY = FALSE) 
  }
  dataLabels <- setLabels(inhibTarget, desc)
  
  trainingInhibitor1 <- "GSK690693" # targets akt
  trainingInhibitor2 <- "GSK690693_GSK1120212" # targets akt and mek
  desc1 <- getDesc(trainingInhibitor1)
  desc2 <- getDesc(trainingInhibitor2)
  
  inhibTarget <- c("Akt_pS473", "Akt_pT308")
  descAkt <- mapply(function(x, y) {
    desc <- x & (is.na(y) | y) # agreement from both inhibitors (just use 1st inhibitor if 2nd is missing)
    desc[x & !y] <- NA
    return(desc)
  
    }, desc1, desc2, SIMPLIFY = FALSE)
  dataLabels <- setLabels(inhibTarget, descAkt)
  
  inhibTarget <- "MEK1_pS217_S221"
  descMek <- mapply(function(x, y) {
    desc <- !x & y # TRUE for 2nd inhibitor and FALSE for 1st inhibitor
    desc[is.na(y)] <- NA
    return(desc)
  }, desc1, desc2, SIMPLIFY = FALSE)
  dataLabels <- setLabels(inhibTarget, descMek)
  
  return(list(dataTrain = dataTrain, dataLabels = dataLabels, 
              testInhibitor = testInhibitor, testInhibitorTarget = "mTOR_pS2448"))
}


generateTrainTestLabelData <- function(labels, rho, M = NULL, propnEdgesInclDiag = TRUE,
                                       method = "random", propnEdgesRowLB = NULL) {
  if(is.list(labels)) {
    M <- length(labels)
    out <- lapply(labels, generateTrainTestLabelData, rho, M = 1, 
                  propnEdgesInclDiag, method, propnEdgesRowLB)
    Gobs <- sapply(out, "[[", "Gobs")
    Gtest <- sapply(out, "[[", "Gtest")
  } else {
    if(is.null(M)) {
      stop("M needs to be specified")
    }
    
    Gobs <- vector("list", M)
    Gtest <- vector("list", M)
    
    for(m in 1:M) {
      
      if(method == "random") {  # stratified sample of binary labels, sampling the edges and non-edges separately
        
        if(propnEdgesInclDiag) { # include diagonal (self-loops) in edge counts
          
          edges <- which(labels == 1)
          nonEdges <- which(labels == 0)
          trainEdges <- c(sample(edges,floor(rho*length(edges))),
                          sample(nonEdges,floor(rho*length(nonEdges))))
          
        } else {
          
          labelsNonDiag <- labels
          diag(labelsNonDiag) <- NA
          edgesNonDiag <- which(labelsNonDiag == 1)
          nonEdgesNonDiag <- which(labelsNonDiag == 0)
          diagIdx <- seq(1, ncol(labels) ^ 2, by = ncol(labels)+1)
          trainEdges <- c(sample(edgesNonDiag, floor(rho * length(edgesNonDiag))),
                          sample(nonEdgesNonDiag, floor(rho * length(nonEdgesNonDiag))),
                          sample(diagIdx, floor(rho * length(diagIdx))))
        }
        
      } else if(method == "rows") { # whole rows sampled together, stratified by rows with and without at least propnEdgesRowLB edges
        
        if(is.null(propnEdgesRowLB)) {
          stop("need input propnEdgeRowLB for method=rows")
        }
        
        nRow <- nrow(labels)
        nCol <- ncol(labels)
        
        if(propnEdgesInclDiag) { # include diagonal (self-loops) in edge counts

          nEdgesRow <- apply(labels, 1, sum) # number of edges in each row
          nNonEdgesRow <- apply(labels == 0, 1, sum) # number of non-edges in each row
          
          aboveLB <- which(nEdgesRow >= propnEdgesRowLB * nCol) # index for rows that meet the #edges threshold
          belowLB <- setdiff(1:nRow, aboveLB) # index for rows that do not meet the threshold

        } else {

          labelsNonDiag <- labels
          diag(labelsNonDiag) <- NA
          nEdgesRow <- apply(labelsNonDiag, 1, sum, na.rm = TRUE) # number of edges in each row
          nNonEdgesRow <- apply(labelsNonDiag == 0, 1, sum, na.rm = TRUE) # number of non-edges in each row
          
          aboveLB <- which(nEdgesRow >= propnEdgesRowLB * (nCol - 1)) # index for rows that meet the #edges threshold
          belowLB <- setdiff(1:nRow, aboveLB) # index for rows that do not meet the threshold

        }
        
        finalPropnEdgesTrain <- -1
        finalPropnNonEdgesTrain <- -1
        count <- 0
        
        while(abs(finalPropnEdgesTrain - rho) > 0.025 || abs(finalPropnNonEdgesTrain - rho) > 0.025) {
          # resample labels until sample has proportion rho (to within 0.025) of the edges and non-edges
          
          count <- count + 1
          trainRows <- c(sample(aboveLB, ceiling(rho * length(aboveLB))),
                         sample(belowLB, floor(rho * length(belowLB)))) # sample rows, stratified by whether rows meet #edges thrshold or not
          
          finalPropnEdgesTrain <- sum(nEdgesRow[trainRows]) / sum(nEdgesRow) # update propn edges that have been sampled
          finalPropnNonEdgesTrain <- sum(nNonEdgesRow[trainRows]) / sum(nNonEdgesRow) # update propn non-edges that have been sampled
        }
        
        trainEdges <- which(row(labels) %in% trainRows) # index for labels sampled to use as observed labels
        
        
      } else {
        
        stop("invalid input 'method': should be 'random' or 'rows'")
      }
      
      Gobs[[m]] <- labels
      Gobs[[m]][-trainEdges] <- NA # observed label matrix
      
      Gtest[[m]] <- labels
      Gtest[[m]][trainEdges] <- NA # unobserved label matrix to use as test data
    }
  }
  return(list(Gobs = Gobs,Gtest = Gtest))
} 


generateSubsampledData <- function(dat, n, geneIdxInt = NULL, M = NULL, fixed = NULL,
                                   exclude = NULL, scaled = FALSE, interventions = "non-nodes") {
  
  if(!is.list(dat) || is.data.frame(dat)) {
    if(is.null(M)) {
      stop("M needs to be specified")
    }
    dat <- rep(list(dat), M)
  } else {
    M <- length(dat)
  }
  if(!is.list(fixed)) {
    fixed <- list(fixed)
  }
  if(!is.list(exclude)) {
    exclude <- list(exclude)
  }
  if(!is.list(geneIdxInt)) {
    geneIdxInt <- list(geneIdxInt)
  }

  datSubsample <- mapply(function(d, f, e, g, n) {
    if(length(setdiff(f, 1:nrow(d))) > 0) {
      stop("Fixed sample index out of range")
    }
    if(length(setdiff(e, 1:nrow(d))) > 0) {
      stop("Excluded sample index out of range")
    }
    
    nSubsample <- n - length(f)
    if(nSubsample < 0) {
      stop(paste0("n must be larger than or equal to ", length(f),
                  ", the number of fixed samples that are always in the subsample."))
    }
    if(n > nrow(d) - length(e)) {
      stop(paste0("n must be smaller than or equal to ", nrow(d) - length(e),
                  ", the number of samples minus the number excluded from sampling."))
    }
    
    if(length(intersect(f, e)) > 1) {
      stop("fixed and exclude arguments should not have any overlap.")
    }
    
    dFixed <- d[f, ]
    nf <- setdiff(1:nrow(d), c(f, e))
    dNonfixed <- d[nf, ]
    if(length(nf) > nSubsample) {
      sampleIdx <- sample(nrow(dNonfixed), nSubsample)
    } else if(length(nf) == nSubsample) {
      sampleIdx <- 1:nSubsample  
    } else {
      stop("nSubsample is larger than number of samples to sample from")
    }
    
    return(list(dat = rbind(dFixed, dNonfixed[sampleIdx, ]),
                geneIdxInt=g[c(f, nf[sampleIdx])]))
  }, d = dat, f = fixed, e = exclude, g = geneIdxInt, n = n, SIMPLIFY = FALSE)
  
  if(scaled) {
    datSubsample <- lapply(datSubsample, function(x) {
      x$dat <- scale(x$dat)
      return(x)
    })
  }
  
  return(datSubsample)
}


initRes <- function(params) {
  
  initResSingle <- function(M, nIterMethod = NULL, pointEst = FALSE) {
    # initialization of results list for a single method
    
    Ghat <- vector("list", M)
    if(pointEst) { # methods that return a binary adjacency matrix
      pointEval <- vector("list", M)
      return(list(Ghat = Ghat, pointEval = pointEval))
    } else { # methods that return a weighted adjacency matrix
      roc <- vector("list", M)
      if(is.null(nIterMethod)) {
        auc <- vector(length = M)
      } else {
        auc <- matrix(NA, nrow = nIterMethod, ncol = M)
      }
      return(list(Ghat = Ghat, roc = roc, auc = auc))
    }
  }
   
  res <- list()
  
  if(params$meth$doMRCL) {
    res$MRCL <- initResSingle(params$meth$M, params$meth$M_MRCL)
    res$MRCL$feat <- vector("list", length = params$meth$M)
  } 
  
  if(params$meth$doCor) {
    res$CorP <- initResSingle(params$meth$M)
    res$CorK <- initResSingle(params$meth$M)
  } 
  
  if(params$meth$doLasso) {
    res$Lasso <- initResSingle(params$meth$M,params$meth$M_Lasso)
  } 
  
  if(params$meth$doIDA) {
    
    methodNames <- c("IDA","IDAcons")
    for(i in 1:length(methodNames)) {
      res[[methodNames[i]]] <- initResSingle(params$meth$M)
    }
    methodNames <- c("PCstrict", "PCloose", "PCconsStrict", "PCconsLoose",
                     "PCstrictTC", "PClooseTC", "PCconsStrictTC", "PCconsLooseTC")
    for(i in 1:length(methodNames)) {
      res[[methodNames[i]]] <- initResSingle(params$meth$M, pointEst=TRUE)
    }
    
  } 
  
  if(params$meth$doRFCI) {
    
    methodNames <- c("RFCIstrict", "RFCIloose", "RFCIconsStrict", "RFCIconsLoose",
                     "RFCIstrictTC", "RFCIlooseTC", "RFCIconsStrictTC", "RFCIconsLooseTC")
    for(i in 1:length(methodNames)) {
      res[[methodNames[i]]] <- initResSingle(params$meth$M, pointEst = TRUE)
    }
    res$RFCIstrict$PAG <- vector("list", length = params$meth$M)
    res$RFCIconsStrict$PAG <- vector("list", length = params$meth$M)
    
  }
    
  if(params$meth$doGIES) {
    
    methodNames <- c("GIESstrict", "GIESloose", "GIESconsStrict", "GIESconsLoose",
                     "GIESstrictTC", "GIESlooseTC", "GIESconsStrictTC", "GIESconsLooseTC")
    for(i in 1:length(methodNames)) {
      res[[methodNames[i]]] <- initResSingle(params$meth$M, pointEst = TRUE)
    }
  }
  
  if(params$meth$doKNN) {
    res$KNN <- initResSingle(params$meth$M)
    res$KNN$modData <- vector("list", params$meth$M)
  }
  
  return(list(res = res, params = params))
}
  

runMRCL <- function(X, Gobs, M = 1, lambda1 = 0.001, lambda2 = 0.001, outOfRange = "ignore") {
  
  Ghat <- vector("list", M)
  for(m in 1:M) {
    Gobs[Gobs == 0] <- -1
    out <- mrcl(X = X, Gobs = Gobs, lambda1 = lambda1, lambda2 = lambda2, outOfRange = outOfRange)
    Ghat[[m]] <- out$Ghat[[1]]$Ghat
    dimnames(Ghat[[m]]) <- list(colnames(X), colnames(X))
    if(m == 1) {
      feat <- out$feat
    }
  }
  return(list(Ghat = Ghat, feat = feat))
}


runCor <- function(X, method) {
  Ghat <- abs(cor(X, method = method))
  return(Ghat)
}


runLasso <- function(X, M) {
  Ghat <- vector("list", M)
  for(m in 1:M) {
    Ghat[[m]] <- matrix(nrow = ncol(X), ncol = ncol(X))
    dimnames(Ghat[[m]]) <- list(colnames(X), colnames(X))
    for(j in 1:ncol(X)) {
      cvg <-cv.glmnet(as.matrix(X[, -j]), X[, j])
      Ghat[[m]][-j, j] <- abs(coef(cvg)[-1])
    }
  }
  return(Ghat)
}


runIDA <- function(X, sigLevel, Gobs = NULL) {
  GhatIDA <- matrix(nrow = ncol(X), ncol = ncol(X))
  dimnames(GhatIDA) <- list(colnames(X), colnames(X))
  
  suffStat <- list(C = cor(X), n = nrow(X))
  
  if(!is.null(Gobs)) {
    fixedGaps <- (Gobs == 0) & t(Gobs == 0)
    fixedGaps[is.na(fixedGaps)] <- FALSE
    dimnames(fixedGaps) <- NULL
    fixedEdges <- (Gobs == 1) | t(Gobs == 1)
    fixedEdges[is.na(fixedEdges)] <- FALSE
    dimnames(fixedEdges) <- NULL
  } else {
    fixedGaps <- NULL
    fixedEdges <- NULL
  }
  pc.fit <- pc(suffStat, indepTest = gaussCItest, alpha = sigLevel, p = ncol(suffStat$C),
               fixedGaps = fixedGaps, fixedEdges = fixedEdges, skel.method = "stable.fast")
  GhatPCloose <- t(as(pc.fit, "matrix"))
  dimnames(GhatPCloose) <- list(colnames(X), colnames(X))
  GhatPCstrict <- cpdag2adj(GhatPCloose,method = "strict")
  GhatPClooseTC <- adj2transitiveClosure(GhatPCloose)
  GhatPCstrictTC <- adj2transitiveClosure(GhatPCstrict)
  
  for (j in 1:ncol(X)){
    eff <- idaFast(j, 1:ncol(X), cov(X), pc.fit@graph) 
    GhatIDA[j, ] <- apply(abs(eff), 1, min)
  }
  return(list(GhatIDA = GhatIDA,
              GhatPCloose = GhatPCloose,GhatPClooseTC = GhatPClooseTC,
              GhatPCstrict = GhatPCstrict,GhatPCstrictTC = GhatPCstrictTC))
}


runRFCI <- function(X, sigLevel, Gobs = NULL) {
  # Based on code from Umberto Noe 
  
  if(!is.null(Gobs)) {
    fixedGaps <- (Gobs == 0) & t(Gobs == 0)
    fixedGaps[is.na(fixedGaps)] <- FALSE
    dimnames(fixedGaps) <- NULL
    fixedEdges <- (Gobs == 1) | t(Gobs == 1)
    fixedEdges[is.na(fixedEdges)] <- FALSE
    dimnames(fixedEdges) <- NULL
  } else {
    fixedGaps <- NULL
    fixedEdges <- NULL
  }
  
  suffStat  <- list(C = cor(X), n = nrow(X))
  indepTest <- gaussCItest
  res       <- rfci(suffStat, indepTest, alpha = sigLevel, p = ncol(suffStat$C),
                    fixedGaps = fixedGaps, fixedEdges = fixedEdges, skel.method = "stable.fast")
  
  PAG <- as(res, "matrix")
  dimnames(PAG) <- list(colnames(X), colnames(X))
  
  GhatStrict <- pag2adj(PAG, "strict")
  GhatStrictTC <- adj2transitiveClosure(GhatStrict)
  GhatLoose  <- pag2adj(PAG, "loose")
  GhatLooseTC  <- adj2transitiveClosure(GhatLoose)
  
  Output <- list( PAG = PAG,
                  GhatStrict = GhatStrict, GhatStrictTC = GhatStrictTC,
                  GhatLoose  = GhatLoose, GhatLooseTC = GhatLooseTC )
  return(Output)
  
}


runGIES <- function(X, targets, Gobs = NULL, maxDegree = integer(0)) {
  # Based on code from Umberto Noe
  
  if(!is.null(Gobs)) {
    fixedGaps <- (Gobs == 0) & t(Gobs == 0)
    fixedGaps[is.na(fixedGaps)] <- FALSE
    dimnames(fixedGaps) <- NULL
  } else {
    fixedGaps <- NULL
  }
  
  if (all(is.na(targets))) {
    # obs data only
    targets      <- list( integer(0) )
    target.index <- rep(1, nrow(X))
  } else {
    # obs and int data
    obsIdx <- is.na(targets)
    targetsNew <- c(list(integer(0)), as.list(targets[!obsIdx]))
    target.index <- targets
    target.index[obsIdx] <- 1
    target.index[!obsIdx] <- 1 + seq_len(sum(!obsIdx))
    
  }
  
  score <- new("GaussL0penIntScore", data = X, targets = targetsNew,
               target.index = target.index, nodes = colnames(X))
  
  gies.fit <- gies(score, fixedGaps = fixedGaps, maxDegree = maxDegree)
  
  CPDAG <- as(gies.fit$essgraph, "matrix") * 1
  dimnames(CPDAG) <- list(colnames(X), colnames(X))
  
  GhatStrict <- cpdag2adj(CPDAG, "strict")
  GhatStrictTC <- adj2transitiveClosure(GhatStrict)
  GhatLoose  <- CPDAG
  GhatLooseTC  <- adj2transitiveClosure(GhatLoose)
  
  Output <- list(cpdag = CPDAG, GhatStrict = GhatStrict, GhatStrictTC = GhatStrictTC,
                 GhatLoose  = GhatLoose, GhatLooseTC = GhatLooseTC)
  return( Output )
  
}


runKNN <- function(X, Gobs, q, Ktune, feat = NULL) {
  
  Ghat <- matrix(nrow = ncol(X), ncol = ncol(X))
  dimnames(Ghat) <- list(colnames(X), colnames(X))
  
  # bivariate featurisation using MRCL featurization
  if(is.null(feat) || ncol(feat) < q) {
    feat <- compute_featurization(X, Gobs, binsRange = c(-3, 3), outOfRange = "ignore", q = q)$postPCA
  } else {
    feat <- feat[, 1:q]
  }
  Y <- as.matrix( c(t(Gobs)) )
  Yidx <- cbind(c(t(row(Gobs))), c(t(col(Gobs))))
  
  Ytrain <- Y[!is.na(Y)]
  Ytrain <- factor(Ytrain, levels = c(1, 0), labels = c("edge", "nonEdge"))
  featTrain <- data.frame(feat[!is.na(Y), ])
  featTest <- data.frame(feat[is.na(Y), ])
  ctrl <- trainControl(method = "repeatedcv", number = 5, repeats = 5, returnResamp="all",
                       savePredictions = TRUE, classProbs = TRUE, 
                       summaryFunction = twoClassSummary, allowParallel = FALSE)
  
  mod <- train(featTrain, Ytrain, method = "knn", trControl = ctrl, preProcess = c("center", "scale"),
               tuneGrid = Ktune, metric = "ROC")
  preds <- predict(mod, newdata = featTest, type = "prob")
  Yidx <- Yidx[is.na(Y), ]
  for(i in 1:nrow(Yidx)) {
    Ghat[Yidx[i, 1], Yidx[i, 2]] <- preds[i, 1]
  }
  
  return(list(Ghat = Ghat, mod = mod, feat = feat))
}


pag2adj <- function(G, method = "strict") {
  # Based on code from Umberto Noe 
  method <- match.arg(method, c("strict", "loose"))
  
  p <- ncol(G)
  Gtilde <- G * 0
  
  if (method == "strict") {
    #
    # only keep -->
    #
    for (i in 1:p) {
      for (j in 1:p) {
        arrow <- G[i, j] == 2 & G[j, i] == 3
        if (arrow) { 
          Gtilde[i, j] <- 1 
        }
      }
    }
  } else {
    #
    # --> or o-> or o-o
    #
    for (i in 1:p) {
      for (j in 1:p) {
        arrow <- G[i, j] == 2 & G[j, i] == 3
        circleArrow <- G[i, j] == 2 & G[j, i] == 1
        circleCircle <- G[i, j] == 1 & G[j, i] == 1
        
        if (arrow | circleArrow | circleCircle) { 
          Gtilde[i, j] <- 1
        }
      }
    }
  }
  
  return(Gtilde)
}


cpdag2adj <- function(G, method = "strict") {
  #  Based on code from Umberto Noe 
  method <- match.arg(method, c("strict", "loose"))
  
  if (method == "strict") {
    #
    # only keep -->
    #
    p <- ncol(G)
    Gtilde <- G * 0
    
    for (i in 1:p) {
      for (j in 1:p) {
        arrow <- G[i,j] == 1 & G[j,i] == 0
        if (arrow) { 
          Gtilde[i,j] <- 1 
        }
      }
    }
    
  } else {
    #
    # --> or --- (return the CPDAG)
    #
    Gtilde <- G
  }
  
  return(Gtilde)
}


adj2transitiveClosure <- function(G) {
    n <- nrow(G)
    Gtc <- G
    Gexp <- G
    for(i in 1:(n - 2)) {
      Gexp <- Gexp %*% G 
      Gtc <- Gtc + Gexp
    }
    diag(Gtc) <- 0
    return((Gtc > 0) * 1)
}


runMethods <- function(dataTrain, methPar, res, Gobs, Gtest, geneIdxNode = NULL) {
  
  for(m in 1:methPar$M) {
    
    if(is.matrix(dataTrain[[m]]) || is.data.frame(dataTrain[[m]])) {
      X <- dataTrain[[m]]
    } else if(is.list(dataTrain[[m]]) && !is.null(dataTrain[[m]]$dat)) {
      X <- dataTrain[[m]]$dat
    } else {
      stop("invalid dataTrain structure")
    }
    
    # MRCL -----------------------------
    if(methPar$doMRCL) { 
      
      # run MRCL M_MRCL times
      out <- runMRCL(X = X, Gobs = Gobs[[m]], M = methPar$M_MRCL, 
                     lambda1 = methPar$lambda1_MRCL, lambda2 = methPar$lambda2_MRCL)
      res$MRCL$Ghat[[m]] <- out$Ghat
      res$MRCL$feat[[m]] <- out$feat
      
      # evaluation 
      res$MRCL$roc[[m]] <- lapply(res$MRCL$Ghat[[m]], evaluateScores, Gtest[[m]])
      res$MRCL$auc[, m] <- sapply(res$MRCL$roc[[m]], "[[", "auc")
    }
    
    # correlations -----------------------------
    if(methPar$doCor) { 
      
      # calculate correlations
      res$CorP$Ghat[[m]] <- runCor(X, method = "pearson")
      res$CorK$Ghat[[m]] <- runCor(X, method = "kendall")
      
      # evaluation
      res$CorP$roc[[m]] <- evaluateScores(res$CorP$Ghat[[m]], Gtest[[m]])
      res$CorK$roc[[m]] <- evaluateScores(res$CorK$Ghat[[m]], Gtest[[m]])
      
      res$CorP$auc[m] <- res$CorP$roc[[m]]$auc
      res$CorK$auc[m] <- res$CorK$roc[[m]]$auc
    }
    
    # Lasso -------------------------------
    if(methPar$doLasso) { 
      
      # calculate regression coefficients
      res$Lasso$Ghat[[m]] <- runLasso(X, methPar$M_Lasso)
      
      # evaluation
      res$Lasso$roc[[m]] <- lapply(res$Lasso$Ghat[[m]], evaluateScores, Gtest[[m]])
      res$Lasso$auc[,m] <- sapply(res$Lasso$roc[[m]], "[[", "auc")
    }
    
    # IDA ---------------------------
    if(methPar$doIDA) {  
    
      # calculate causal effects
      tmp <- runIDA(X, methPar$sigLevelIDA)
      res$IDA$Ghat[[m]] <- tmp$GhatIDA
      methodNames1 <- c("PCstrict", "PCloose","PCstrictTC", "PClooseTC")
      for(i in 1:length(methodNames1)) {
        res[[methodNames1[i]]]$Ghat[[m]] <- tmp[[paste0("Ghat", methodNames1[i])]]
      }
      
      tmp <- runIDA(X, methPar$sigLevelIDA, Gobs = Gobs[[m]])
      res$IDAcons$Ghat[[m]] <- tmp$GhatIDA
      methodNames2 <- c("PCconsStrict", "PCconsLoose","PCconsStrictTC", "PCconsLooseTC")
      for(i in 1:length(methodNames2)) {
        res[[methodNames2[i]]]$Ghat[[m]] <- tmp[[paste0("Ghat", methodNames1[i])]]
      }
      rm(tmp)
      
      # evaluation
      res$IDA$roc[[m]] <- evaluateScores(res$IDA$Ghat[[m]], Gtest[[m]])
      res$IDA$auc[m] <- res$IDA$roc[[m]]$auc
      res$IDAcons$roc[[m]] <- evaluateScores(res$IDAcons$Ghat[[m]], Gtest[[m]])
      res$IDAcons$auc[m] <- res$IDAcons$roc[[m]]$auc
      
      methodNames <- c("PCstrict", "PCloose", "PCconsStrict", "PCconsLoose",
                       "PCstrictTC", "PClooseTC", "PCconsStrictTC", "PCconsLooseTC")
      for(i in 1:length(methodNames)) {
        res[[methodNames[i]]]$pointEval[[m]] <- evaluatePoint(res[[methodNames[i]]]$Ghat[[m]], Gtest[[m]])
      }
    }
    
    # RFCI ---------------------------
    if(methPar$doRFCI) {  

      tmp <- runRFCI(X, methPar$sigLevelRFCI)
      res$RFCIstrict$Ghat[[m]] <- tmp$GhatStrict
      res$RFCIstrictTC$Ghat[[m]] <- tmp$GhatStrictTC
      res$RFCIloose$Ghat[[m]] <- tmp$GhatLoose
      res$RFCIlooseTC$Ghat[[m]] <- tmp$GhatLooseTC
      res$RFCIstrict$PAG[[m]] <- tmp$PAG
      
      tmp <- runRFCI(X, methPar$sigLevelRFCI, Gobs=Gobs[[m]])
      res$RFCIconsStrict$Ghat[[m]] <- tmp$GhatStrict
      res$RFCIconsStrictTC$Ghat[[m]] <- tmp$GhatStrictTC
      res$RFCIconsLoose$Ghat[[m]] <- tmp$GhatLoose
      res$RFCIconsLooseTC$Ghat[[m]] <- tmp$GhatLooseTC
      res$RFCIconsStrict$PAG[[m]] <- tmp$PAG
      rm(tmp)
      
      # evaluation
      methodNames <- c("RFCIstrict", "RFCIloose", "RFCIconsStrict", "RFCIconsLoose",
                       "RFCIstrictTC", "RFCIlooseTC", "RFCIconsStrictTC", "RFCIconsLooseTC")
      for(i in 1:length(methodNames)) {
        res[[methodNames[i]]]$pointEval[[m]] <- evaluatePoint(res[[methodNames[i]]]$Ghat[[m]], Gtest[[m]])
      }
    }
    
    # GIES -----------------------------
    if(methPar$doGIES) {  
      
      targets <- match(dataTrainSubsample[[m]]$geneIdxInt, geneIdxNode[[m]])
      
      # calculate interventional essential graph
      if(methPar$maxDegreeGIES == 0) {
        tmp1 <- runGIES(X, targets)
        tmp2 <- runGIES(X, targets, Gobs = Gobs[[m]])
      } else {
        tmp1 <- runGIES(X, targets, maxDegree = methPar$maxDegreeGIES)
        tmp2 <- runGIES(X, targets, maxDegree = methPar$maxDegreeGIES, Gobs = Gobs[[m]])  
      }
      
      res$GIESstrict$Ghat[[m]] <- tmp1$GhatStrict
      res$GIESloose$Ghat[[m]] <- tmp1$GhatLoose
      res$GIESstrictTC$Ghat[[m]] <- tmp1$GhatStrictTC
      res$GIESlooseTC$Ghat[[m]] <- tmp1$GhatLooseTC
      
      res$GIESconsStrict$Ghat[[m]] <- tmp2$GhatStrict
      res$GIESconsLoose$Ghat[[m]] <- tmp2$GhatLoose
      res$GIESconsStrictTC$Ghat[[m]] <- tmp2$GhatStrictTC
      res$GIESconsLooseTC$Ghat[[m]] <- tmp2$GhatLooseTC
      rm(tmp1, tmp2)
      
      # evaluation
      methodNames <- c("GIESstrict", "GIESloose", "GIESconsStrict", "GIESconsLoose",
                       "GIESstrictTC", "GIESlooseTC", "GIESconsStrictTC", "GIESconsLooseTC")
      for(i in 1:length(methodNames)) {
        res[[methodNames[i]]]$pointEval[[m]] <- evaluatePoint(res[[methodNames[i]]]$Ghat[[m]], Gtest[[m]])
      }
    }
    
    # KNN ----------------------------
    if(methPar$doKNN) {  

      res$KNN$modData[[m]] <- runKNN(X,Gobs[[m]],methPar$nPrCmpKNN,methPar$Ktune,feat=res$MRCL$feat[[m]])
      res$KNN$Ghat[[m]] <- res$KNN$modData[[m]]$Ghat
      
      # evaluation
      res$KNN$roc[[m]] <- evaluateScores(res$KNN$Ghat[[m]],Gtest[[m]])
      res$KNN$auc[m] <- res$KNN$roc[[m]]$auc
    }
  }
  
  return(res)
  
}


evaluateScores <- function(Ghat, Gtrue, direction = "<", algorithm = 1) {
  diag(Ghat) <- NA
  roc <- roc(c(Gtrue), c(data.matrix(Ghat)), direction = direction, algorithm = algorithm)
  # algorithm argument set to 1 for legacy reasons. This was previously the default,
  # but has been changed in the pROC package since this analysis was run.
}


evaluatePoint <- function(Ghat, Gtrue, nEdge = NULL) {
  
  diag(Ghat) <- NA
  diag(Gtrue) <- NA
  
  if(!is.null(nEdge)) {
    
    sortedScores <- sort(Ghat, decreasing = TRUE)
    if(nEdge > 0) {
      threshold <- sortedScores[nEdge]
    } else {
      threshold <- Inf
    }
    
    if(sum(sortedScores == threshold) > 1) {
      nEdgeActual <- sum(sortedScores >= threshold)
      if(nEdgeActual - nEdge > nEdge - sum(sortedScores > threshold)) {
        nEdge <- sum(sortedScores > threshold)
      } else {
        nEdge <- nEdgeActual
      }
    }
    if(nEdge>0) {
      threshold <- sortedScores[nEdge]
    } else {
      threshold <- Inf
    }
    Ghat <- Ghat>=threshold
  } else {
    threshold <- NULL
  }
  
  naIdx <- is.na(Ghat) | is.na(Gtrue)
  
  ntp <- sum(Ghat[!naIdx] == 1 & Gtrue[!naIdx] == 1)
  nfp <- sum(Ghat[!naIdx] == 1 & Gtrue[!naIdx] == 0)
  ntn <- sum(Ghat[!naIdx] == 0 & Gtrue[!naIdx] == 0)
  nfn <- sum(Ghat[!naIdx] == 0 & Gtrue[!naIdx] == 1)
  
  sens <- ntp / (ntp + nfn)
  spec <- ntn / (ntn + nfp)
  prec <- ntp / (ntp + nfp)

  return(list(sens = sens, spec = spec, prec = prec, threshold = threshold))
  
}


collateAUCs <- function(experiment, rhoSet = NULL, nSet = NULL, single_n_per_rho = FALSE) {
  
  if(is.null(rhoSet) && is.null(nSet)) { # cell line experiment
    load(file.path("output", experiment, "results", "results.RData"))
    aucs <- lapply(res, function(x) {
      if(!is.null(x$auc)) {
        if(is.matrix(x$auc)) {
          data.frame(cellLine = colnames(x$auc), auc = colMeans(x$auc), stringsAsFactors = FALSE)
        } else {
          data.frame(cellLine = names(x$auc), auc = x$auc, stringsAsFactors = FALSE)
        }
      }
    })
    
    cellLine <- aucs[[1]]$cellLine
    avgROCs <- lapply(cellLine, function(cl) {
      lapply(res, function(x) {
        if(!is.null(x$roc)) {
          if(class(x$roc[[cl]]) == "list") {
            avgROCs <- pROCavgROC(x$roc[[cl]], plot = FALSE)
          } else {
            avgROCs <- x$roc[[cl]] # no average to take
          }
          return(avgROCs)
        } else {
          return(NA)
        }
      })
    })
    names(avgROCs) <- cellLine

    return(list(aucs = aucs, avgROCs = avgROCs, params = params))
    
  } else { # yeast and TCPA experiments
  
    aucs <- NULL
    
    if(single_n_per_rho){
      nSetAll <- nSet
      avgROCs <- vector("list", length(rhoSet))
      names(avgROCs) = paste0("rho=", rhoSet, "_n=", nSetAll)
    } else {
      avgROCs <- matrix(list(), nrow = length(rhoSet), ncol = length(nSet),
                        dimnames = list(rhoSet, nSet))
    }
    
    for(rho in rhoSet) {
      if(single_n_per_rho) {
        nSet <- nSetAll[which(rho == rhoSet)]
      }
      for(n in nSet) {
        load(file.path("output", experiment, "results", paste0("n=", n, "_rho=", rho, ".RData")))
        
        # collate aucs, averaging over within-method iterations
        auc_rho_n <- lapply(res, function(x) {
          if(!is.null(x$auc)) {
            if(is.matrix(x$auc)) {
              data.frame(rho = rho, n = params$dataSamp$n, iteration = 1:params$meth$M,
                         auc = colMeans(x$auc), stringsAsFactors = FALSE)
            } else {
              data.frame(rho = rho, n = params$dataSamp$n, iteration = 1:params$meth$M,
                         auc = x$auc, stringsAsFactors = FALSE)
            }
          }
        })
        
        if(is.null(aucs)) {
          aucs <- auc_rho_n
        } else {
          aucs <- mapply(rbind, aucs, auc_rho_n, SIMPLIFY=FALSE)
        }
        
        # avg ROC
        avgROCs_rho_n <- lapply(res, function(x) {
          if(!is.null(x$roc)) {
            if(class(x$roc[[1]]) == "list") {
              rocData <- unlist(x$roc, recursive = FALSE)
            } else {
              rocData <- x$roc
            }
            return(pROCavgROC(rocData, plot = FALSE))
          } else {
            return(NA)
          }
        })
        if(single_n_per_rho) {
          avgROCs[[paste0("rho=", rho, "_n=", n)]] <- avgROCs_rho_n
        } else {
          avgROCs[[as.symbol(rho), as.symbol(n)]] <- avgROCs_rho_n
        }      
      }
    }
    
    if(single_n_per_rho) {
      params$dataSamp$n <- as.numeric(nSetAll)
    } else {
      params$dataSamp$n <- as.numeric(unique(aucs[[1]]$n))  
    }
    
    params$labelSplit$rho <- rhoSet
  
    return(list(aucs = aucs, avgROCs = avgROCs, params = params))
  }
}


pROCavgROC <- function(x, plot = TRUE, auc = TRUE) {
  controlsCollated <- unlist(lapply(x, "[[", "controls"))
  casesCollated <- unlist(lapply(x, "[[", "cases"))
  newThresholds <- pROC:::roc.utils.thresholds(c(controlsCollated, casesCollated),x[[1]]$direction)
  
  newMetrics <- lapply(x,function(xi) {
    pROC:::roc.utils.perfs.all.safe(newThresholds, xi$controls, xi$cases, xi$direction)
  })
  
  avgSe <- rowMeans(sapply(newMetrics, "[[", "se"))
  avgSp <- rowMeans(sapply(newMetrics, "[[", "sp"))
  
  avgROC <- list()
  class(avgROC) <- "roc"
  avgROC$percent <- x[[1]]$percent
  avgROC$sensitivities <- avgSe
  avgROC$specificities <- avgSp
  avgROC$thresholds <- newThresholds
  
  if(plot) {
    pROC:::plot.roc.roc(x[[1]])
    for(i in 2:length(x)) {
      pROC:::plot.roc.roc(x[[i]], add = TRUE)
    }
    pROC:::plot.roc.roc(avgROC, add = TRUE, col = "red")
  }
  
  if(auc) {
    avgROC$auc <- as.numeric(pROC:::auc.roc(avgROC))
  }
  
  return(avgROC)
  
}


aucMeanSd <- function(aucs) { 
  # function to calculate means and standard deviations for each combination of n and rho
  aucs <- aucs[!sapply(aucs, is.null)]
  aucMeans <- lapply(aucs, function(x) aggregate(auc ~ rho + n, data = x, FUN = mean))
  aucMeans <- lapply(aucMeans, setNames, c("rho", "n", "aucMean"))
  aucMeans <- bind_rows(aucMeans, .id = "method")
  
  aucSDs <- lapply(aucs, function(x) aggregate(auc ~ rho + n, data = x, FUN = sd))
  aucSDs <- lapply(aucSDs, setNames, c("rho", "n", "aucSD"))
  aucSDs <- bind_rows(aucSDs, .id = "method")
  
  aucSummary <- merge(aucMeans, aucSDs)
}


collatePointEstimateScores <- function(experiment, listEntryName, rhoSet = NULL, 
                                       nSet = NULL, single_n_per_rho = FALSE) {
  
  if(is.null(rhoSet) & is.null(nSet)) { # cellLine experiment
    
    load(file.path("output", experiment, "results", "results.RData"))
    
    # collate scores
    scores <- lapply(res, function(x) {
      tmp <- do.call("rbind", lapply(x$pointEval, unlist))
      if(is.null(tmp) && all(sapply(x, is.list))) {
        tmp <- lapply(x, function(y) do.call("rbind", lapply(y$pointEval, unlist)))
      }
      if(!is.null(tmp)) {
        if(is.list(tmp)) {
          lapply(tmp, function(y) data.frame(cellLine = rownames(y), y, 
                                             stringsAsFactors = FALSE))
        } else {
          data.frame(cellLine = rownames(tmp), tmp, stringsAsFactors = FALSE)
        }
      }
    })
    
    scores <- scores[!sapply(scores,is.null)]
    
  } else { # other experiments
    
    scores <- NULL
    
    if(single_n_per_rho){
      nSetAll <- nSet
    }
    
    for(rho in rhoSet) {
      if(single_n_per_rho) {
        nSet <- nSetAll[which(rho == rhoSet)]
      }
      for(n in nSet) {
        load(file.path("output", experiment, "results", paste0("n=", n, "_rho=", rho, ".RData")))
        
        # collate scores
        scores_rho_n <- lapply(res, function(x) {
          tmp <- do.call("rbind", lapply(x$pointEval, unlist))
          if(is.null(tmp) && all(sapply(x, is.list))) {
            tmp <- lapply(x, function(y) do.call("rbind", lapply(y$pointEval, unlist)))
          }
          if(!is.null(tmp)) {
            if(is.list(tmp)) {
              lapply(tmp, function(y) data.frame(rho = rho, n = params$dataSamp$n, 
                                                 iteration = 1:params$meth$M, y, 
                                                 stringsAsFactors = FALSE))
            } else {
              data.frame(rho = rho, n = params$dataSamp$n, 
                         iteration = 1:params$meth$M, tmp, stringsAsFactors = FALSE)
            }
          }
        })
        
        expandOut <- !(sapply(scores_rho_n, is.null) | sapply(scores_rho_n, is.data.frame))
        scores_rho_n <- do.call(c, c(list(scores_rho_n[!expandOut]), scores_rho_n[expandOut]))
        
        if(is.null(scores)) {
          scores <- scores_rho_n
        } else {
          scores <- mapply(rbind, scores, scores_rho_n, SIMPLIFY = FALSE)
        }
        
      }
    }
    
    scores <- scores[!sapply(scores, is.null)]
    
    params$dataSamp$n <- unique(scores[[1]]$n)
    params$labelSplit$rho <- rhoSet
  }    
  return(list(scores = scores, params = params))
}


pointEstimateScoresMeanSd <- function(scores) {
  # function to calculate means and standard deviations for each combination of n and rho
  scoreMeans <- lapply(scores, function(x) aggregate(cbind(sens, spec, prec) ~ rho + n,
                                                     data = x, FUN = mean, 
                                                     na.rm = TRUE, na.action = NULL))
  scoreMeans <- lapply(scoreMeans, setNames, c("rho", "n", 
                                               paste0(colnames(scoreMeans[[1]][-(1:2)]), "Mean")))
  scoreMeans <- bind_rows(scoreMeans, .id="method")
  
  scoreSDs <- lapply(scores, function(x) aggregate(cbind(sens, spec, prec) ~ rho + n,
                                                   data = x, FUN = sd,
                                                   na.rm = TRUE, na.action = NULL))
  scoreSDs <- lapply(scoreSDs, setNames, c("rho", "n",
                                           paste0(colnames(scoreSDs[[1]][-(1:2)]), "SD")))
  scoreSDs <- bind_rows(scoreSDs, .id="method")
  
  scoreSummary <- merge(scoreMeans, scoreSDs)
}


makeAUClinePlotYeastTCPA <- function(experiment, ylims) {
  
  load(file.path("output", experiment, "results", "collatedResults.RData"))

  dat <- aucsMeanSd
  dat$nLabel <- factor(dat$n, labels = paste0("bold(italic(n)[train] == ", levels(factor(dat$n)), ")"))
  
  dat <- dat[dat$method!="IDAcons", ]
  dat$Method <- factor(dat$method, levels = c("MRCL", "CorP", "CorK", "Lasso", "IDA", "KNN"), 
                       labels = c("MRCL", "Pearson", "Kendall", "Lasso", "IDA", "k-NN"))
 dat$method <- NULL
  
  pd <- position_dodge(0.06) 
  gg <- ggplot(dat, aes(x = rho, y = aucMean, colour = Method)) + 
    geom_errorbar(aes(ymin = aucMean - aucSD / sqrt(params$meth$M), 
                      ymax = aucMean + aucSD / sqrt(params$meth$M)), width = .1, position = pd) +
    geom_line(position = pd) +
    geom_point(position = pd) +
    facet_wrap(~nLabel, labeller = label_parsed) +
    scale_x_continuous(name = expression(rho), breaks = params$labelSplit$rho) +
    ylab("AUC") +
    ylim(ylims) +
    scale_colour_manual(name = "Method",
                        values = c(brewer.pal(5, "Set1")[c(1, 2, 4, 5, 3)], "#EEEE00")) +
    theme_bw(base_size = 15) +
    theme(text = element_text(color = "black"),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(face = "plain", size = rel(1.3)),
          axis.text = element_text(size = rel(1), colour = "black"),
          legend.title = element_text(face = "plain"),
          legend.text = element_text(size = rel(1)),
          legend.background = element_rect(fill = "transparent"),
          legend.key.size = unit(1, 'lines'),
          legend.position = "bottom",
          legend.direction = "horizontal",
          panel.border = element_rect(colour = "black"),
          panel.grid.minor = element_line(colour = NA),
          strip.text = element_text(face = "bold")
    )
  
  return(gg)
}


makeROCcurvePlotYeastTCPA <- function(experiment, ylims, xlims = NULL, pointMethodSet = "final") {
  
  load(file.path("output", experiment, "results", "collatedResults.RData"))
  nSet <- params$dataSamp$n
  if(experiment == "TCPA") {
    colnames(avgROCs)[colnames(avgROCs) == "max"] <- max(nSet)
  }
  
  if(is.null(xlims)) {
    xlims <- c(0, 1)
  }
  
  pointMethods <- c("PCloose", "PClooseTC", "PCstrict", "PCstrictTC",
                    "PCconsLoose", "PCconsLooseTC", "PCconsStrict", "PCconsStrictTC", 
                    "RFCIloose", "RFCIlooseTC", "RFCIstrict", "RFCIstrictTC", 
                    "RFCIconsLoose", "RFCIconsLooseTC", "RFCIconsStrict", "RFCIconsStrictTC", 
                    "GIESloose", "GIESlooseTC", "GIESstrict", "GIESstrictTC", 
                    "GIESconsLoose", "GIESconsLooseTC", "GIESconsStrict", "GIESconsStrictTC")
  
  pointMethodsLabels <- c("PC", "PC (TC)", "PC (strict)", "PC (strict;TC)", 
                    "PC (cnstrnts)", "PC (cnstrnts;TC)", "PC (cnstrnts;strict)", "PC (cnstrnts;strict;TC)", 
                    "RFCI", "RFCI (TC)", "RFCI (strict)", "RFCI (strict;TC)", 
                    "RFCI (cnstrnts)", "RFCI (cnstrnts;TC)", "RFCI (cnstrnts;strict)", "RFCI (cnstrnts;strict;TC)", 
                    "GIES", "GIES (TC)", "GIES (strict)", "GIES (strict;TC)", 
                    "GIES (cnstrnts)", "GIES (cnstrnts;TC)", "GIES (cnstrnts;strict)", "GIES (cnstrnts;strict;TC)")
  
  pointMethodsColor <- c(rep(rgb(0, 0, 0), 8), 
                         rep("#663600", 8), 
                         rep(rgb(0.5, 0.5, 0.5), 8))
  
  pointMethodsShape <- rep(c(1, 4, 3, 5, 0, 2, 6, 8), 3)
                    
  if(experiment == "yeast_row-wise_GIES") {
    
    rhoSet <- params$labelSplit$rho
    datROC <- lapply(avgROCs, function(x) {
      x <- x[sapply(x, class) == "roc"]
      x <- lapply(x, "[", c("sensitivities", "specificities"))
      bind_rows(x, .id = "Method")
    })
    datROC <- mapply(function(x, y, z) {
      cbind(n = y, rho = z, x)
    }, datROC, nSet, rhoSet, SIMPLIFY = FALSE)
    datROC <- bind_rows(datROC)
    colnames(datROC)[4:5] <- c("sensitivity", "specificity")
    
    if(pointMethodSet == "loose") {
      inclPointMethodsIdx <- c(1, 2, 9, 10, 17, 18)
    } else if(pointMethodSet == "strict") {
      inclPointMethodsIdx <- c(3, 4, 11, 12, 19, 20)
    } else if(pointMethodSet == "consLoose") {
      inclPointMethodsIdx <- c(5, 6, 13, 14, 21, 22)
    } else if(pointMethodSet == "consStrict") {
      inclPointMethodsIdx <- c(7, 8, 15, 16, 23, 24)
    } else if(pointMethodSet == "final") {
      if(substr(experiment, 1, 5) == "yeast") {
        inclPointMethodsIdx <- c(1, 2, 5, 6, 17, 18)
      } else if(experiment == "TCPA") {
        inclPointMethodsIdx <- c(1, 2, 5, 6, 17, 18)
      }
    }
  } else {
    
    datROC <- lapply(as.character(nSet), function(n) {
      dat_n <- lapply(avgROCs[, n], function(x) {
        x <- x[sapply(x, class) == "roc"]
        x <- lapply(x, "[", c("sensitivities", "specificities"))
        bind_rows(x, .id = "Method")
      })
      bind_rows(dat_n, .id = "rho")
    })
    names(datROC) <- nSet
    datROC <- bind_rows(datROC, .id = "n")
    colnames(datROC)[4:5] <- c("sensitivity", "specificity")
    
    if(pointMethodSet == "loose") {
      inclPointMethodsIdx <- c(1, 2, 9, 10)
    } else if(pointMethodSet == "strict") {
      inclPointMethodsIdx <- c(3, 4, 11, 12)
    } else if(pointMethodSet == "consLoose") {
      inclPointMethodsIdx <- c(5, 6, 13, 14)
    } else if(pointMethodSet == "consStrict") {
      inclPointMethodsIdx <- c(7, 8, 15, 16)
    } else if(pointMethodSet == "final") {
      if(substr(experiment, 1, 5) == "yeast") {
        inclPointMethodsIdx <- c(1, 2, 5, 6)
      } else if(experiment == "TCPA") {
        inclPointMethodsIdx <- c(1, 2, 5, 6)
      }
    }
  }
  
  datROC <- datROC[datROC$Method != "IDAcons", ]
  datROC$Method <- factor(datROC$Method, levels = c("MRCL", "CorP", "CorK", "Lasso", "IDA", "KNN", pointMethods[inclPointMethodsIdx]), 
                          labels = c("MRCL", "Pearson", "Kendall", "Lasso", "IDA", "k-NN", pointMethodsLabels[inclPointMethodsIdx]))
  
  datPoint <- scoresMeanSd[scoresMeanSd$method %in% pointMethods[inclPointMethodsIdx], c("n", "rho", "method", "sensMean", "specMean")]
  datPoint$method <- factor(datPoint$method, levels = pointMethods[inclPointMethodsIdx], labels = pointMethodsLabels[inclPointMethodsIdx])#,levels = c("MRCL","CorP","IDA","CorK","Lasso","KNN"),labels = c("MRCL","Pearson","IDA","Kendall","Lasso","KNN"))
  
  dat <- rbind(datROC, setNames(datPoint, names(datROC)))
  dat$n <- factor(dat$n, levels = nSet)
  dat$nLabel <- factor(dat$n, labels = paste0("bold(italic(n)[train] == ", levels(factor(dat$n)), ")"))
  dat$rhoLabel <- factor(dat$rho, labels = paste0("bold(rho == ", levels(factor(dat$rho)), ")"))
  
  
  gg <- ggplot(data = dat, aes(x = 1-specificity, y = sensitivity, colour = Method, shape = Method, linetype = Method)) +
    geom_line(size = 0.65) +
    geom_point(size = 2, stroke = 1.1) +
    scale_x_continuous(labels = c(0, "", 0.5, "", 1), limits = xlims) +
    scale_y_continuous(labels = c(0, "", 0.5, "", 1), limits = ylims) +
    scale_colour_manual(name = "Method", 
                        values = c(brewer.pal(5, "Set1")[c(1, 2, 4, 5, 3)], "#EEEE00", 
                                 pointMethodsColor[inclPointMethodsIdx])) +
    scale_shape_manual(name = "Method",
                       values = c(rep(NA, 6), pointMethodsShape[inclPointMethodsIdx])) +
    scale_linetype_manual(name = "Method", 
                       values = c(rep(1, 6), rep(0, length(inclPointMethodsIdx)))) +
    theme_bw(base_size = 15) +
    theme(text = element_text(color = "black"), 
          axis.line = element_line(colour = "black"), 
          axis.title = element_text(face = "plain", size = rel(1.1)),
          axis.text = element_text(size  =  rel(1), colour ="black"),
          legend.title = element_text(face = "plain"),
          legend.text = element_text(size = rel(1)),
          legend.background = element_rect(fill = "transparent"),
          legend.key.size = unit(1, 'lines'),
          legend.position = "bottom",
          legend.direction = "horizontal",
          panel.border = element_rect(colour = "black"),
          panel.grid.minor = element_line(colour = NA),
          strip.text  =  element_text(face = "bold")
    ) + 
    guides(colour = guide_legend(nrow = 4))
  
  if(experiment == "yeast_row-wise_GIES") {
    gg <- gg + facet_wrap(~rhoLabel+nLabel, labeller = label_parsed)    
  } else {
    gg <- gg + facet_grid(rhoLabel~nLabel, labeller = label_parsed)  
  }

  return(gg)
}


makeAUCbarPlotCellLine <- function(ylims) {
  
  load(file.path("output","cellLine", "results", "collatedResults.RData"))
  
  dat <- aucs
  dat$Method <- factor(dat$Method, levels = c("MRCL", "CorP", "CorK", "Lasso", "IDA", "KNN"), 
                       ordered = TRUE, labels = c("MRCL", "Pearson", "Kendall", "Lasso", "IDA", "k-NN"))
  dat <- dat[!is.na(dat$Method), ]
  
  gg <- ggplot(dat,  aes(x = Method, y = auc, fill = Method)) + 
    geom_col(position = "dodge") +
    ylab("AUC") +
    xlab("") +
    ylim(ylims) +
    facet_wrap(~cellLine) +
    scale_fill_manual(values = c(brewer.pal(5, "Set1")[c(1, 2, 4, 5, 3)], "#EEEE00")) +
    theme_bw(base_size = 15) +
    theme(text = element_text(color = "black"), 
          axis.line = element_line(colour = "black"), 
          axis.title.y = element_text(face = "plain",  size = rel(1.3)), 
          axis.title.x = element_blank(), 
          axis.text = element_text(size = rel(1),  colour = "black"), 
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
          legend.title = element_text(face = "plain"), 
          legend.text = element_text(), 
          legend.background = element_rect(fill = "transparent"), 
          legend.key.size = unit(1, 'lines'), 
          legend.position = "none",
          legend.direction = "horizontal",
          panel.border = element_rect(colour = "black"),
          panel.grid.minor = element_line(colour = NA),
          strip.text = element_text(face = "bold")
    )
  
  return(gg)
}