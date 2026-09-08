


EnumerateNAs <- function(quantTbl, idxVal=1) {
  
  # Enumerate NaNs per sample
  freqNA  <- do.call("rbind",apply(quantTbl, idxVal, function(vec) data.frame(N_MISSING=length(vec[is.na(vec) | is.nan(vec)]))))
  
  # Set feature names
  if (idxVal == 1) { freqNA["ID"] <- rownames(quantTbl) }
  if (idxVal == 2) { freqNA["ID"] <- colnames(quantTbl) }
  
  # Enumerate proportion
  if (idxVal == 1) { freqNA["PROP_MISSING"] <- freqNA$N_MISSING / ncol(quantTbl) }
  if (idxVal == 2) { freqNA["PROP_MISSING"] <- freqNA$N_MISSING / nrow(quantTbl) }
  
  # Return
  return(freqNA[,c("ID","N_MISSING","PROP_MISSING")])
  
}


TestProteinDropThresholds <- function(perProtNA, quantTbl, protThresholds=c(0.05, 0.1, 0.15, 0.2, 0.25, 1), idxVal=1) {
  
  # Get data
  dropProts <- data.frame()
  
  # Loop over thresholds
  for (protVal in protThresholds) {
    
    # Subset proteins
    subProts <- perProtNA[perProtNA$PROP_MISSING <= protVal,]$ID
    print(paste0("//// With N missing <= ", protVal, " we get ",length(subProts)," proteins"))
      
    # If rows or columns
    if (idxVal == 1 & length(subProts) > 0) {

      # Enumerate at threshold
      dropProts <- rbind(dropProts,
                         data.frame(PROT_THRESHOLD=protVal, EnumerateNAs(quantTbl[subProts,], idxVal=idxVal)))
      
    } else if (idxVal == 2 & length(subProts) > 0) {
      
      # Enumerate at threshold
      dropProts <- rbind(dropProts,
                         data.frame(PROT_THRESHOLD=protVal, EnumerateNAs(quantTbl[,subProts, drop=FALSE], idxVal=idxVal)))

    }
    
  }
  
  # Make %
  dropProts["PERCENTAGE"] <- paste0(dropProts$PROT_THRESHOLD * 100, "%")
  if (any(dropProts$PERCENTAGE == "100%")) { dropProts[dropProts$PERCENTAGE == "100%",]["PERCENTAGE"] <- "Unfiltered" }
  
  # return
  return(dropProts)
  
}


TestMissingByFactor <- function(metaData, ids, factorCol) {
  
  # Report if Chi or T test (ie factorCol is categorical or numeric)
  if (is.character(metaData[[factorCol]]) | is.factor(metaData[[factorCol]])) { print("///// Running Chi-Squared (Categorical Factor Provided)"); testPerformed <- "Chi-Sq Test" }
  if (is.numeric(metaData[[factorCol]])) { print("///// Running T test (Numerical Covariate Provided)"); testPerformed <- "T Test" }
  
  # Make ids easier to work with
  ids <- gsub(";","_",ids)
  colnames(metaData) <- gsub(";","_",colnames(metaData))
  rownames(metaData) <- gsub(";","_",rownames(metaData))
  
  # Loop over IDs
  resData <- data.frame()
  for (id in ids) {
    
    # If two groups and > 1 per gorup (ie are some missing values)
    if (length(unique(metaData[[id]])) == 2 & 
        length(metaData[[id]][metaData[[id]] == 0]) >= 2 &
        length(metaData[[id]][metaData[[id]] == 1]) >= 2)
      
    {
      
      # Run chi-squared
      if (is.character(metaData[[factorCol]]) | is.factor(metaData[[factorCol]])) { res <- chisq.test(metaData[[factorCol]], metaData[[id]]) }
      if (is.numeric(metaData[[factorCol]])) { metaData[[id]] <- as.character(metaData[[id]]); res <- t.test(as.formula(paste0(factorCol, " ~ ", id)), data=metaData[,c(factorCol,id)]) }
      
      # Combine
      resData <- rbind(resData,
                       data.frame(ID=id, PVAL=res$p.value, STAT=res$statistic, PARAM=res$parameter, METHOD=res$method, TEST=testPerformed))
      
    }
    
  }
  
  # FDR correct
  resData["PADJ"] <- p.adjust(resData$PVAL, method="BH")
  
  # Add ID
  resData <- cbind(Factor=factorCol, resData)
  
  # Return
  return(resData)
  
}


TestFactors <- function(metaData, factorCols, ids) {
  
  # Loop over factors
  outData <- do.call("rbind",lapply(factorCols, function(fctr) {
    
    # Run  test
    print(paste0("##### RUNNING FOR :- ",fctr))
    testBy <- suppressWarnings(TestMissingByFactor(metaData=metaData, ids=ids, factorCol=fctr))
    
    # Sig check
    print(paste0("// There are ",nrow(testBy[testBy$PVAL < 0.05,]), " sig proteins by Nominal P-Value")) # Nominal P Value
    print(paste0("// There are ",nrow(testBy[testBy$PADJ < 0.05,]), " sig proteins by FDR")) # FDR-corrected
    
    # Return
    return(testBy)
    
  }))
  
  # Return
  return(outData)
  
}



GenerateBinaryMissingValTable <- function(mat, meta) {
  
  # Set any non-NA values to 1
  mat_naBinarised <- log2(as.matrix(mat[,meta$Sample_ID]))
  mat_naBinarised[!is.na(mat_naBinarised)] <- 1
  mat_naBinarised[is.na(mat_naBinarised)] <- 0
  mat_naBinarised[is.nan(mat_naBinarised)] <- 0
  
  # Drop any proteins without missing values
  mat_naBinarised <- mat_naBinarised[,colSums(mat_naBinarised) < nrow(mat_naBinarised)]
  
  # Return
  return(mat_naBinarised)
  
}


sim_trKNN_wrapper <- function(data) {
  result <- data %>% as.matrix %>% t %>% imputeKNN(., k=10, distance='truncation', perc=0) %>% t
  return(result)  
}





CalculateKSStats <- function(inputData, chosenSamples, baseLine, scaleData=FALSE) {
  
  # Loop over normalised matrices and perform KS test
  ksCompare <- do.call("rbind",lapply(names(inputData), function(method) 
    
    data.frame(ImputeMethod=method, PairwiseKSTest(inputData[[method]], ids=chosenSamples, scaleData=scaleData))
    
  ))
  
  # Calculate Delta values relative to baseline (non-imputed, non-normalised)
  ksCompare["DELTA_STAT"] <- ksCompare$statistic - ksCompare[ksCompare$ImputeMethod == baseLine,]$statistic
  ksCompare["DELTA_P"] <- ksCompare$p.value - ksCompare[ksCompare$ImputeMethod == baseLine,]$p.value
  
  # Return
  return(ksCompare)
  
}



# Function
PairwiseKSTest <- function(dataTable, ids, scaleData=FALSE, ...) {
  
  # Dependenceis
  require(broom)
  
  # List of performed comparisons
  seenList <- c()
  
  # Make dataframe
  dataTable <- data.frame(dataTable)
  
  # Tidy ID names
  ids <- gsub("[;]",".",ids)
  
  # Loop over QCs
  ksCompare <- do.call("rbind",lapply(ids, function(qcID_1) {
    
    ksSub <- do.call("rbind",lapply(setdiff(ids, qcID_1), function(qcID_2) {
      
      # If comparison has been seen
      compareKey <- paste(sort(c(qcID_1, qcID_2)), collapse="-") 
      if (!compareKey %in% seenList) {
        
        # Scale data
        if (scaleData) { 
          
          v1 <- scale(dataTable[[qcID_1]], center=TRUE, scale=TRUE) 
          v2 <- scale(dataTable[[qcID_2]], center=TRUE, scale=TRUE)  
          
        } else {
          
          v1 <- dataTable[[qcID_1]]
          v2 <- dataTable[[qcID_2]]
          
        }
    
        # Run KS test
        ksRes <- suppressWarnings(ks.test(v1, v2)) #, ...))
        
        # Combine
        outData <- data.frame(SAMPLE1=qcID_1, SAMPLE2=qcID_2, broom::tidy(ksRes))
        
        # Update seen list
        seenList <- c(seenList, compareKey)
        
        # Return
        return(outData)
        
      }
    
    }))
    
    # Return
    return(ksSub)
    
  }))
  
  # Return
  return(ksCompare)
  
}


# Function
PairwiseNRMSE <- function(dataTable, ids, originalData=NULL, ids2=NULL, scaleData=FALSE, ...) {
  
  # Dependenceis
  require(INDperform)
  
  # If second set of ids
  if (is.null(ids2)) { ids2 <- ids }
  
  # List of performed comparisons
  seenList <- c()
  
  # Loop over QCs
  nrmseCompare <- do.call("rbind",lapply(ids, function(qcID_1) {
    
    nrmseSub <- do.call("rbind",lapply(setdiff(ids2, qcID_1), function(qcID_2) {
      
      # If comparison has been seen
      compareKey <- paste(sort(c(qcID_1, qcID_2)), collapse="-") 
      if (!compareKey %in% seenList) {
        
        # Scale data
        if (scaleData) { 
          
          v1 <- scale(dataTable[[qcID_1]], center=TRUE, scale=TRUE) 
          v2 <- scale(dataTable[[qcID_2]], center=TRUE, scale=TRUE)  
          
        } else {
          
          v1 <- dataTable[[qcID_1]]
          v2 <- dataTable[[qcID_2]]
          
        }
        
        # Name as proteins
        names(v1) <- rownames(dataTable)
        names(v2) <- rownames(dataTable)
        
        # Choose whether to only keep features with NAs
        if (!is.null(originalData)) {
          
          # Subset to keep only NA prots
          subsetProts <- rownames(originalData[is.na(originalData[,qcID_1]),qcID_1,drop=FALSE])
        
        # Else just use all proteins
        } else {
          
          # Get all names
          subsetProts <- names(v1)
          
        }
    
        # Combine
        outData <- data.frame(SAMPLE1=qcID_1, SAMPLE2=qcID_2, NRMSE=INDperform::nrmse(v1[subsetProts], v2[subsetProts], ...))
        
        # Update seen list
        seenList <- c(seenList, compareKey)
        
        # Return
        return(outData)
        
      }
      
    }))
    
    # Return
    return(nrmseSub)
  
  }))
  
  # Return
  return(nrmseCompare)
  
}


RunAllImputationMethods <- function(inputData, metaData, groupCol, ...) {
  
  # Order data
  inputData <- inputData[,metaData$Sample_ID]
  
  # Make matrix
  inputData <- as.matrix(inputData)
  
  ##### MCAR/EITHER #####
  
  ## Global similarity methods
  ### SLSA
  print(paste0("///// PERFORMING SLSA"))
  datag_slsa <- as.data.frame(impute.slsa(tab=as.matrix(inputData), conditions=factor(metaData[[groupCol]])))
  
  ### MLE
  print(paste0("///// PERFORMING MLE"))
  datag_mle <- as.data.frame(impute.mle(tab=inputData, conditions=factor(metaData[[groupCol]])))
  
  ### ImpSeq
  datag_impSeq <- as.data.frame(impSeq(inputData))
  
  ### ImpSeq-prob
  print(paste0("///// PERFORMING IMPSEQ"))
  datag_impSeq_prob <- impSeqRob(inputData, alpha=0.9)
  datag_impSeq_prob <- datag_impSeq_prob$x
  datag_impSeq_prob <- as.data.frame(datag_impSeq_prob)
  
  ### BPCA
  ### Note BPCA and SVD expect the samples to be rows, have to transpose
  print(paste0("///// PERFORMING BPCA"))
  datag_bPCA <-pcaMethods::pca(t(as.matrix(inputData)), nPcs = ncol(inputData)-1, method = "bpca", maxSteps =100)
  datag_bPCA <-completeObs(datag_bPCA)
  datag_bPCA <- t(datag_bPCA)
  datag_bPCA <- as.data.frame(datag_bPCA)
  
  ### SVD 
  print(paste0("///// PERFORMING SVD"))
  datag_SVD <-pcaMethods::pca(t(as.matrix(inputData)), nPcs = ncol(inputData)-1, method = "svdImpute")
  datag_SVD <-completeObs(datag_SVD)
  datag_SVD <- as.data.frame(t(datag_SVD))

  
  ## Local similarity methods
  ### Random forest
  print(paste0("///// PERFORMING RF"))
  datag_rf <- as.data.frame(impute.RF(tab=inputData, conditions=factor(metaData[[groupCol]])))
  
  
  ### Mice-cart
  #minum<-5
  #datareadmi <-mice(inputData, m=minum,seed = 1234, method ="cart")
  #newdatareadmi<-0
  
  #for (i in 1:minum) {
  #  newdatareadmi<-complete(datareadmi,action = i)+newdatareadmi
  #}
  
  #datag_mice_cart <-newdatareadmi/minum
  #rownames(datag_mice_cart) <-rownames(inputData)
  
  
  
  ### Mice-norm
  #minum<-5
  #datareadmi <-mice(inputData, m=minum,seed = 1234, method ="norm")
  #newdatareadmi<-0
  
  #for (i in 1:minum) {
  #  newdatareadmi<-complete(datareadmi,action = i)+newdatareadmi
  #}
  #
  #datag_mice_norm <-newdatareadmi/minum
  #rownames(datag_mice_norm) <-rownames(inputData)
  
  ### Seq-KNN
  print(paste0("///// PERFORMING SEQKNN"))
  datag_seqKNN <- tryCatch({as.data.frame(seqKNNimp(inputData, k = 10))},
                           
                            # Error out as NULL
                            error = function(cond) { warning("// Error with method, returning NULL "); NULL })
  
  
  ### trKNN
  print(paste0("///// PERFORMING TRKNN"))
  datag_trKNN <- tryCatch({as.data.frame(t(sim_trKNN_wrapper(t(inputData))))},
      
                           # Error out as NULL
                           error = function(cond) { warning("// Error with method, returning NULL "); NULL })
  
  
  
  ##### MAR ######
  
  ## Single value methods
  ### MinDet
  print(paste0("///// PERFORMING MINDET"))
  data_mat <-as.matrix(inputData)
  datag_mindet <- as.data.frame(imputeLCMD::impute.MinDet(data_mat, q = 0.01))
  
  
  ### MinProb
  print(paste0("///// PERFORMING MINPROB"))
  data_mat <-as.matrix(inputData)
  datag_minprob <- as.data.frame(imputeLCMD::impute.MinProb(data_mat, q = 0.01, tune.sigma = 1))
  
  # Add rownames
  rownames(datag_slsa) <- rownames(inputData)
  rownames(datag_mle) <- rownames(inputData)
  rownames(datag_impSeq_prob) <- rownames(inputData)
  rownames(datag_bPCA) <- rownames(inputData)
  rownames(datag_SVD) <- rownames(inputData)
  rownames(datag_seqKNN) <- rownames(inputData)
  if (!is.null(datag_trKNN)) { rownames(datag_trKNN) <- rownames(inputData) }
  rownames(datag_mindet) <- rownames(inputData)
  rownames(datag_minprob) <- rownames(inputData)
  
  # Combine to list
  imputData <- list(inputData,
                    datag_slsa,
                    datag_mle,
                    datag_impSeq_prob,
                    datag_bPCA,
                    datag_SVD,
                   datag_rf,
                    datag_seqKNN,
                    datag_trKNN,
                    datag_mindet,
                    datag_minprob)
  
  # Tidy up 
  names(imputData) <- c("Un-Imputed", 
                        "SLSA", 
                        "MLE", 
                        "impSeq", 
                        "bPCA", 
                        "SVD", 
                        "RF", 
                        "seqKNN", 
                        "trKNN", 
                        "minDet", 
                        "minProb")
  
  # Return
  return(imputData)
  
}




crossValidateImpute_trainTestSplit <- function(data2Impute, metaData, fileID, groupCol, chosenN_CVs=10, outDir="./", propTrainSplit=0.7, propProts=1, regenFile=TRUE) {
  
  # Initialise outputs
  cvRes <- list()
  
  # Loop over cross-validations
  for (i in 1:chosenN_CVs) {
    
    # If output file doesn't exot
    fileName <- here::here(outDir, paste0(fileID, "_imputationCrossValdiation_CV",i,"_trainProp",propTrainSplit,".rds"))
    if (!file.exists(fileName) | regenFile) {
      
      print(paste0("########## WORKING ON CV = ",i))
      
      set.seed(i)
      
      # Randomise metadata order
      metaData_random <- metaData[sample(1:nrow(metaData), nrow(metaData)),]
      data2Impute <- data2Impute[sample(1:nrow(data2Impute), nrow(data2Impute)),]
      
      # Define train/test amounts
      trainN <- floor(nrow(metaData_random)*propTrainSplit)
      
      # Get train/test metadata
      meta_train <- metaData_random[1:trainN,]
      meta_test  <- metaData_random[(trainN+1):nrow(metaData_random),]
      impTrainData  <- suppressWarnings(RunAllImputationMethods(inputData=data2Impute[1:floor(nrow(data2Impute)*propProts),meta_train$Sample_ID], 
                                                                metaData=meta_train, 
                                                                groupCol=groupCol))
      impTestData   <- suppressWarnings(RunAllImputationMethods(data2Impute[1:floor(nrow(data2Impute)*propProts),meta_test$Sample_ID], meta_test, groupCol))
      
      # Append
      cvRes <- list.append(cvRes, list(TRAIN=impTrainData, TEST=impTestData))
      
      # Save
      saveRDS(list(TRAIN=impTrainData, TEST=impTestData), fileName)
      
      # Else read in file if it already exists
    } else {
      
      # Read in
      cvRes <- list.append(cvRes, readRDS(fileName))
      
    }
    
  }
  
  # Return
  return(cvRes)
  
}




# Loop over CVs
imputeOn_simMissingVals <- function(inData, missingData, metaData, missingMethod, maxPropMissingRows=NULL, maxPropMissingCols=NULL, groupFactor=NULL) {
  
  # Dependencies
  require(missMethods)
  require(INDperform)
  
  # Create output
  simMissM <- matrix(0, nrow(inData[!is.na(rowSums(inData)),]), 0)
  
  # Get max proportion of missing rows
  if (is.null(maxPropMissingRows)) { maxPropMissingRows <- max(EnumerateNAs(inData, idx=1)$PROP_MISSING) }
  print(paste0("/ Dropping proportion of missing rows for whole matrix = ",maxPropMissingRows))
  
  # If group NULL, don't subset
  if (is.null(groupFactor)) {
    
    # Set proportion missing to total matrix
    if (is.null(maxPropMissingCols)) { maxPropMissingCols <- max(EnumerateNAs(inData, idx=2)$PROP_MISSING) }
    print(paste0("/ Dropping proportion of missing columns for whole matrix = ",maxPropMissingRows))
    
    # Drop NAs
    inData_noNA <- inData[!is.na(rowSums(inData)), metaData$Sample_ID]
    
    # CHoose
    switch(missingMethod,
           
           # Missing not at-random
           MCAR={
             
             # Generate missing data
             missData <- delete_MCAR_censoring(inData_noNA, maxPropMissingRows, colnames(inData))
             
           },
           
           # Missing not at-random
           MAR={
             
             # Generate missing data
             missData <- delete_MAR_censoring(inData_noNA, maxPropMissingRows, colnames(inData))
             
           },
           
           # Missing not at-random
           MNAR={
             
             # Generate missing data
             missData <- delete_MNAR_censoring(inData_noNA, maxPropMissingRows, colnames(inData))
             
           })
    
    # Else split 
  } else {
    
    # Loop over groups
    for (group in unique(metaData[[groupFactor]])) {
      
      print(paste0("// Filtering by ",groupFactor," level = ",group))
      
      # Get group IDs
      groupIDs <- metaData[metaData[[groupFactor]] == group,]$Sample_ID
      
      # Set proportion missing to total matrix
      maxPropMissingCols <- max(EnumerateNAs(inData[,groupIDs], idx=2)$PROP_MISSING)
      print(paste0("/ Dropping proportion of missing columns for group = ",maxPropMissingCols))
      
      # Drop NAs
      inData_noNA <- inData[!is.na(rowSums(inData)), metaData$Sample_ID]
      
      # Determine chosen source of missingess
      switch(missingMethod,
             
             # Missing not at-random
             MCAR={
               
               # Generate missing data
               missData <- delete_MAR_censoring(inData_noNA[,groupIDs], maxPropMissingRows, groupIDs)
               
             },
             
             # Missing not at-random
             MAR={
               
               # Generate missing data
               missData <- delete_MAR_censoring(inData_noNA[,groupIDs], maxPropMissingRows, groupIDs)
               
             },
             
             # Missing not at-random
             MNAR={
               
               # Generate missing data
               missData <- delete_MNAR_censoring(inData_noNA[,groupIDs], maxPropMissingRows, groupIDs)
               
             })
      
      # Add imputed data
      simMissM <- cbind(simMissM, missData)
      
    }
    
  }
  
  # Filter by missingness
  rowMissing <- EnumerateNAs(simMissM, 1)
  
  # Ensure filtering to match true data
  simMissM <- simMissM[rowMissing[rowMissing$PROP_MISSING < maxPropMissingRows,]$ID,]
  
  # Test imputation
  simImputM <- suppressWarnings(RunAllImputationMethods(simMissM[,metaData$Sample_ID], metaData, groupFactor))
  
  # Drop unimputed
  simImputM <- simImputM[!grepl("Un-Imputed",names(simImputM))]
  
  # Initialise nrmse output
  nrmseOut <- as.list(rep(NA, length(names(simImputM))))
  names(nrmseOut) <- names(simImputM)
  
  # Retain only proteins with missing data
  simOnlyMissM <- simMissM[is.na(rowSums(simMissM)),]
  
  # Drop NULLs
  simImputM <- simImputM[!sapply(simImputM, is.null)]
  
  # Loop over imputation methods
  print("///// Calculating NRMSE")
  for (impMethod in names(simImputM)) {
    
    # Update imputed data name
    colnames(simImputM[[impMethod]]) <- paste0(colnames(simImputM[[impMethod]]), "_Imputed")
    
    # If not null
    print(paste0("// On :- ",impMethod))
    if (!is.null(simImputM[[impMethod]])) {
      
      # Calculate NRMSE
      nrmse <- PairwiseNRMSE(dataTable=cbind(inData[rownames(simOnlyMissM),], simImputM[[impMethod]][rownames(simOnlyMissM),]),
                             ids=colnames(inData),
                             ids2=colnames(simImputM[[impMethod]]))
      
      # Append
      nrmseOut[[impMethod]] <- data.frame(nrmse, IMP_METHOD=impMethod, MiSSING_METHOD=missingMethod)
      
    }
    
    # Update imputed data name
    colnames(simImputM[[impMethod]]) <- gsub("_Imputed","",colnames(simImputM[[impMethod]]))
    
  }
  
  # Combine to table
  nrmseOut <- do.call("rbind", nrmseOut)
  
  # Return
  return(list(MISSING_DATA=simMissM, IMPUT_DATA=simImputM, NRMSE=nrmseOut))
  
}


calcMeanAbundance_missingAndNonMissingProteins <- function(impData, unImpData, metaData) {
  
  # Generate output list
  impToPlot <- data.frame()
  
  # Impute data
  for (impMethod in names(impData)) { 
    
    # Subset matrix
    unImputeMatrix <- unImpData[rownames(impData[[impMethod]]),metaData$Sample_ID]
    impMatrix <- impData[[impMethod]][,metaData$Sample_ID]
    impMatrix <- impMatrix[,colnames(unImputeMatrix)]
    
    # Get NA values
    naCheck <- reshape2::melt(data.frame(Protein=rownames(unImputeMatrix), unImputeMatrix))
    impTall <- reshape2::melt(data.frame(Protein=rownames(impMatrix), impMatrix))
    
    # Add protein/sample
    impTall["ProtSample"] <- paste(impTall$variable, impTall$Protein, sep="_")
    naCheck["ProtSample"] <- paste(naCheck$variable, naCheck$Protein, sep="_")
    
    # Ensure sorting
    impTall <- impTall[order(impTall$ProtSample, decreasing=TRUE),]
    naCheck <- naCheck[order(naCheck$ProtSample, decreasing=TRUE),]
    
    # Ensure cosnistent order
    if (!all(impTall$Protein == naCheck$Protein)) { warning("Measurements aren't consistently ordered between imputed and unimputed data") }
    
    # Define NAs
    impTall["Missing"] <- "Complete"
    impTall[which(naCheck$Protein %in% unique(naCheck[is.na(naCheck$value),]$Protein)),]["Missing"] <- "Missing"
    
    # Drop QCs
    impTall <- impTall[!grepl("QC",impTall$variable),]
    
    # Aggregate
    aggImpTall <- aggregate(. ~ Missing + Protein, data=impTall[,!grepl("variable|ProtSample",colnames(impTall))], FUN=mean)
    
    # add
    impToPlot <- rbind(impToPlot, data.frame(Method=impMethod, aggImpTall))
    
  }
  
  # Return
  return(impToPlot)
  
}

