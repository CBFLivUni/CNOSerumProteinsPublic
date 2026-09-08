

splitTrainData <- function(meta, propSplit=0.7, seedVal=1, sampleMethod="Random", considerGroup=TRUE, dateCol=NULL, ...) {
  
  # If we dont want to split by group, create a dummy variable with one level
  if (!considerGroup | is.null(considerGroup)) { meta["Condition"] <- "XXXX" }
  
  # Sub-sample by group
  trainData <- do.call("rbind",lapply(unique(meta$Condition), function(cond) {
    
    # Subset to condition
    condMeta <- meta[meta$Condition == cond,]
    
    # If randomly subsampling
    if (sampleMethod == "Random") {
      
      # Subset by train proportion
      set.seed(seedVal)
      subMeta <- condMeta[sample(1:nrow(condMeta), 
                                 ceiling(propSplit*nrow(condMeta))),]
      
      # Else if sampling by date
    } else if (sampleMethod == "Date") {
      
      # Order by date
      condMeta <- condMeta[order(condMeta[[dateCol]], decreasing=FALSE),]
      
      # Subsample
      subMeta <- condMeta[seq_len(ceiling(propSplit*nrow(condMeta))),]
      
    }
    
    # Return
    return(subMeta)
    
  }))
  
  # Return
  return(trainData)
  
}



splitOverIters <- function(meta, qntTable, classVar, splitProp=0.8, method="Random", scaleData=TRUE, nIters=10, underSampleN=NULL, ...) {
  
  # If not undersampling
  if (is.null(underSampleN)) {
    
    # Loop 10 test/train splits
    dataSplits <- lapply(seq_len(nIters), function(i) {
      
      # Get train IDs
      trainIDs <- splitTrainData(meta=meta, propSplit=splitProp, seedVal=i, sampleMethod="Random", ...)$Sample_ID
      testIDs <- meta[!meta$Sample_ID %in% trainIDs,]$Sample_ID
      
      # Get train/test IDs
      factorTrain <- meta[trainIDs,][[classVar]]
      names(factorTrain) <- trainIDs
      factorTest <- meta[testIDs,][[classVar]]
      names(factorTest) <- testIDs
      
      trainRes <- t(as.matrix(qntTable[,trainIDs]))
      testRes <- t(as.matrix(qntTable[,testIDs]))
      
      # If wanting to scale
      if (scaleData) {
        
        # Scale data
        trainResScale <- apply(trainRes, 2, scale, scale=TRUE, center=TRUE)
        testResScale <- apply(testRes, 2, scale, scale=TRUE, center=TRUE)
        rownames(trainResScale) <- rownames(trainRes)
        rownames(testResScale) <- rownames(testRes)
        
        # Assemble output split - scaled
        out <- list("train" = trainResScale, 
                    "test"  = testResScale,
                    "factor_train" = factorTrain,
                    "factor_test"  = factorTest,
                    "idx" = which(!meta$Sample_ID %in% testIDs))
      
      # Else don't scale
      } else {
        
        # Assemble output split
        out <- list("train" = trainRes, 
                    "test"  = testRes,
                    "factor_train" = factorTrain,
                    "factor_test"  = factorTest,
                    "idx" = which(!meta$Sample_ID %in% testIDs))
        
      }
      
    })
    
  # Else perform nuder-sampling
  } else {
    
    # Get output
    dataSplit_underSamples <- vector(mode="list",
                                     length=underSampleN)
    
    # Loop over undersamples
    for (j in seq_len(underSampleN)) {
      
      # Loop 10 test/train splits
      dataSplits <- lapply(seq_len(nIters), function(i) {
        
        # Get train IDs
        trainIDs <- splitTrainData(meta=meta, propSplit=splitProp, seedVal=i, sampleMethod="Random", ...)$Sample_ID
        testIDs <- meta[!meta$Sample_ID %in% trainIDs,]$Sample_ID
        
        # Get train/test IDs
        factorTrain <- meta[trainIDs,][[classVar]]
        names(factorTrain) <- trainIDs
        factorTest <- meta[testIDs,][[classVar]]
        names(factorTest) <- testIDs
        
        trainRes <- t(as.matrix(qntTable[,trainIDs]))
        testRes <- t(as.matrix(qntTable[,testIDs]))
        
        # If wanting to scale
        if (scaleData) {
          
          # Scale data
          trainResScale <- apply(trainRes, 2, scale, scale=TRUE, center=TRUE)
          testResScale <- apply(testRes, 2, scale, scale=TRUE, center=TRUE)
          rownames(trainResScale) <- rownames(trainRes)
          rownames(testResScale) <- rownames(testRes)
          
          # Assemble output split - scaled
          out <- list("train" = trainResScale, 
                      "test"  = testResScale,
                      "factor_train" = factorTrain,
                      "factor_test"  = factorTest,
                      "idx" = which(!meta$Sample_ID %in% testIDs))
          
          # Else don't scale
        } else {
          
          # Assemble output split
          out <- list("train" = trainRes, 
                      "test"  = testRes,
                      "factor_train" = factorTrain,
                      "factor_test"  = factorTest,
                      "idx" = which(!meta$Sample_ID %in% testIDs))
          
        }
        
      })
      
      # Add under samples
      dataSplit_underSamples[[j]]
      
    }
    
    
  }
  
  # Return
  return(dataSplits)
  
}




splitLOO <- function(meta, 
                     qntTable, 
                     classVar, 
                     covs2Select=NULL,
                     downSampleNReps=NULL,
                     underSampleToLevel="CNO") {
  
  # If dont want to undersample
  if (is.null(downSampleNReps)) {
    
    # Loop 10 test/train splits
    dataSplits <- lapply(seq_len(nrow(meta)), function(i) {
      
      set.seed(i)
      
      # Get train IDs
      trainIDs <- meta[-i,]$Sample_ID
      testIDs <- meta[i,]$Sample_ID
      
      # Get Validation IDs
      factorTrain <- meta[trainIDs,][[classVar]]
      names(factorTrain) <- trainIDs
      factorTest <- meta[testIDs,,drop=FALSE][[classVar]]
      names(factorTest) <- testIDs
      
      # Get matrix
      trainRes <- t(as.matrix(qntTable[,trainIDs]))
      testRes <- t(as.matrix(qntTable[,testIDs,drop=FALSE]))
      
      # If any covariates to select
      if (!is.null(covs2Select)) {
        
        # Add covariates
        trainRes <- cbind(trainRes, meta[trainIDs,covs2Select, drop=FALSE])
        testRes <- cbind(testRes, meta[testIDs,covs2Select, drop=FALSE])
        
      }
      
      # Assemble output split
      out <- list("train" = trainRes, 
                  "test"  = testRes,
                  "factor_train" = factorTrain,
                  "factor_test"  = factorTest,
                  "idx" = which(!meta$Sample_ID %in% testIDs))
      
    })
    
    # Return
    return(dataSplits)
    
    # Else under-sample
  } else {
    
    # Get output
    dataSplit_underSamples <- vector(mode="list",
                                     length=downSampleNReps)
    
    # Get IDs of reference group to downsample to
    downSampleGroupIDs <- meta[meta[[classVar]] == "CNO",]$Sample_ID
    
    # Get number to undersample
    nSamples2Sample <- length(downSampleGroupIDs)
    
    # Loop over undersamples
    for (j in seq_len(downSampleNReps)) {
      
      # Re-sample metadata
      set.seed(j)
      meta <- meta[sample(seq_len(nrow(meta))),]
      
      # Sample a list of samples, that excludes the left-out sample
      sampledIDs <- lapply(seq_len(nrow(meta)), function(j) { 
        
        # Sample test set
        set.seed(j) 
        sample(meta[meta[[classVar]] != "CNO" & 
                      meta$Sample_ID != meta[j,]$Sample_ID,]$Sample_ID,
               nSamples2Sample) 
        
      })
      
      # Loop over samples
      dataSplits <- lapply(seq_len(nrow(meta)), function(i) {
        
        # Set seed
        set.seed(i)
        
        # Get test ID
        testIDs <- meta[i,]$Sample_ID # Test sample is drawn from the NON-down-sampled data
        
        # Get if test sample is in the reference group
        testIDInGroup <- ifelse(testIDs %in% downSampleGroupIDs, 1, 0)
        
        # test ID is down-sampled set
        if (testIDInGroup == 1) {
          
          # Get idx to drop
          dropIdx <- sample(seq_len(length(sampledIDs[[i]])), testIDInGroup)
          
          # Subset meta
          subMeta <- meta[c(sampledIDs[[i]][-dropIdx], # Randomly drop 1 to match reference lost below 
                            setdiff(downSampleGroupIDs, testIDs)),] # Need to drop cases where test is reference data
          
          # Else just keep set
        } else {
          
          # Subset meta
          subMeta <- meta[c(sampledIDs[[i]], downSampleGroupIDs),]
          
        }
        
        # Get train IDs
        trainIDs <- c(subMeta$Sample_ID) # Training set is drawn from the down-sampled data
        
        # Get Validation IDs
        factorTrain <- meta[trainIDs,][[classVar]]
        names(factorTrain) <- trainIDs
        factorTest <- meta[testIDs,,drop=FALSE][[classVar]]
        names(factorTest) <- testIDs
        
        # Get matrix
        trainRes <- t(as.matrix(qntTable[,trainIDs]))
        testRes <- t(as.matrix(qntTable[,testIDs,drop=FALSE]))
        
        # If any covariates to select
        if (!is.null(covs2Select)) {
          
          # Add covariates
          trainRes <- cbind(trainRes, subMeta[trainIDs,covs2Select, drop=FALSE])
          testRes <- cbind(testRes, subMeta[testIDs,covs2Select, drop=FALSE])
          
        }
        
        # Assemble output split
        out <- list("train" = trainRes, 
                    "test"  = testRes,
                    "factor_train" = factorTrain,
                    "factor_test"  = factorTest,
                    "idx" = which(!meta$Sample_ID %in% testIDs))
        
      })
      
      # Add to final output
      dataSplit_underSamples[[j]] <- dataSplits
      
    }
    
    # Return
    return(dataSplit_underSamples)
    
  }
  
}


cleanVars <- function(metricData) {
  
  # Tidy conditions
  metricData$Group <- gsub("Class: ", "", metricData$Group)
  
  # Fix metric names
  colnames(metricData) <- mgsub(colnames(metricData), 
                                c("Pos.Pred.Value", 
                                  "Neg.Pred.Value",
                                  "Detection",
                                  "Prevalence",
                                  "Balanced.Accuracy"),
                                c("PPV", 
                                  "NPV",
                                  "Det",
                                  "Prev",
                                  "Bal Acc"))
  colnames(metricData) <- gsub("[.]", " ", colnames(metricData))
  
  # Remove redundant recall & Precition
  metricData <- metricData[,!colnames(metricData) %in% c("Precision", "Recall", "Prev", "Det Prev", "Det Rate"), drop=FALSE]
  
  # Return
  return(metricData)
  
}


addCovs_toSplitData <- function(splitList, meta, covs2Select) {
  
  # Loop over split lists
  for (i in seq_len(length(splitList))) {
    
    # Add covariates
    splitList[[i]][["train"]] <- data.frame(cbind(splitList[[i]]$train, 
                                                  meta[rownames(splitList[[i]]$train),covs2Select]))
    splitList[[i]][["test"]] <- data.frame(cbind(splitList[[i]]$test,
                                                 meta[rownames(splitList[[i]]$test),covs2Select]))
    
  }
  
  # Return
  return(splitList)
  
}



# Function to find inconsistent elements between two named lists
getInconsistentClasses <- function(list1, list2) {
  
  # Check if the names are identical
  if (!all(names(list1) == names(list2))) {
    stop("The names (IDs) in the lists do not match!")
  }
  
  # If mismatched
  if (any(list1 != list2)) {
    
    # Compare the classes in both lists and identify mismatches
    mismatches <- list1 != list2
    
    # Extract the names and classes where mismatches occur
    inconsistent_names <- names(list1)[mismatches]
    inconsistent_list1 <- list1[inconsistent_names]
    inconsistent_list2 <- list2[inconsistent_names]
    
    # Combine and return the result as a data frame
    result <- data.frame(
      ID = inconsistent_names,
      Predicted = inconsistent_list1,
      Actual = inconsistent_list2,
      stringsAsFactors = FALSE
    )
    
  # Else if not mismatched
  } else {
    
    # Combine and return the result as a data frame
    result <- data.frame(
      ID = names(list1),
      Predicted = list1,
      Actual = list2,
      stringsAsFactors = FALSE
    )
    
  }
  
  return(result)
  
}


testModel <- function(qnt, 
                      trainData, 
                      testData, 
                      metaData,
                      responseVar, 
                      selectProts, 
                      covs2Include=NULL,
                      ...) {
  
  # Dependencies
  require(ranger)
  allLevels <- unique(c(as.factor(testData[[responseVar]]), 
                        as.factor(trainData[[responseVar]]))) # Need this as sometimes levels missing in either
  
  # Sort wonky protein names
  rownames(qnt) <- gsub("[;]","_",rownames(qnt))
  selectProts <- gsub("[;]","_",selectProts)
  
  # Fix data 
  trainData <- cbind(trainData[,c(responseVar, covs2Include), drop=FALSE], 
                     (t(qnt[,rownames(trainData)])))
  trainData <- cbind(trainData[,c(responseVar,covs2Include), drop=FALSE],
                     data.frame(apply(trainData[,!colnames(trainData) %in% c(responseVar,covs2Include)], 2, as.numeric)))
  testData <- cbind(testData[,c(responseVar, covs2Include), drop=FALSE], 
                    data.frame(t(qnt[,rownames(testData)])))
  testData <- cbind(testData[,c(responseVar,covs2Include), drop=FALSE],
                    data.frame(apply(testData[,!colnames(testData) %in% c(responseVar,covs2Include), drop=FALSE], 2, as.numeric)))
  # testData <- cbind(testData[,c(responseVar,covs2Include)],
  #                   data.frame(apply(t(testData[,!colnames(testData) %in% c(responseVar,covs2Include), drop=FALSE], 2, as.numeric))))
  
  # Fix colnames
  colnames(trainData)[1:length(c(responseVar,covs2Include))] <- c(responseVar,covs2Include)
  colnames(testData)[1:length(c(responseVar,covs2Include))] <- c(responseVar,covs2Include)
  
  # Build final model
  trainData[[responseVar]] <- factor(trainData[[responseVar]],
                                     levels=allLevels)
  finalModel <- ranger(data = trainData,
                       formula = as.formula(paste0(responseVar, " ~ ", paste(c(selectProts, covs2Include), collapse=" + "))),
                       always.split.variable=covs2Include,
                       respect.unordered.factors = "order",
                       ...)
  
  # Model prediction
  finalPred <- predict(object = finalModel, 
                       data = testData)
  
  # Sort levels
  predicted <- factor(finalPred$predictions, 
                      levels=allLevels)
  actual <- factor(testData[[responseVar]], 
                   levels=allLevels)
  
  # If just one sample
  if (length(actual) > 1) {
    
    # Calculate MCC
    allMetrics <- calcOneVsAllMetrics(model=list(list(rfmodelRes=list(Predictions=finalPred))),
                                      splitData=list(list(factor_test=setNames(testData[[responseVar]],
                                                                               rownames(testData)))),
                                      metaData=testData,
                                      multiGroup=ifelse(length(unique(testData[[responseVar]])) > 2, TRUE, FALSE),
                                      predCol=responseVar)
    
  # Else just NA
  } else { allMetrics <- NA }
    
  # Return
  return(list(Model=finalModel,
              rfmodelRes=list(Predictions=finalPred),
              InputProts=selectProts,
              additionalCovs=covs2Include,
              Metrics=allMetrics,
              trainData=trainData,
              testData=testData))
  
}


testModel_validation <- function(qnt, 
                                 trainData, 
                                 testData, 
                                 responseVar, 
                                 selectProts, 
                                 covs2Include=NULL,
                                 ...) {
  
  # Dependencies
  require(ranger)
  allLevels <- unique(c(as.factor(testData[[responseVar]]), 
                        as.factor(trainData[[responseVar]]))) # Need this as sometimes levels missing in either
  
  # Sort wonky protein names
  rownames(qnt) <- gsub("[;]","_",rownames(qnt))
  colnames(qnt) <- gsub("[;]","_",colnames(qnt))
  selectProts <- gsub("[;]","_",selectProts)
  
  # Fix data 
  # trainData <- cbind(trainData[,c(responseVar, covs2Include), drop=FALSE], 
  #                    (t(qnt[,trainData$Sample_ID])))
  # trainData <- cbind(trainData[,c(responseVar,covs2Include)],
  #                    data.frame(apply(trainData[,!colnames(trainData) %in% c(responseVar,covs2Include)], 2, as.numeric)))
  # testData <- cbind(testData[,c(responseVar, covs2Include), drop=FALSE], 
  #                   data.frame(t(qnt[,testData$Sample_ID])))
  # testData <- cbind(testData[,c(responseVar,covs2Include)],
  #                   data.frame(t(apply(testData[,!colnames(testData) %in% c(responseVar,covs2Include), drop=FALSE], 2, as.numeric))))
  # testData <- cbind(testData[,c(responseVar,covs2Include)],
  
  # Fix colnames
  colnames(trainData)[1:length(c(responseVar,covs2Include))] <- c(responseVar,covs2Include)
  colnames(testData)[1:length(c(responseVar,covs2Include))] <- c(responseVar,covs2Include)
  
  # Build final model
  trainData[[responseVar]] <- factor(trainData[[responseVar]],
                                     levels=allLevels)
  finalModel <- ranger(data = trainData,
                       formula = as.formula(paste0(responseVar, " ~ ", paste(c(selectProts, covs2Include), collapse=" + "))),
                       always.split.variable=covs2Include,
                       respect.unordered.factors = "order",
                       ...)
  
  # Model prediction
  finalPred <- predict(object = finalModel, 
                       data = testData)
  
  # Sort levels
  predicted <- factor(finalPred$predictions, 
                      levels=allLevels)
  actual <- factor(testData[[responseVar]], 
                   levels=allLevels)
  if (!is.null(names(finalPred$Prediction$predictions))) { names(finalPred$Prediction$predictions) <- rownames(testData) }
  if (!is.null(names(finalPred$predictions))) { names(finalPred$predictions) <- rownames(testData) }
  
  # Get confusion matrix
  ConfMat <- caret::confusionMatrix(predicted,
                                    actual)
  
  # If just one sample
  if (length(actual) > 1) {
    
    # Calculate MCC
    rfMCC <- mcc(predicted,
                 actual)
    
    # Else just NA
  } else { rfMCC <- NA }
  
  # if (rownames(testData) == "C456") { print(finalModel) }
  
  # if (rownames(testData) == "C456") { print(paste0(testData[[responseVar]], " // ", predicted)) }
  
  # Return
  return(list(rfmodelRes=finalModel,
              Prediction=finalPred,
              InputProts=selectProts,
              additionalCovs=covs2Include,
              MCC=rfMCC,
              trainData=trainData,
              testData=testData))
  
}



testDecoyModel <- function(qnt, 
                           trainData,
                           testData,
                           nProts, 
                           nReps, 
                           covs2Include=NULL, 
                           ...) {
  
  # Dependencies
  require(rlist)
  
  # Fix names
  rownames(qnt) <- gsub("[;]","_",rownames(qnt))
  
  # Define output
  outList <- list()
  
  # Loop over splits
  for (i in seq_len(nReps)) {
    
    # Subsample proteins
    set.seed(i)
    subProts <- sample(rownames(qnt[,!colnames(qnt) %in% covs2Include]), 
                       nProts)
    
    # Define train/test factor splits
    validRes <- testModel(qnt=qnt, 
                          trainData=trainData,
                          testData=testData,
                          selectProts=subProts, 
                          ...)
    
    # Append
    outList <- rlist::list.append(outList, 
                                  validRes)
    
  }
  
  # Return
  names(outList) <- seq_len(nReps)
  return(outList)
  
}


testDecoyProteins_onValidData <- function(allProts,
                                          nProts,
                                          dataSplits, 
                                          metaData,
                                          responseVar,
                                          nReps=10,
                                          ...) {
  
  # Generate output
  outList <- list()
  metricList <- list()
  
  # Loop over reps
  for (i in seq_len(nReps)) {
    
    print(i)
    
    # Set seed
    set.seed(i)
    
    # Sample N prots
    sampledProts <- sample(allProts, 
                           nProts)
    sampledProts <- gsub("[;]|[.]", "_", sampledProts)
    
    # Initialise output
    repList <- list()
    
    # Loop over splits
    for (j in seq_along(dataSplits)) {
      
      # Copy splits
      splitCopy <- dataSplits
      
      # Combine all samples
      qnt <- t(rbind(splitCopy[[j]]$train,
                     splitCopy[[j]]$test))
      
      # Sort wonky protein names
      colnames(splitCopy[[j]]$train) <- gsub("[;]|[.]","_",colnames(splitCopy[[j]]$train))
      colnames(splitCopy[[j]]$test) <- gsub("[;]|[.]","_",colnames(splitCopy[[j]]$test))
      rownames(qnt) <- gsub("[;]|[.]","_",rownames(qnt))
      
      # Filter proteins in test/train
      splitCopy[[j]]$train <- splitCopy[[j]]$train[,sampledProts]
      splitCopy[[j]]$test <- splitCopy[[j]]$test[,sampledProts,drop=FALSE]
      
      # Get train/test data
      trainData <- cbind(X=splitCopy[[j]]$factor_train,
                         data.frame(splitCopy[[j]]$train))
      testData <- cbind(X=splitCopy[[j]]$factor_test,
                        data.frame(splitCopy[[j]]$test))
      rownames(trainData) <- names(splitCopy[[j]]$factor_train)
      rownames(testData) <- names(splitCopy[[j]]$factor_test)
      colnames(trainData)[[1]] <- responseVar
      colnames(testData)[[1]] <- responseVar
      
      # Test on validation data
      validRes <- testModel_validation(qnt=qnt[sampledProts,],
                                       selectProts=sampledProts,
                                       trainData=trainData,
                                       testData=testData,
                                       responseVar=responseVar) #,
      #...)
      
      # Add real labels
      validRes <- append(validRes,
                         list(RandomProteins=sampledProts))
      
      # Append
      repList <- rlist::list.append(repList, 
                                    validRes)
      
    }
    
    # Calculate metrics
    allMetrics <- calcOneVsAllMetrics(model=list(list(rfmodelRes=list(Predictions=list(predictions=sapply(seq_along(repList), function(i) repList[[i]]$Prediction$predictions))))),
                                      splitData=list(list(factor_test=unlist(lapply(seq_along(dataSplits), function(i) dataSplits[[i]]$factor_test)))),
                                      metaData=metaData,
                                      multiGroup=ifelse(length(unique(metaData[[responseVar]])) > 2, TRUE, FALSE),
                                      predCol=responseVar)
    # allMetrics <- calcOneVsAllMetrics(model=sapply(seq_along(repList), function(i) list(list(rfmodelRes=list(Predictions=repList[[i]]$rfmodelRes)))),
    #                                   splitData=list(list(factor_test=sapply(seq_along(repList), function(i) setNames(repList[[i]]$testData[[responseVar]], rownames(repList[[i]]$testData))))),
    #                                   metaData=metaData,
    #                                   multiGroup=ifelse(length(unique(metaData[[responseVar]])) > 2, TRUE, FALSE),
    #                                   predCol=responseVar)
    
    # Add to final output
    outList <- rlist::list.append(outList,
                                  repList)
    
    # Add metrics
    metricList <- rlist::list.append(metricList,
                                     allMetrics)
    
  }
  
  # Set names
  return(list(Models=outList,
              Metrics=metricList))
  
}


getPredictsPerShuff <- function(shuffData) {
  
  # Generate output
  outData <- data.frame()
  
  # Loop over reps
  for (i in seq_along(shuffData$Metrics)) {
    
    # Combined
    outData <- rbind(outData,
                     cbind(Rep=i,
                           shuffData$Metrics[[i]]))
    
  }
  
  # Return
  return(outData)
  
}




shuffleValidData <- function(dataSplits, 
                             metaData,
                             selectProts,
                             responseVar,
                             nReps=10,
                             ...) {
  
  # Generate output
  outList <- list()
  metricList <- list()
  
  # Loop over reps
  for (i in seq_len(nReps)) {
    
    print(i)
    
    # Set seed
    set.seed(i)
    
    # Get actual
    realLabels <- metaData[[responseVar]]
    names(realLabels) <- metaData$Sample_ID
    
    # Sample labels
    shuffLabels <- sample(metaData[[responseVar]])
    names(shuffLabels) <- metaData$Sample_ID
    
    # Get rep-specific list
    repList <- list()
    
    # Loop over splits
    for (j in seq_along(dataSplits)) {
      
      # Copy splits
      splitCopy <- dataSplits
      
      # Combine all samples
      qnt <- t(rbind(splitCopy[[j]]$train,
                     splitCopy[[j]]$test))
      
      # test/train sample IDs
      trainIDs <- names(splitCopy[[j]]$factor_train)
      testIDs <- names(splitCopy[[j]]$factor_test)
      
      # Update train/test labels
      splitCopy[[j]]$factor_train <- metaData[trainIDs,
                                              responseVar]
      splitCopy[[j]]$factor_test <- shuffLabels[testIDs]
      names(splitCopy[[j]]$factor_train) <- trainIDs
      names(splitCopy[[j]]$factor_test) <- testIDs
      
      # Get train/test data
      trainData <- cbind(X=splitCopy[[j]]$factor_train,
                         data.frame(splitCopy[[j]]$train))
      testData <- cbind(X=splitCopy[[j]]$factor_test,
                        data.frame(splitCopy[[j]]$test))
      colnames(trainData)[[1]] <- responseVar
      colnames(testData)[[1]] <- responseVar
      
      # Test on validation data
      validRes <- testModel_validation(qnt=qnt,
                                       selectProts=selectProts,
                                       trainData=trainData,
                                       testData=testData,
                                       responseVar=responseVar,
                                       # mtry=bestModel_biModal$OptimalParams$mTry,
                                       # num.trees=bestModel_biModal$OptimalParams$nRFTrees,
                                       # splitrule=bestModel_biModal$OptimalParams$SplitRule,
                                       # min.node.size=bestModel_biModal$OptimalParams$MinNodeSize,
                                       # max.depth=bestModel_biModal$OptimalParams$MaxDepth,
                                       ...)
      
      # Add real labels
      validRes <- append(validRes,
                         list(RealLabels=realLabels[testIDs]))
      
      # Append
      repList <- rlist::list.append(repList, 
                                    validRes)
      
    }
    
    # print("/////////////")
    # print(as.character(unname(sapply(seq_along(repList), function(i) repList[[i]]$Prediction$predictions))))
    # print(as.character(unname(sapply(seq_along(repList), function(i) repList[[i]]$RealLabels))))
    
    # Calculate metrics
    allMetrics <- calcOneVsAllMetrics(model=list(list(rfmodelRes=list(Predictions=list(predictions=sapply(seq_along(repList), function(i) repList[[i]]$Prediction$predictions))))),
                                      splitData=list(list(factor_test=shuffLabels)),
                                      metaData=metaData,
                                      multiGroup=ifelse(length(unique(realLabels)) > 2, TRUE, FALSE),
                                      predCol=responseVar)
    
    # Add to final output
    outList <- rlist::list.append(outList,
                                  repList)
    
    # Add metrics
    metricList <- rlist::list.append(metricList,
                                     allMetrics)
    
  }
  
  # Set names
  return(list(Models=outList,
              Metrics=metricList))
  
}



testShuffleModel <- function(qnt, 
                             trainData,
                             testData, 
                             selectProts, 
                             responseVar, 
                             nReps,
                             ...) {
  
  # Dependencies
  require(rlist)
  
  # Fix names
  rownames(qnt) <- gsub("[;]","_",rownames(qnt))
  selectProts <- gsub("[;]","_",selectProts)
  
  # Define output
  outList <- list()
  
  # Loop over splits
  for (i in seq_len(nReps)) {
    
    # Shuffle labels
    set.seed(i*7) # Multiplying this  as seeds 5 & 8 gave reorders of split 9 that were the same as before
    # if (i == 2) { set.seed(111) } # Seed 2 is the same as one of the splits, so change here
    
    # Get old labels
    oldLabels <- testData[[responseVar]]
    
    # If just LOOCV, ie one test sample
    if (nrow(testData) == 1) { 
      
      # get all possible labels
      possibleLabels <- sort(union(as.character(testData[[responseVar]]), 
                                   as.character(unique(trainData[[responseVar]]))))
      
      # Just sample from all possible classes
      newLabel <- sample(possibleLabels, 1)
      testData[[responseVar]] <- newLabel
      
      # Fix cases where 
      # while (all(oldLabels) == newLabel)
    
    # Else shuffle multiple samples
    } else {
      
      # Shuffle order
      testData[[responseVar]] <- sample(testData[[responseVar]])
      
    }
    
    # Test on shuffled data
    validRes <- testModel(qnt=qnt,
                          selectProts=selectProts,
                          testData=testData,
                          trainData=trainData,
                          responseVar=responseVar,
                          ...)
    
    # Append
    outList <- rlist::list.append(outList, 
                                  validRes)
    
  }
  
  # Return
  names(outList) <- seq_len(nReps)
  return(outList)
  
}


testModel_onDataSplits <- function(dataSplits, metaData, selectProts, ...) {
  
  # Get output
  outList <- list()  
  
  # Loop over splits
  for (i in seq_len(length(dataSplits))) {
    
    # Test on validation data
    res <- testModel(qnt=t(rbind(dataSplits[[i]]$train,
                                            dataSplits[[i]]$test)),
                                trainData=metaData[rownames(dataSplits[[i]]$train),],
                                testData=metaData[rownames(dataSplits[[i]]$test),],
                                selectProts=selectProts,
                                ...)
    
    # Add
    outList <- rlist::list.append(outList,
                                  res)
    
  }
  
  # Rename
  names(outList) <- seq_len(length(dataSplits))
  
  # Return
  return(outList)
  
}


testDecoys_onSplits <- function(dataSplits, 
                                metaData, 
                                selectProts, 
                                covs2Include=NULL,
                                ...) {
  
  # Generate output
  outList <- list()
  
  # Loop over splits
  for (i in seq_len(length(dataSplits))) {
    
    # Flip data back
    dataSplits[[i]]$train <- t(dataSplits[[i]]$train)
    dataSplits[[i]]$test <- t(dataSplits[[i]]$test)
    
    # Drop selected proteins from split
    qnt <- cbind(dataSplits[[i]]$train[!rownames(dataSplits[[i]]$train) %in% selectProts,],
                 dataSplits[[i]]$test[!rownames(dataSplits[[i]]$test) %in% selectProts,,drop=FALSE])
    
    # Perform decoy on biModal model
    decoyRes <- testDecoyModel(qnt=qnt,
                               trainData=metaData[names(dataSplits[[i]]$factor_train),],
                               testData=metaData[names(dataSplits[[i]]$factor_test),],
                               covs2Include=covs2Include,
                               ...)
    
    # Append
    outList <- rlist::list.append(outList, 
                                  decoyRes)
    
  }
  
  # Set names
  names(outList) <- seq_len(length(dataSplits))
  
  # Return
  return(outList)
  
}



shuffleDiscoverySplits <- function(dataSplits,
                                   metaData, 
                                   selectProts, 
                                   ...) {
  
  # Generate output
  outList <- list()
  
  # Loop over splits
  for (i in seq_len(length(dataSplits))) {
    
    set.seed(i)
    
    # Combine all samples
    qnt <- t(rbind(dataSplits[[i]]$train,
                   dataSplits[[i]]$test))
    
    # Perform decoy on multiGroup model
    shuffRes <- testShuffleModel(qnt=qnt,
                                 selectProts=selectProts,
                                 trainData=metaData[names(dataSplits[[i]]$factor_train),],
                                 testData=metaData[names(dataSplits[[i]]$factor_test),],
                                 metaData=metaData,
                                 ...)
    
    # Append
    outList <- rlist::list.append(outList, 
                                  shuffRes)
    
  }
  
  # Set names
  names(outList) <- seq_len(length(dataSplits))
  return(outList)
  
}


decoyDiscoverySplits <- function(dataSplits,
                                 metaData, 
                                 nProts,
                                 ...) {
  
  # Generate output
  outList <- list()
  
  # Loop over splits
  for (i in seq_len(length(dataSplits))) {
    
    set.seed(i)
    
    # Combine all samples
    qnt <- t(rbind(dataSplits[[i]]$train,
                   dataSplits[[i]]$test))
    
    # Perform decoy on multiGroup model
    shuffRes <- testDecoyModel(qnt=qnt,
                               nProts=nProts,
                               trainData=metaData[names(dataSplits[[i]]$factor_train),],
                               testData=metaData[names(dataSplits[[i]]$factor_test),],
                               ...)
    
    # Append
    outList <- rlist::list.append(outList, 
                                  shuffRes)
    
  }
  
  # Set names
  names(outList) <- seq_len(length(dataSplits))
  return(outList)
  
}


getRankedModels <- function(modelData,
                            metricName="Model") {
  
  # Order by variance
  aggVar <- aggregate(as.formula(paste0(metricName, " ~ Model")),
                      data=modelData[,c("Model", metricName)],
                      FUN=var)
  colnames(aggVar)[[2]] <- "Var"
  aggVar <- aggVar[order(aggVar$Var, decreasing=FALSE),]
  
  # Aggregate to mean
  rankedModels <- aggregate(as.formula(paste0(metricName, " ~ Model")), 
                       data=modelData[,c("Model", metricName)], 
                       FUN=mean) # mean
  # Combine
  rankedModels <- merge(rankedModels, aggVar, by="Model")
  
  # Order by variance then mean
  rankedModels <- rankedModels[order(rankedModels$Var, decreasing=FALSE),]
  rankedModels <- rankedModels[order(rankedModels[[metricName]], decreasing=ifelse(metricName == "BrierScore", FALSE, TRUE)),]
  
  # Add model data
  rankedModels <- merge(rankedModels, 
                        modelData[!duplicated(modelData$Model), c("Model","Condition","nRFTrees","MaxDepth","mTry","SplitRule","MinNodeSize","SampleFrac")],
                        by="Model")
  
  # Return
  return(rankedModels)
  
}


calcOneVsAllMCC <- function(predData,
                            allGroups) {
  
  # Compute metrics for each pair
  mccs <- sapply(unique(predData$Actual), function(cond) {
    
    # Copy data
    pairData <- predData
    pairData$Actual <- as.character(pairData$Actual)
    pairData$Predicted <- as.character(pairData$Predicted)
    
    # # Replace other groups with "Other"
    if (sum(pairData$Predicted != cond) > 1 ) { pairData[pairData$Predicted != cond,]["Predicted"] <- "Other" }
    if (sum(pairData$Actual != cond) > 1 ) { pairData[pairData$Actual != cond,]["Actual"] <- "Other" }
    
    # Calculate MCC
    mccVal <- mltools::mcc(pairData$Predicted,
                           pairData$Actual)
    
    # Return
    return(mccVal)
    
  })
  names(mccs) <- unique(predData$Actual)
  
  # If missing
  presentGroups <- names(mccs)
  missingGroups <- setdiff(allGroups,
                           presentGroups)
  mccs <- c(mccs, rep(NA, length(missingGroups)))
  names(mccs) <- c(presentGroups, missingGroups)
  
  # Sort
  sortedNames <- sort(names(mccs))
  
  # Return
  return(mccs[sortedNames])
  
}


calcOneVsAllAUROC <- function(predData,
                              allGroups) {
  
  # Compute metrics for each pair
  conds2Test <- setdiff(unique(predData$Actual), "Other")
  aucrocs <- sapply(conds2Test, function(cond) {
    
    # Copy data
    pairData <- predData
    pairData$Actual <- as.character(pairData$Actual)
    pairData$Predicted <- as.character(pairData$Predicted)
    
    # # Replace other groups with "Other"
    if (sum(pairData$Predicted != cond) > 1 ) { pairData[pairData$Predicted != cond,]["Predicted"] <- "Other" }
    if (sum(pairData$Actual != cond) > 1 ) { pairData[pairData$Actual != cond,]["Actual"] <- "Other" }
    
    # Calculate AUC ROC
    rocVal <- suppressMessages(roc(as.numeric(factor(pairData$Actual, levels=c("Other", as.character(cond)))) - 1,
                                   as.numeric(factor(pairData$Predicted, levels=c("Other", as.character(cond)))) - 1))
    aucVal <- pROC::auc(rocVal)
    
    # Return
    return(aucVal)
    
  })
  names(aucrocs) <- conds2Test
  
  # If missing
  presentGroups <- names(aucrocs)
  missingGroups <- setdiff(allGroups,
                           presentGroups)
  aucrocs <- c(aucrocs, rep(NA, length(missingGroups)))
  names(aucrocs) <- c(presentGroups, missingGroups)
  
  # Sort
  sortedNames <- sort(names(aucrocs))
  
  # Return
  return(aucrocs[sortedNames])
  
}



getMisMatchedClasses <- function(lassoRes,
                                 dataSplits,
                                 metaData,
                                 predCol="Condition") {
  
  # Get mismatches
  ## Group model
  misMatches <- do.call("rbind",lapply(seq_len(length(lassoRes)), function(i) {
    
    # Name predictionsp
    predList <- setNames(lassoRes[[i]]$rfmodelRes$Predictions$predictions, 
                         names(dataSplits[[i]]$factor_test))
    
    # Get inconsistent
    inconsistentClasses <- getInconsistentClasses(list1=predList,# [names(groupSplits[[i]]$factor_test)], 
                                                  list2=dataSplits[[i]]$factor_test)
    
    # Get same IDs
    sameIDs <- setdiff(names(predList),
                       rownames(inconsistentClasses))
    
    # Get same
    sameData <- data.frame(ID=sameIDs,
                           Predicted=metaData[sameIDs,][[predCol]],
                           Actual=metaData[sameIDs,][[predCol]])
    rownames(sameData) <- sameData$ID
    
    # Check
    if (!all(sameData$Actual == sameData$Predicted)) { stop("Mismatches improperly labelled as the same") }
    
    # Combine
    allClasses <- cbind(Split=i,
                        rbind(sameData,
                              inconsistentClasses))
    
    # Return
    return(allClasses)    
    
  }))
  
  # Add actual class
  misMatches["ActualClass"] <- metaData[misMatches$ID,]$Condition
  misMatches["Mismatch"] <- "Yes"
  
  # Add mistmatched class to bimodal
  subMeta <- metaData[misMatches$ID,]
  if (sum(misMatches$Predicted == misMatches$Actual) > 0) { misMatches[misMatches$Predicted == misMatches$Actual,]["Mismatch"] <- "No" }
  
  # Get missing, untested samples
  unTestedSamples <- setdiff(metaData$Sample_ID,
                             misMatches$ID)
  
  # If any missing samples
  if (length(unTestedSamples) > 0) {
    
    # Create dataframe
    unTestedData <- data.frame(Split=NA,
                               ID=unTestedSamples,
                               Actual=metaData[unTestedSamples,][[predCol]],
                               Predicted=metaData[unTestedSamples,][[predCol]],
                               ActualClass=metaData[unTestedSamples,]$Condition,
                               Mismatch="Not Tested")
    
    # Combined
    misMatches <- rbind(misMatches,
                        unTestedData)
    
  }
  
  # Merge
  metaData["Sample_ID"] <- rownames(metaData)
  misMatches <- merge(misMatches, 
                      metaData,
                      by.x="ID",
                      by.y="Sample_ID")
  
  # Return
  return(misMatches)
  
}

calculateMetrics <- function(misMatches,
                             condName=NULL) {
  
  # If no condition name set
  if (is.null(condName)) {
    
    # Define condition name as not "Other"
    condName <- setdiff(union(unique(misMatches$Actual),
                              unique(misMatches$Predicted)),
                        "Other")
    
  }
  
  # Get confusion matrix
  confMatrix <- table(misMatches[,c("Predicted", "Actual")])
  
  # If any missing
  missingCol <- setdiff(c(condName, "Other"),
                        colnames(confMatrix))
  if (ncol(confMatrix) == 1) {
    
    # Add in dummy column
    confMatrix <- cbind(rep(0, nrow(confMatrix)),
                        confMatrix)
    colnames(confMatrix)[[1]] <- missingCol
    
  }
  
  # If any missing
  missingRow <- setdiff(c(condName, "Other"),
                        rownames(confMatrix))
  if (nrow(confMatrix) == 1) {
    
    # Add in dummy column
    confMatrix <- rbind(rep(0, ncol(confMatrix)),
                        confMatrix)
    rownames(confMatrix)[[1]] <- missingRow
    
  }
  
  # Match order
  confMatrix <- confMatrix[c(condName, "Other"), 
                           c(condName, "Other")]
  
  # Prediction
  misMatches$Predicted <- as.character(misMatches$Predicted)
  metrics <- caret::confusionMatrix(confMatrix)
  
  # If list returned
  if (is.null(nrow(metrics$byClass))) { 
    
    # Get to matrix
    metrics$byClass <- t(as.matrix(metrics$byClass))
    rownames(metrics$byClass) <- condName
    
  }
  
  # Get metrics, make tall and clean
  metricsData <- data.frame(Group=rownames(metrics$byClass),
                            metrics$byClass)
  metricsData <- cleanVars(metricsData)
  
  # Calculate MCC
  mccVal <- mltools::mcc(ifelse(misMatches$Predicted == condName, 1, 0),
                         ifelse(misMatches$Actual == condName, 1, 0))
  metricsData["MCC"] <- mccVal
  
  # If any
  if (length(unique(intersect(misMatches$Predicted, misMatches$Actual))) > 1) {
    
    # Calculate AUC ROC
    rocVal <- suppressMessages(roc(ifelse(misMatches$Predicted == condName, 1, 0),
                                   ifelse(misMatches$Actual == condName, 1, 0)))
    aucVal <- pROC::auc(rocVal)
    
  } else {aucVal <- NA}
  
  # Add AUCROC
  metricsData["AUCROC"] <- aucVal
  
  # Return
  return(metricsData)
  
}


calcOneVsAllMetrics <- function(model,
                                splitData,
                                metaData,
                                multiGroup=FALSE,
                                predCol="BiModal") {
  
  # Loop over splits
  misMatch <- getMisMatchedClasses(lassoRes=model,
                                   dataSplits=splitData,
                                   metaData=metaData,
                                   predCol=predCol)
  
  # If multiple groups
  if (multiGroup) {
    
    # Get groups
    groups <- unique(metaData[[predCol]])
    
    # Loop over 
    splitMetrics <- do.call("rbind", lapply(seq_along(groups), function(i) {
      
      # Loop over splits
      splitMetrics <- do.call("rbind", lapply(seq_along(splitData), function(j) {
        
        # Get split-level data  
        subMisMatch <- misMatch[misMatch$Split == j &
                                  !is.na(misMatch$Split) & 
                                  misMatch$Mismatch != "Not Tested",]

        # Set group-level dummy factor
        subMisMatch$Predicted <- ifelse(subMisMatch$Predicted == groups[[i]],
                                        groups[[i]],
                                        "Other")
        subMisMatch$Actual <- ifelse(subMisMatch$Actual == groups[[i]],
                                     groups[[i]],
                                     "Other")
        
        # Calculate metrics
        metricData <- calculateMetrics(misMatches=subMisMatch,
                                       condName=groups[[i]])
        
        # Return
        return(metricData)
        
      }))
      
    }))
    
    # Else if binary
  } else {
    
    # Loop over splits
    splitMetrics <- do.call("rbind", lapply(seq_along(splitData), function(i) {
      
      # Get split-level data  
      subMisMatch <- misMatch[misMatch$Split == i &
                                !is.na(misMatch$Split) & 
                                misMatch$Mismatch != "Not Tested",]
      
      # Calculate metrics
      metricData <- calculateMetrics(misMatches=subMisMatch)
      
      # Return
      return(metricData)
      
    }))
    
  }
  
  # Return
  return(splitMetrics)
  
}


cleanVars <- function(metricData) {
  
  # Tidy conditions
  metricData$Group <- gsub("Class: ", "", metricData$Group)
  
  # Fix metric names
  colnames(metricData) <- mgsub(colnames(metricData), 
                                c("Pos.Pred.Value", 
                                  "Neg.Pred.Value",
                                  "Detection",
                                  "Prevalence",
                                  "Balanced.Accuracy"),
                                c("PPV", 
                                  "NPV",
                                  "Det",
                                  "Prev",
                                  "Bal Acc"))
  colnames(metricData) <- gsub("[.]", " ", colnames(metricData))
  
  # Remove redundant recall & Precition
  metricData <- metricData[,!colnames(metricData) %in% c("Precision", "Recall", "Prev", "Det Prev", "Det Rate"), drop=FALSE]
  
  # Return
  return(metricData)
  
}


getBestModel <- function(modelFits, 
                         metaData,
                         useDefaultIfNoneBetter=TRUE,
                         baseLine=NULL,
                         chosenConditions=NULL) {
  
  # Get RF models
  rfModels <- lapply(modelFits, function(subData) subData$rfmodelRes$Models)
  
  # Reverse order of nesting
  rfModels <- simplify2array(rfModels)
  rfModels <- split(rfModels, seq_len(nrow(rfModels)))
  
  # Get parameters
  params <- do.call("rbind",lapply(rfModels, function(mod) mod[[1]]$Params))
  
  # Split data
  allData <- do.call("rbind",lapply(seq_along(rfModels), function(i) {
      
      scores <- do.call("rbind",lapply(seq_along(rfModels[[i]]), function(j) {
        
        # Replace sample ID with condii
        splitData <- sapply(rfModels[[i]][[j]]$PerCondAUC, function(cond) cond[[1]])
        data.frame(Model=i,
                   Split=j,
                   AUC=unname(splitData),
                   Condition=names(splitData))
        
      }))
      
  }))
  
  # Add params
  allData <- merge(allData,
                   params,
                   by.x="Model",
                   by.y="row.names")
  
  # Drop NAs (these can happen when not one of a condition in a split)
  allData <- allData[!is.na(allData$AUC),]
  
  # Filter to keep only desired conditions
  if (!is.null(chosenConditions)) { 

    # Filter to desired condition
    allData <- allData[which(allData$Condition %in% chosenConditions),]
    baseLine <- baseLine[which(baseLine$Group %in% chosenConditions),]
  
  }
  
  # Rank by coefficient of variation first - we favour this as the differences in means is so small and we favour low variance across folds,
  # As this most likely corresponds to least variance across datasets
  # By being divisible by mean, CV catches cases with consistently low performance
  coefVar <- function(x) sd(x) / mean(x) # Coef var not sd as it helps removes cases with low means
  cvVal <- aggregate(as.formula(AUC ~ Model),
                     data=allData[,c("Model", "AUC")],
                     FUN=coefVar)
  colnames(cvVal)[[2]] <- "CV"
  
  # Get mean
  meanVal <- aggregate(as.formula(AUC ~ Model),
                       data=allData[,c("Model", "AUC")],
                       FUN=mean)
  colnames(meanVal)[[2]] <- "Mean"
  
  # Merge
  cvVal <- merge(cvVal,
                 meanVal,
                 by="Model")

  # Remove values for CV = 0, ie where SD = 0 (usually associated with consistently bad performance)
  cvVal <- cvVal[cvVal$CV > 0,]
  
  # If baseline provide
  if (!is.null(baseLine)) {
    
    # Keep only cases better than baseline
    cvVal_betterThanDefault <- cvVal[cvVal$Mean > mean(baseLine$AUCROC),]
    
    # If none left-over, return ranger defaults
    if (nrow(cvVal_betterThanDefault) == 0 & useDefaultIfNoneBetter) {
      
      # Get defaults
      out <- list(OptModel="default",
                  OptimalParams=data.frame(Model="default",
                                           nRFTrees=500,
                                           SplitRule="gini",
                                           maxDepth=0,
                                           minNodeSize=1,
                                           mTry="default", # Default to floor(sqrt(p))
                                           SampleFrac=1),
                  ModelRankings=NA ,
                  ModelStats=NA,
                  AllData=allData) 
      
      # Return
      return(out)
    
    # Else if none better than default but we want to proceed anyway
    } else if (nrow(cvVal_betterThanDefault) == 0) { 
      
    cvVal <- cvVal 
    
    # Else use filtered
    } else { cvVal <- cvVal_betterThanDefault } 
    
  
  }

  # Extract smallest CV 
  cvVal <- cvVal[cvVal$CV == base::min(cvVal$CV, na.rm=TRUE),]
  
  # Order by mean for remaining
  cvVal <- cvVal[order(cvVal$Mean, decreasing=TRUE),]
  
  # Get models by metrics
  modelByMetric <- cvVal[,c("Model", "CV")]
  
  # Get optimal model
  optModel <- modelByMetric[1,]$Model
  
  # Get optimal model
  params <- lapply(modelFits, function(subData) subData$rfmodelRes$Params)[[1]][optModel,]
  params["Model"] <- optModel
  
  # Define output
  out <- list(OptModel=optModel,
              OptimalParams=params,
              ModelRankings=modelByMetric,
              ModelStats=cvVal,
              AllData=allData) 
  
  # Return
  return(out)
  
}


