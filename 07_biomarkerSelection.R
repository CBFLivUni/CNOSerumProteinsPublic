
# Source parameters
source("params/runParameters_default.R") # Source default parameters 

# Source directories
source("params/directories.R")



# Install
source("install/installRLibs.R")

# Load libraries
libs_df <- read.table(here::here("install","allLibs2Load.txt"), sep="\t", header=TRUE)
x <- suppressMessages(sapply(libs_df$V1, library, character.only = TRUE)); rm(x)



# Source script that sources files
source(here::here("install", "source_rscripts.R"))

# Source relevant functions
SourceExternalScripts(here::here("utilities","R","functions"), "*.R$", ignore.case=FALSE)



# Source custom parameters
if (file.exists(here::here("params","runParameters_biomarkerRunSpecific.R"))) { source(here::here("params","runParameters_biomarkerRunSpecific.R")) }



# Sort numeric parameters
nIters <- as.numeric(nIters)
freqLassoRuns <- as.numeric(freqLassoRuns)
freqAcrossModels <- as.numeric(freqAcrossModels)
minChosenVars <- as.numeric(minChosenVars)

# Sort boolean parameters
forceCovariatesInModel <- as.logical(ifelse(forceCovariatesInModel == "1" | forceCovariatesInModel == "TRUE" | forceCovariatesInModel == TRUE, TRUE, FALSE))



# Metadata
source(here::here("utilities", "R", "codeModules", "data", "readMetadata_filtered.R"))

# If using unimputed data
if (imputMethod == "None") {
  
  # Read in cyclic Loess-normalised data
  quantTable <- read.table(here::here(diaDir, "normalised", 
                                      "metaCrossSectional_NormalyzerDE_byCondition", "CycLoess-normalized.txt"), sep="\t", header=TRUE, row.names=1)

# Else if using imputed data
} else {

  # Read in chosen imputed data
  quantTable <- read.table(here::here(impDir, paste0("quantTable_imputedWith",imputMethod,".tsv")), sep="\t", header=TRUE, row.names=1)
  
}



# Quantification table
rownames(quantTable) <- gsub("[;]","_",rownames(quantTable))



# Metadata
source(here::here("utilities", "R", "codeModules", "data", "readMetadata_filtered.R"))

# Read in pre-subsetted train data IDs
trainMetaData <- read.table(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom",
                                       paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_allSites.tsv")), header=TRUE, sep="\t")

# Get relevant samples
metaCrossSectional_child <- metaCrossSectional[metaCrossSectional$Age_Group == "Child",]

# Get test data
testMetaData <- metaCrossSectional_child[!metaCrossSectional_child$Sample_ID %in% trainMetaData$Sample_ID,]

# Add years
addMetaCross_child["Years"] <- ( max(as.numeric(ymd(addMetaCross_child$Date.of.sample)), na.rm=TRUE) - as.numeric(ymd(addMetaCross_child$Date.of.sample)) ) / 365 



# Exclude specific sites/conditions
metaCrossSectional_child <- metaCrossSectional_child[!metaCrossSectional_child$Condition %in% excludeConditions,]
metaCrossSectional_child <- metaCrossSectional_child[!metaCrossSectional_child$Site %in% excludeSites,]
metaCrossSectional_child <- metaCrossSectional_child[!metaCrossSectional_child$CondSite %in% excludeCondSites,]

trainMetaData <- trainMetaData[!trainMetaData$Condition %in% excludeConditions,]
trainMetaData <- trainMetaData[!trainMetaData$Site %in% excludeSites,]
trainMetaData <- trainMetaData[!trainMetaData$CondSite %in% excludeCondSites,]

testMetaData <- testMetaData[!testMetaData$Condition %in% excludeConditions,]
testMetaData <- testMetaData[!testMetaData$Site %in% excludeSites,]
testMetaData <- testMetaData[!testMetaData$CondSite %in% excludeCondSites,]

# Fix levels
metaCrossSectional_child$Condition <- as.factor(metaCrossSectional_child$Condition)
trainMetaData$Condition <- as.factor(trainMetaData$Condition)
testMetaData$Condition <- as.factor(testMetaData$Condition)



# If no defined name
if (is.null(mlRunName)) {
  
  # Configure output
  outDir <- here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "biomarkers", "models",
                       imputMethod,
                       # ifelse(batchAdjust, "batchAdj", "unAdj"),
                       ifelse(addCovariates, ifelse(forceCovariatesInModel, "covsPred", "covsSel"), "justProts"),
                       ifelse(univarFSFiltP == "padj", paste0("padj",candidatePThr,"_",univarFiltMethod), paste0("p",candidatePThr,"_",univarFiltMethod)),
                       protSetName,
                       paste0("lfc",log2ThrVal),
                       paste0("lr",freqLassoRuns),
                       paste0("fs",freqAcrossModels),
                       ifelse(useCommonProts_toPredict, "commProts", "diffProts"))
  
  # Define bimodal/group output directory
  biModalOutDir  <- here::here(outDir, "biModal")
  multiGrpOutDir <- here::here(outDir, "multiGroup")
 
# Else 
} else {
  
  # Configure output
  outDir <- here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "biomarkers", "models",
                       mlRunName,
                       imputMethod,
                       # ifelse(batchAdjust, "batchAdj", "unAdj"),
                       ifelse(addCovariates, ifelse(forceCovariatesInModel, "covsPred", "covsSel"), "justProts"),
                       ifelse(univarFSFiltP == "padj", paste0("padj",candidatePThr,"_",univarFiltMethod), paste0("p",candidatePThr,"_",univarFiltMethod)),
                       protSetName,
                       paste0("lfc",log2ThrVal),
                       paste0("lr",freqLassoRuns),
                       paste0("fs",freqAcrossModels),
                       ifelse(useCommonProts_toPredict, "commProts", "diffProts"))
  
  # Define bimodal/group output directory
  multiGrpOutDir <- here::here(outDir, "multiGroup")
  biModalOutDir  <- here::here(outDir, "biModal")
  
}



# Create biomarker directory output
dir.create(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "biomarkers"),
           showWarnings=FALSE,
           recursive=TRUE)

# Save metadata
write.table(trainMetaData, here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "biomarkers",
                                      "allSamplesMetadata_takenForSelection.tsv"), sep="\t", quote=FALSE, row.names=FALSE)



source(here::here("params", "colors.R"))



# Subset to child samples
quantTable <- quantTable[,metaCrossSectional_child$Sample_ID]



# Subset by ID
quantTableTrain <- quantTable[,trainMetaData$Sample_ID]
quantTableTest <- quantTable[,testMetaData$Sample_ID]



# Adjust batches
if (batchAdjust) {
  
  # Generate matrices
  combatDesignM_train <- model.matrix(~ Sex + Age + Condition, data=trainMetaData)
  combatDesignM_test  <- model.matrix(~ Sex + Age + Condition, data=testMetaData)

  # Adjust separately
  quantTableTrain <- ComBat(dat=quantTableTrain[,trainMetaData$Sample_ID], batch=trainMetaData$PLATE, mod=combatDesignM_train, par.prior=TRUE)
  quantTableTest <- ComBat(dat=quantTableTest[,testMetaData$Sample_ID], batch=testMetaData$PLATE, mod=combatDesignM_test, par.prior=TRUE)
  
}



# If we want to perform the analysis excluding proteins with ANY NAs
if (imputMethod == "None") { quantTable <- quantTable[!is.na(rowSums(quantTable)),] }

# Match train/test (have to be the same proteins)
quantTableTrain <- quantTableTrain[which(rownames(quantTableTrain) %in% rownames(quantTable)),]
quantTableTest <- quantTableTest[which(rownames(quantTableTest) %in% rownames(quantTable)),]



# Read in filtered proteins from differential abundance
protSets <- readRDS(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance", 
                               paste0("chosenProteins_featureSelection_pThr",candidatePThr,"_treatLog2Thr",log2ThrVal,".rds")))

# Define filtered proteins
if (protSetName == "All") {
  
  selectedProts <- protSets[[univarFiltMethod]][[protSetName]]
  
} else {
  
  selectedProts <- protSets[[univarFiltMethod]][[univarFSFiltP]][[protSetName]]
  
}

# Report
print(paste0("##### Performing multivariate feature selection using as input:"))
print(paste0("   // ", protSetName, " ( ",length(selectedProts)," Proteins, ",univarFSFiltP,"-filtered ) ",ifelse(batchAdjust, "Batch-adjusted", "Un-adjusted")))

# Clean protein names
protSets <- gsub("[;]","_",protSets)

# Filter proteins
quantTableTrain_filtProts <- quantTableTrain[which(rownames(quantTableTrain) %in% selectedProts),]
quantTableTest_filtProts <- quantTableTest[which(rownames(quantTableTest) %in% selectedProts),]



# If add smaple-level covariates as well as proteins
if (addCovariates) {
  
  # Make relevant column names categorical
  trainMetaData["Sex_asNumber"] <- as.numeric(as.factor(trainMetaData$Sex)) - 1
  trainMetaData["PLATE_asNumber"] <- as.numeric(as.factor(trainMetaData$PLATE)) - 1
  testMetaData["Sex_asNumber"] <- as.numeric(as.factor(testMetaData$Sex)) - 1
  testMetaData["PLATE_asNumber"] <- as.numeric(as.factor(testMetaData$PLATE)) - 1
  
  # Scale age
  trainMetaData["Z_Age"] <- scale(trainMetaData$Age, center=TRUE, scale=TRUE)
  testMetaData["Z_Age"] <- scale(testMetaData$Age, center=TRUE, scale=TRUE)
  
  # If adjusting batches, exclude plate
  if (batchAdjust) {
    
    # Define covariates to add
    covs2Add <- c("Z_Age","Sex_asNumber")
    
  } else {
    
    # Define covariates to add
    covs2Add <- c("Z_Age","Sex_asNumber","PLATE_asNumber")
    
    
  }
    
  # Add covariates
  quantTableTrain_filtProts <- rbind(quantTableTrain_filtProts[,trainMetaData$Sample_ID], t(trainMetaData[,covs2Add]))
  quantTableTest_filtProts <- rbind(quantTableTest_filtProts[,testMetaData$Sample_ID], t(testMetaData[,covs2Add]))
  
} else { covs2Add <- NULL }



# Make sure condition is character
trainMetaData$Condition <- as.character(trainMetaData$Condition)
rownames(trainMetaData) <- trainMetaData$Sample_ID

# Split data train/test
groupSplits <- splitOverIters(trainMetaData, quantTableTrain_filtProts[,trainMetaData$Sample_ID], "Condition", splitProp=trainSplitProp, nIters=nIters, method="Random", considerGroup=FALSE)
biModalSplits <- splitOverIters(trainMetaData, quantTableTrain_filtProts[,trainMetaData$Sample_ID], "BiModal", splitProp=trainSplitProp, nIters=nIters, method="Random", considerGroup=FALSE)



# Decide if we select covariates for prediction
if (forceCovariatesInModel) { covs4Selection <- covs2Add } else { covs4Selection <- NULL }

  # Apply feature selection
  lassoRes_multiGrp <- applyMultiVarSelectionForCrossVal_withModel(inputListwSplits = groupSplits, 
                                                                   covs2Select = covs4Selection,
                                                                   functionSelection = multimodalLasso, 
                                                                   thresMod = as.numeric(freqLassoRuns))
  lassoRes_biModal  <- applyMultiVarSelectionForCrossVal_withModel(inputListwSplits = biModalSplits, 
                                                                   covs2Select = covs4Selection,
                                                                   functionSelection = TwoClassLasso, 
                                                                   thresMod = as.numeric(freqLassoRuns))

  # Get all selected proiteins
  lassoRes_multiGrpRes <- as.data.frame(table(concatenateVarSelec(lassoRes_multiGrp)))
  lassoRes_biModalRes <- as.data.frame(table(concatenateVarSelec(lassoRes_biModal)))
  
  # Get top hits
  chosenVars_group <- as.character(lassoRes_multiGrpRes[lassoRes_multiGrpRes$Freq >= as.numeric(freqAcrossModels),]$Var1)
  chosenVars_biModal <- as.character(lassoRes_biModalRes[lassoRes_biModalRes$Freq >= as.numeric(freqAcrossModels),]$Var1)
  
  # Get chosen covariatees
  chosenCovs_group <- chosenVars_group[which(chosenVars_group %in% covs2Add)]
  chosenCovs_biModal <- chosenVars_biModal[which(chosenVars_biModal %in% covs2Add)]
  
  # Remove covariates
  chosenVars_group <- chosenVars_group[!chosenVars_group %in% covs2Add]
  chosenVars_biModal <- chosenVars_biModal[!chosenVars_biModal %in% covs2Add]
  
  # If using common protein sets fo random forest models
  if (useCommonProts_toPredict) {
    
    if (length(chosenVars_group) > as.numeric(0)) {
      
      # Loop over splits
      for (j in seq_len(length(lassoRes_multiGrp))) {
      
          # Define train/test factor splits
          factor_train <- trainMetaData[rownames(groupSplits[[j]]$train),]$Condition
          factor_test <- trainMetaData[rownames(groupSplits[[j]]$test),]$Condition
          names(factor_train) <- rownames(groupSplits[[j]]$train)
          names(factor_test) <- rownames(groupSplits[[j]]$test)
          
          # Run random forest prediction on split
          rfmodelRes <- RandomForestWithStats(train = groupSplits[[j]]$train[,colnames(groupSplits[[j]]$train) 
                                                                             %in% c(chosenVars_group, covs2Add), drop=FALSE],
                                              test  = groupSplits[[j]]$test[,colnames(groupSplits[[j]]$test) 
                                                                            %in% c(chosenVars_group, covs2Add), drop=FALSE],
                                              factor_train = as.factor(factor_train),
                                              factor_test  = as.factor(factor_test),
                                              covs2Select  = covs2Add,
                                              splitIdx = j)
          
          # Calculate MCC
          rfMCC <- mcc(as.character(rfmodelRes$Predictions$predictions[names(factor_test)]), 
                       as.character(factor_test))
          
          # Remake output
          lassoRes_multiGrp[[j]] <- list("data"=groupSplits[[j]], 
                                         "varSel"=chosenVars_group,
                                         "SampleConditions"=list(Train=factor_train,
                                                                 Test=factor_test),
                                         "VarsModelling"=chosenVars_group,
                                         "rfmodelRes"=rfmodelRes,
                                         "rfMCC"=rfMCC,
                                         "allSelVars"=lassoRes_multiGrp[[j]]$allSelVars,
                                         "finalSelVars"=lassoRes_multiGrp[[j]]$finalSelVars,
                                         "lassoModel"=lassoRes_multiGrp[[j]]$out,
                                         "RFModel"=rfmodelRes$Model)
          
          
      }
        
    } else { 
      
      # Condition if less than a minimum number of chosen proteins
      print(paste0("Minimum chosen variables for Group Model is < ",minChosenVars,": Switching to differential selected proteins (see 'diffSplitProts/' as output)" ))
      multiGrpOutDir <- gsub("commProts", "diffSplitProts", multiGrpOutDir)
      useCommonProts_toPredict <- FALSE
      
      # Set MCC to NA
      lassoRes_multiGrp[[j]]$rfMCC <- mcc(as.character(rfmodelRes$Predictions$predictions[names(factor_test)]), 
                                          as.character(factor_test))
      
    }
    
    if (length(chosenVars_biModal) > as.numeric(minChosenVars)) {
      
      # Loop over splits
      for (j in seq_len(length(lassoRes_biModalRes))) {
        
          # Define train/test factor splits
          factor_train <- ifelse(trainMetaData[rownames(biModalSplits[[j]]$train),]$Condition == "CNO", "CNO", "Other")
          factor_test <- ifelse(trainMetaData[rownames(biModalSplits[[j]]$test),]$Condition == "CNO", "CNO", "Other")
          names(factor_train) <- rownames(biModalSplits[[j]]$train)
          names(factor_test) <- rownames(biModalSplits[[j]]$test)
          
          # Run random forest prediction on split
          rfmodelRes <- RandomForestWithStats(train = biModalSplits[[j]]$train[,colnames(biModalSplits[[j]]$train) %in% c(chosenVars_biModal, covs2Add), drop=FALSE],
                                              test  = biModalSplits[[j]]$test[,colnames(biModalSplits[[j]]$test) %in% c(chosenVars_biModal, covs2Add), drop=FALSE],
                                              factor_train = as.factor(factor_train),
                                              factor_test = as.factor(factor_test),
                                              covs2Select = covs2Add,
                                              splitIdx = j)
          
          # Calculate MCC
          rfMCC <- mcc(as.character(rfmodelRes$Predictions$predictions[names(factor_test)]), 
                       as.character(factor_test))
          
          # Remake output
          lassoRes_biModal[[j]] <- list("data"=biModalSplits[[j]], 
                                        "varSel"=chosenVars_biModal,
                                        "SampleConditions"=list(Train=factor_train,
                                                                Test=factor_test),
                                        "VarsModelling"=chosenVars_biModal,
                                        "rfmodelRes"=rfmodelRes,
                                        "rfMCC"=rfMCC,
                                        "allSelVars"=lassoRes_biModal[[j]]$allSelVars,
                                        "finalSelVars"=lassoRes_biModal[[j]]$finalSelVars,
                                        "lassoModel"=lassoRes_biModal[[j]]$out,
                                        "RFModel"=rfmodelRes$Model)
        
      }
      
    } else {
      
      # Condition if less than a minimum number of chosen proteins
      print(paste0("Minimum chosen variables for Bimodal Model is < ",minChosenVars,": Switching to differential selected proteins (see 'diffSplitProts/' as output)" ))
      biModalOutDir <- gsub("commProts", "diffSplitProts", biModalOutDir)
      useCommonProts_toPredict <- FALSE
      
      # Set MCC to NA
      lassoRes_biModal[[j]]$rfMCC <- mcc(as.character(rfmodelRes$Predictions$predictions[names(factor_test)]), 
                                         as.character(factor_test))
      
    }
    
  }

# Subset split data if necessary
biModalSplits <- biModalSplits[sapply(lassoRes_biModal, function(model) !is.null(model))]
groupSplits <- groupSplits[sapply(lassoRes_multiGrp, function(model) !is.null(model))]

# Drop NULLs
lassoRes_biModal  <- lassoRes_biModal[sapply(lassoRes_biModal, function(model) !is.null(model))]
lassoRes_multiGrp <- lassoRes_multiGrp[sapply(lassoRes_multiGrp, function(model) !is.null(model))]



# Create output
dir.create(here::here(multiGrpOutDir, "pca"), showWarnings=FALSE, recursive=TRUE)
dir.create(here::here(biModalOutDir, "pca"), showWarnings=FALSE, recursive=TRUE)

# Print
print(paste0("##### Saving to \\"))
print(paste0(multiGrpOutDir))
print(paste0(biModalOutDir))

# Save metadata
write.table(trainMetaData, here::here(outDir,
                                      "allSampleMeta_taken4Select.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Save split data
saveRDS(groupSplits, here::here(multiGrpOutDir, "mg_dataSplits.rds"))
saveRDS(biModalSplits,  here::here(biModalOutDir, "bm_dataSplits.rds"))



# Get 
modelInfo <- c(imputMethod,
               ifelse(batchAdjust, "batchAdj", "unAdj"),
               ifelse(addCovariates, ifelse(forceCovariatesInModel, "forceCovsPredict", "incCovsSelect"), "justProts"),
               ifelse(useCommonProts_toPredict, "commProts", "diffSplitProts"),
               ifelse(univarFSFiltP == "padj", "padj", "pval"),
               freqLassoRuns,
               freqAcrossModels,
               candidatePThr,
               log2ThrVal,
               univarFiltMethod,
               protSetName)



# Loop over models
for (i in seq_len(length(lassoRes_multiGrp))) {
  
  # Loop over splits
  lassoRes_multiGrp[[i]][["Parameterisation"]] <- modelInfo
  
}
for (i in seq_len(length(lassoRes_biModal))) {

  lassoRes_biModal[[i]][["Parameterisation"]] <- modelInfo
  
}



# Save model fits
saveRDS(lassoRes_multiGrp, here::here(multiGrpOutDir, "mg_lassoRFModel.rds"))
saveRDS(lassoRes_biModal,  here::here(biModalOutDir,  "bm_lassoRFModel.rds"))

