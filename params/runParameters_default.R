
# Global run parameters
reGenerateFiles <- TRUE
runName <- "publicationVersion" # 
excludeAdultSamples <- TRUE

# Fiiltration parameters
maxPropMissing <- 0.2
maxPropMissing_samples <- 0.2
stringentMissValPPCA <- FALSE
maxPropMissingPerGroup <- 0.65 # This is lax, and we still get nothing
maxPropMissingOverall <- 0.4

# Differential Abundance
runVarAnalysis <- FALSE
diffAbundance_nCVs <- 10
cvProtThr <- 1

# Data subsetting
splitTrainTest <- TRUE
subsetTrainTest <- TRUE
validSplitSeed <- 1
trainSplitProp <- 0.8
subsetSeed <- 12
splitMethod <- "Random"
normAfterSplit <- FALSE # TRUE Might not be viable as it makes train/test uncomparable

# ML & Biomarker selection
mlRunName <- "dropDeplProts"
log2ThrVal <- 0
candidatePThr <- 0.05
univarFSFiltP <- "padj"
univarFiltMethod <- "limma"
protSetName <- "DresdenOnly"
imputMethod <- "bPCA"
batchAdjust <- FALSE
addCovariates <- FALSE
forceCovariatesInModel <- FALSE
covs2Select <- NULL
handleCrohnCoDiagnosis <- "keep" # "flip" # "keep" # "remove"
nIters <- 10
minChosenVars <- 0 
useCommonProts_toPredict <- TRUE
excludeSites <- c()
excludeConditions <- c("JIA")
excludeCondSites <- c("JIA_Liverpool", "JIA_Dresden")
excludeDeplProts <- TRUE

# Model evaluation
nReps_forDecoyShuffle <- 25
nUnderSamples <- 25

# Model tuning
tuneModelParams <- FALSE
useDiscoveryTunedParams <- TRUE
stepSizeParams <- 2
nRFTrees <- c(1500, 2000, 2500)
maxDepths <- seq(0, 4)
minNodeSize <- seq(1, 4, 1)
sampleFrac <- 0.632 # Use 0.632 - to ensure full dataset is not used
splitRule <- c("extratrees", "gini") # Extra-trees less prone to overfitting due to random splitting. So should reduce risk on small datasets

# Validation data
includeHPP <- FALSE
dateFilterProts <- FALSE
sigFilterProts <- FALSE
siteFiltProts <- FALSE
covs2Select4Valid <- NULL
dropConditions4Valid <- NULL
dropSites4Valid <- NULL
nReps2DownSample <- 25
chosenMetrics <- c("MCC") # For random subsamples
