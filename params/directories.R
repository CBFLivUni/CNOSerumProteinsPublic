
# Dependencies
if (!require("here", quietly = TRUE)){ install.packages("here") }
require(here)

# Initialise main directory
here::i_am("params/directories.R")

# Installation and utilities
installDir <- here::here("install")
utilDir <- here::here("utilities")

# Data directories
baseDataDir <- here::here("data", "base")
baseMetaDir <- here::here(baseDataDir, "metadata")
baseDiaDir <- here::here(baseDataDir, "discovery", "DIAProteomics", "directDIA_spectronautProcessed")
baseValidDir <- here::here(baseDataDir, "validation")

dataDir <- here::here("data", runName)
metaDir <- here::here(dataDir, "metadata")
validDir <- here::here(dataDir, "validation")
diaDir <- here::here(dataDir, "discovery", "DIAProteomics", "directDIA_spectronautProcessed")
diaValidDir <- here::here(dataDir, "validation", "DIAProteomics", "directDIA_spectronautProcessed")
impDir <- here::here(diaDir, "imputed")

# Results directories
baseResDir <- here::here("results", "base")
resDir <- here::here("results", runName)
figDir <- here::here(resDir, "figs")
resValidDir <- here::here(resDir, "crossSectional", "validation")
resCrossDiscoDir <- here::here(resDir, "crossSectional", "discovery")
finalDiscoDir <- here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "biomarkers", "finalModel")
resValidPredDir <- here::here(resValidDir, "prediction")

# Create
for (makeDir in c(installDir, 
                  utilDir,
                  metaDir, 
                  diaDir, 
                  impDir,
                  resDir, 
                  figDir,
                  here::here(resDir, "crossSectional", "validation", "eda"),
                  here::here(resCrossDiscoDir, "eda", "missingnessSigTests"),
                  finalDiscoDir,
                  here::here(validDir, "DIAProteomics", "directDIA_spectronautProcessed", "processed","QCs"),
                  here::here(diaDir, "processed","QCs"),
                  here::here(diaDir, "imputed"),
                  here::here(diaDir, "normalised"),
                  here::here(diaValidDir, "processed","QCs"),
                  here::here(diaValidDir, "imputed"),
                  here::here(diaValidDir, "normalised"),
                  here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance"),
                  here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "biomarkers"),
                  here::here(resValidDir, "differentialAbundance"),
                  here::here(validDir, "DIAProteomics", "directDIA_spectronautProcessed", "filtered"),
                  here::here(validDir, "DIAProteomics", "directDIA_spectronautProcessed", "normalised"),
                  here::here(validDir, "DIAProteomics", "directDIA_spectronautProcessed", "imputed"),
                  resValidDir,
                  resValidPredDir)) {
  
  # Create outputs
  dir.create(makeDir,
             recursive=TRUE,
             showWarnings=FALSE)
  
}


