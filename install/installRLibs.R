
# Increase timeout for installing large libraries
options(timeout = 100000)

# Install biocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#function to install libararies from bioconductor
installPackage <- function(libName) {
  if(libName %in% rownames(installed.packages()) == FALSE){
    BiocManager::install(libName,ask = FALSE)
  }}

# Install here package
if (!require("here", quietly = TRUE)){ install.packages("here") }
require(here)

# Define here
here::i_am("README.md")

#read the libraries needed
packagesToInstall <- read.delim(here::here("install", "Rlibs.txt"),header=F,stringsAsFactors = F)

# Install packages
x <- sapply(packagesToInstall[,1],installPackage)

# Install gitHUb
if (!require("devtools", quietly = TRUE)) { install.packages("devtools") }
if (!require("NormalyzerDE", quietly = TRUE)) { devtools::install_github("ComputationalProteomics/NormalyzerDE") }
if (!require("INDperform", quietly = TRUE)) { remotes::install_github("saskiaotto/INDperform") }

# Install pacakges for NAguideR: See https://github.com/wangshisheng/NAguideR
if(!require(pacman)) install.packages("pacman")
# pacman::p_load(devtools, shiny, shinyBS, shinyjs, shinyWidgets, DT, gdata, ggplot2, glmnet, reshape2, ggsci, openxlsx, data.table, DT, raster, Metrics, vegan, tidyverse, ggExtra, cowplot, Amelia, e1071, impute, SeqKnn, pcaMethods, norm, imputeLCMD, VIM, rrcovNA, mice, missForest)
if (!require(DreamAI)) install_github("WangLab-MSSM/DreamAI/Code")

# Manually install NAguideR packages
if(!require(SeqKnn)) install.packages(here::here("install", "SeqKnn_1.0.1.tar.gz"), repos = NULL,type="source")
if(!require(GMSimpute)) install.packages(here::here("install", "GMSimpute_0.0.1.1.tar.gz"), repos = NULL,type="source")

# Install NA guideR
if(!require(NAguideR)) devtools::install_github("wangshisheng/NAguideR")

# Install specific missForest version - imp4p uses missForest for impute.RF. missForest switched to ranger from randomForest > v1.4
detach("package:DreamAI", unload = TRUE)
if ("imp4p" %in% .packages()) detach("package:imp4p", unload = TRUE)
detach("package:missForest", unload = TRUE)
if (!require(missForest) | as.character(packageVersion("missForest")) != "1.4") remotes::install_version("missForest", version = "1.4", repos = "http://cran.r-project.org")

# Add github installs
packagesToInstall <- rbind(packagesToInstall, data.frame(V1=c("NormalyzerDE","SeqKnn", "GMSimpute", "DreamAI", "pacman", "NormalyzerDE", "INDperform")))

# Save
write.table(packagesToInstall, here::here("install", "allLibs2Load.txt"), sep="\t", quote=FALSE)
