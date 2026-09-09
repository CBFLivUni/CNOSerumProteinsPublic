

# Read in metadata
metaCrossSectional <- read.table(here::here(metaDir, "preSchool_crossSectional_tallFormatClean.tsv"), sep="\t", header=TRUE)
addMetaCross <- read.table(here::here(metaDir, "preSchool_additionalMetadata_commonToCrossSectional.tsv"), sep="\t", header=TRUE)
addMetaCross_child <- addMetaCross[addMetaCross$Age_Group == "Child",]
rownames(metaCrossSectional) <- metaCrossSectional$Sample_ID

# Read annotation
proteinAnnotData <- read.table(here::here(baseDiaDir, "spectronaut", "20240716_091619_crossSectionalOnly", "crossSectionalOnly_Report_BGS Analysis Grid View Report (Pivot).tsv"), sep="\t", header=TRUE)
proteinAnnotData <- proteinAnnotData[!duplicated(proteinAnnotData$PG.ProteinAccessions),c("PG.ProteinAccessions","PG.ProteinDescriptions","PG.ProteinNames")]
proteinAnnotData$PG.ProteinNames <- gsub("_HUMAN","",proteinAnnotData$PG.ProteinNames)