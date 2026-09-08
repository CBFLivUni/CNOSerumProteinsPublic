
# Source metadata
here::i_am("06_differentialAbundance.qmd")
source(here::here("utilities", "R", "codeModules", "boilerPlate", "loadParamsLibrariesFunctions.R"))

# Create output
dir.create(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate"),
           recursive=TRUE)



### Externally provide FDR and Log2Thr
args <- commandArgs(trailingOnly = TRUE)
candidatePThr <- as.numeric(args[[1]])
log2ThrVal <- as.numeric(args[[2]])



# Save as RDS
plotList <- readRDS(here::here(figDir, 
                    "05_dataImputation_figs.rds"))
 


# Train/test split output
splitDir <- here::here(resCrossDiscoDir, paste0("trainSplitProp", trainSplitProp))
splitDir <- here::here(splitDir, "sampleRandom")



# Metadata
source(here::here("utilities", "R", "codeModules", "data", "readMetadata_filtered.R"))

# Abundance tables
source(here::here("utilities", "R", "codeModules", "data", "unFiltered_quantTables.R"))

# Read in chosen normalised data
normQuantTable <- read.table(here::here(diaDir, "normalised", "crossSectional_quantTable_normalisedCyclicLoess.tsv"),
                             sep="\t", header=TRUE, row.names=1)

# Drop QC samples
normQuantTable <- normQuantTable[,!grepl("QC", colnames(normQuantTable))]
dim(normQuantTable)



source(here::here("params", "colors.R"))



# Filter cross-sectional data
metaCrossSectional <- metaCrossSectional[colnames(normQuantTable),]
addMetaCross <- addMetaCross[which(addMetaCross$Sample_ID %in% colnames(normQuantTable)),]



metaCrossSectional$Condition <- gsub("Infectious ","Inf_",metaCrossSectional$Condition)



# get child-only IDs
childIDs <- metaCrossSectional[metaCrossSectional$Age_Group == "Child" & !is.na(metaCrossSectional$Age),]$Sample_ID
length(childIDs)

# Subet metadata
metaCrossSectional_child <- metaCrossSectional[childIDs,]

# Change names
trainMetaData <- metaCrossSectional_child
trainIDs <- trainMetaData$Sample_ID
  
# Subset abundance matrix
trainQuantTable <- normQuantTable[,trainIDs]
dim(trainQuantTable)



# Subset metadata
trainMetaData_noVbilt <- trainMetaData[trainMetaData$Site != "Vanderbilt",]
trainMetaData_incVblt <- trainMetaData[trainMetaData$Site == "Vanderbilt",]
trainMetaData_dresden <- trainMetaData[trainMetaData$Site == "Dresden",]
trainMetaData_noLivpl <- trainMetaData[trainMetaData$Site != "Liverpool",]
trainMetaData_liverpl <- trainMetaData[trainMetaData$Site == "Liverpool",]
trainMetaData_lvplH <- rbind(trainMetaData_dresden, # Dresden + Liverpool Healthy
                             trainMetaData_incVblt,
                             trainMetaData[trainMetaData$Site == "Liverpool" & trainMetaData$Condition == "Healthy",])
trainMetaData_lvVbH <- rbind(trainMetaData_dresden, # Dresden + Vanderbilt Healthy + Liverpool Healthy
                           trainMetaData[trainMetaData$Site == "Vanderbilt" & trainMetaData$Condition == "Healthy",],
                           trainMetaData[trainMetaData$Site == "Liverpool" & trainMetaData$Condition == "Healthy",])



# Set to factor
# addMetaCross_child$Condition <- as.factor(addMetaCross_child$Condition)
  
# Get date as years since first
addMetaCross_child["DateOfSample_asNumber"] <- (base::max(as.numeric(ymd(addMetaCross_child$Date.of.sample)), na.rm=TRUE) - as.numeric(ymd(addMetaCross_child$Date.of.sample))) / 365
addMetaCross_child["Date_asCategorical"] <- cut(addMetaCross_child$DateOfSample_asNumber, breaks = seq(0, 20, length.out = 21), include.lowest = TRUE, labels = FALSE)

# Ensure all additional metadata samples are in quant table
addMetaCross_child <- addMetaCross_child[which(addMetaCross_child$Sample_ID %in% colnames(trainQuantTable)),]
rownames(addMetaCross_child) <- addMetaCross_child$Sample_ID

# Remove NAs for date and age
addMetaCross_child <- addMetaCross_child[!is.na(addMetaCross_child$Age) & !is.na(addMetaCross_child$Date.of.sample),]

# Split additional metadata to subsets
addMetaData <- addMetaCross_child # Dresden + Vanderbilt + Liverpool
addMetaData_drsdn <- addMetaCross_child[addMetaCross_child$Site == "Dresden",] # Just Dresden
addMetaData_drsVb <- addMetaCross_child[addMetaCross_child$Site != "Liverpool",] # Dresden + Vanderbilt
addMetaData_drsLp <- addMetaCross_child[addMetaCross_child$Site != "Vanderbilt",] # Dresden + Liverpool
addMetaData_lvplH <- rbind(addMetaData_drsVb, # Dresden + Vanderbilt + Liverpool Healthy
                           addMetaCross_child[addMetaCross_child$Site == "Liverpool" & addMetaCross_child$Condition == "Healthy",])
addMetaData_lvVbH <- rbind(addMetaData_drsdn, # Dresden + Vanderbilt Healthy + Liverpool Healthy
                           addMetaCross_child[addMetaCross_child$Site == "Vanderbilt" & addMetaCross_child$Condition == "Healthy",],
                           addMetaCross_child[addMetaCross_child$Site == "Liverpool" & addMetaCross_child$Condition == "Healthy",])

# Check
table(addMetaData[,c("Site","Condition")])
table(addMetaData_drsdn[,c("Site","Condition")])
table(addMetaData_drsVb[,c("Site","Condition")])
table(addMetaData_drsLp[,c("Site","Condition")])
table(addMetaData_lvplH[,c("Site","Condition")])
table(addMetaData_lvVbH[,c("Site","Condition")])

# Save




# Make Healthy the reference
## Not considering date
trainMetaData$Condition <- as.factor(trainMetaData$Condition)
trainMetaData_noVbilt$Condition <- as.factor(trainMetaData_noVbilt$Condition)
trainMetaData_incVblt$Condition <- as.factor(trainMetaData_incVblt$Condition)
trainMetaData_dresden$Condition <- as.factor(trainMetaData_dresden$Condition)
trainMetaData_noLivpl$Condition <- as.factor(trainMetaData_noLivpl$Condition)
trainMetaData_liverpl$Condition <- as.factor(trainMetaData_liverpl$Condition)
trainMetaData_lvplH$Condition <- as.factor(trainMetaData_lvplH$Condition)
trainMetaData_lvVbH$Condition <- as.factor(trainMetaData_lvVbH$Condition)
if (length(excludeConditions) == 0 | is.null(length(excludeConditions))) { trainMetaData$Condition <- relevel(trainMetaData$Condition, "JIA") }
trainMetaData_noVbilt$Condition <- relevel(trainMetaData_noVbilt$Condition, "Healthy")
trainMetaData_incVblt$Condition <- relevel(trainMetaData_incVblt$Condition, "Healthy")
if (length(excludeConditions) == 0 | is.null(length(excludeConditions))) { trainMetaData_dresden$Condition <- relevel(trainMetaData_dresden$Condition, "JIA") }
trainMetaData_noLivpl$Condition <- relevel(trainMetaData_noLivpl$Condition, "Healthy")
trainMetaData_liverpl$Condition <- relevel(trainMetaData_liverpl$Condition, "Healthy")
trainMetaData_lvplH$Condition <- relevel(trainMetaData_lvplH$Condition, "Healthy")
trainMetaData_lvVbH$Condition <- relevel(trainMetaData_lvVbH$Condition, "Healthy")

## With date
addMetaData$Condition <- relevel(as.factor(addMetaData$Condition), "Healthy")
if (length(excludeConditions) == 0 | is.null(length(excludeConditions))) { addMetaData_drsdn$Condition <- relevel(as.factor(addMetaData_drsdn$Condition), "JIA") }
addMetaData_drsVb$Condition <- relevel(as.factor(addMetaData_drsVb$Condition), "Healthy")
addMetaData_drsLp$Condition <- relevel(as.factor(addMetaData_drsLp$Condition), "Healthy")
addMetaData_lvplH$Condition <- relevel(as.factor(addMetaData_lvplH$Condition), "Healthy")
addMetaData_lvVbH$Condition <- relevel(as.factor(addMetaData_lvVbH$Condition), "Healthy")



# If we want to take a train subset of the data
if (splitTrainTest) {
  
  print(paste0("##### Subsetting metadata to Train = ",trainSplitProp*100,"%"))
  
  # Perform split
  ## Not considering date
  set.seed(validSplitSeed)
  trainMetaData         <- splitTrainData(meta=trainMetaData, 
                                          propSplit=trainSplitProp, seedVal=subsetSeed, sampleMethod="Random")
  trainMetaData_noVbilt <- trainMetaData[trainMetaData$Site != "Vanderbilt",]
  trainMetaData_incVblt <- trainMetaData[trainMetaData$Site == "Vanderbilt",]
  trainMetaData_dresden <- trainMetaData[trainMetaData$Site == "Dresden",]
  trainMetaData_noLivpl <- trainMetaData[trainMetaData$Site != "Liverpool",]
  trainMetaData_liverpl <- trainMetaData[trainMetaData$Site == "Liverpool",]
  trainMetaData_lvplH <- rbind(trainMetaData_dresden, # Dresden + Liverpool Healthy
                               trainMetaData_incVblt,
                               trainMetaData[trainMetaData$Site == "Liverpool" & trainMetaData$Condition == "Healthy",])
  trainMetaData_lvVbH <- rbind(trainMetaData_dresden, # Dresden + Vanderbilt Healthy + Liverpool Healthy
                             trainMetaData[trainMetaData$Site == "Vanderbilt" & trainMetaData$Condition == "Healthy",],
                             trainMetaData[trainMetaData$Site == "Liverpool" & trainMetaData$Condition == "Healthy",])
  
  # Samples with date
  set.seed(validSplitSeed)
  trainAddMetaData       <- splitTrainData(meta=addMetaData, 
                                           propSplit=trainSplitProp, seedVal=subsetSeed, sampleMethod=splitMethod, dateCol="DateOfSample_asNumber")
  trainAddMetaData_drsdn <- trainAddMetaData[trainAddMetaData$Site == "Dresden",] # Just Dresden
  trainAddMetaData_drsVb <- trainAddMetaData[trainAddMetaData$Site != "Liverpool",] # Dresden + Vanderbilt
  trainAddMetaData_drsLp <- trainAddMetaData[trainAddMetaData$Site != "Vanderbilt",] # Dresden + Liverpool
  trainAddMetaData_lvplH <- rbind(addMetaData_drsVb, # Dresden + Vanderbilt + Liverpool Healthy
                                  trainAddMetaData[trainAddMetaData$Site == "Liverpool" & trainAddMetaData$Condition == "Healthy",])
  trainAddMetaData_lvVbH <- rbind(addMetaData_drsdn, # Dresden + Vanderbilt Healthy + Liverpool Healthy
                                  trainAddMetaData[trainAddMetaData$Site == "Vanderbilt" & trainAddMetaData$Condition == "Healthy",],
                                  trainAddMetaData[trainAddMetaData$Site == "Liverpool" & trainAddMetaData$Condition == "Healthy",])

# Else just keep same number of samples
} else {
  
  # Change names
  trainAddMetaData       <- addMetaData
  trainAddMetaData_drsdn <- addMetaData_drsdn
  trainAddMetaData_drsVb <- addMetaData_drsVb
  trainAddMetaData_drsLp <- addMetaData_drsLp
  trainAddMetaData_lvplH <- addMetaData_lvplH
  trainAddMetaData_lvVbH <- addMetaData_lvVbH
  
}



# If normalise after train/test split
if (normAfterSplit) { 
  
  # Re-normalise with Cyclic Loesss
  trainQuantTable <- limma::normalizeCyclicLoess(log2(quantTable_crossSectional[rownames(normQuantTable),trainAddMetaData$sample]), method="fast")
  
}



# Split data
addMetaData["SplitSet"] <- "Test"
addMetaData_drsdn["SplitSet"] <- "Test"
addMetaData_drsVb["SplitSet"] <- "Test"
addMetaData_drsLp["SplitSet"] <- "Test"
addMetaData_lvplH["SplitSet"] <- "Test"
addMetaData_lvVbH["SplitSet"] <- "Test"
addMetaData[which(addMetaData$Sample_ID %in% trainAddMetaData$Sample_ID),]["SplitSet"] <- "Train"
addMetaData_drsdn[which(addMetaData_drsdn$Sample_ID %in% trainAddMetaData_drsdn$Sample_ID),]["SplitSet"] <- "Train"
addMetaData_drsVb[which(addMetaData_drsVb$Sample_ID %in% trainAddMetaData_drsVb$Sample_ID),]["SplitSet"] <- "Train"
addMetaData_drsLp[which(addMetaData_drsLp$Sample_ID %in% trainAddMetaData_drsLp$Sample_ID),]["SplitSet"] <- "Train"
addMetaData_lvplH[which(addMetaData_lvplH$Sample_ID %in% trainAddMetaData_lvplH$Sample_ID),]["SplitSet"] <- "Train"
addMetaData_lvVbH[which(addMetaData_lvVbH$Sample_ID %in% trainAddMetaData_lvVbH$Sample_ID),]["SplitSet"] <- "Train"

# Plot
p <- ggplot(addMetaData[addMetaData$SplitSet == "Train",], 
       aes(x=Condition, y=DateOfSample_asNumber, fill=Condition)) + 
  geom_point(size=10, alpha=0.65, shape=21, position = position_nudge(x = 0.1)) + 
  geom_point(data=addMetaData[addMetaData$SplitSet == "Test",], aes(x=Condition, y=DateOfSample_asNumber),
             size=10, alpha=0.85, shape=21, fill="white", position = position_nudge(x = -0.1)) + 
  stat_summary(data=addMetaData[addMetaData$SplitSet == "Train",],
               fun.y = median, geom = "point", shape=8, size=16, position = position_nudge(x = 0.1),
               aes(ymax = ..y.., ymin = ..y.., group = Condition),
               width = 1, linetype = "solid", show.legend=FALSE) +
  stat_summary(data=addMetaData[addMetaData$SplitSet == "Test",],
               fun.y = median, geom = "point", shape=8, size=16, position = position_nudge(x = -0.1),
               aes(ymax = ..y.., ymin = ..y.., group = Condition),
               width = 1, linetype = "solid", show.legend=FALSE) +
  scale_fill_manual(values=condCols) +
  ylab("Age of Sample\n(Years)") + xlab("") +
  facet_grid(cols=vars(Site), space="free", scale="free") +
  theme_classic(base_size=48) + 
  theme(legend.position="none",
        axis.title=element_text(size=88),
        strip.text=element_text(size=64))
ggplot(addMetaData_drsdn, aes(x=Condition, y=DateOfSample_asNumber, fill=Condition)) + 
  geom_point(size=8, alpha=0.65, shape=21) + 
  geom_point(data=addMetaData_drsdn[addMetaData_drsdn$SplitSet == "Test",], aes(x=Condition, y=DateOfSample_asNumber),
             size=8, alpha=0.85, shape=21, fill="white") + 
  stat_summary(fun.y = median, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y.., group = Condition),
               width = 1, linetype = "solid") +
  scale_fill_manual(values=condCols) +
  ylab("Age of Sample (Years)") +
  facet_grid(cols=vars(Site), space="free", scale="free") +
  theme_classic(base_size=32) + theme(legend.position="none")
ggplot(addMetaData_drsVb, aes(x=Condition, y=DateOfSample_asNumber, fill=Condition)) + 
  geom_point(size=4, alpha=0.65, shape=21) + 
  geom_point(data=addMetaData_drsVb[addMetaData_drsVb$SplitSet == "Test",], aes(x=Condition, y=DateOfSample_asNumber),
             size=4, alpha=0.85, shape=21, fill="white") + 
  stat_summary(fun.y = median, geom = "point", 
               aes(ymax = ..y.., ymin = ..y.., group = Condition),
               width = 1, linetype = "solid") +
  scale_fill_manual(values=condCols) +
  ylab("Age of Sample (Years)") +
  facet_grid(cols=vars(Site), space="free", scale="free") +
  theme_classic(base_size=32) + theme(legend.position="none")
ggplot(addMetaData_drsLp, aes(x=Condition, y=DateOfSample_asNumber, fill=Condition)) + 
  geom_point(size=4, alpha=0.65, shape=21) + 
  geom_point(data=addMetaData_drsLp[addMetaData_drsLp$SplitSet == "Test",], aes(x=Condition, y=DateOfSample_asNumber),
             size=4, alpha=0.85, shape=21, fill="white") + 
  stat_summary(fun.y = median, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y.., group = Condition),
               width = 1, linetype = "solid") +
  scale_fill_manual(values=condCols) +
  ylab("Age of Sample (Years)") +
  facet_grid(cols=vars(Site), space="free", scale="free") +
  theme_classic(base_size=32) + theme(legend.position="none")
ggplot(addMetaData_lvplH, aes(x=Condition, y=DateOfSample_asNumber, fill=Condition)) + 
  geom_point(size=4, alpha=0.65, shape=21) + 
  geom_point(data=addMetaData_lvplH[addMetaData_lvplH$SplitSet == "Test",], aes(x=Condition, y=DateOfSample_asNumber),
             size=4, alpha=0.85, shape=21, fill="white") + 
  stat_summary(fun.y = median, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y.., group = Condition),
               width = 1, linetype = "solid") +
  scale_fill_manual(values=condCols) +
  ylab("Age of Sample (Years)") +
  facet_grid(cols=vars(Site), space="free", scale="free") +
  theme_classic(base_size=32) + theme(legend.position="none")
ggplot(addMetaData_lvVbH, aes(x=Condition, y=DateOfSample_asNumber, fill=Condition)) + 
  geom_point(size=4, alpha=0.65, shape=21) + 
  geom_point(data=addMetaData_lvVbH[addMetaData_lvVbH$SplitSet == "Test",], aes(x=Condition, y=DateOfSample_asNumber),
             size=4, alpha=0.85, shape=21, fill="white") + 
  stat_summary(fun.y = median, geom = "errorbar", 
               aes(ymax = ..y.., ymin = ..y.., group = Condition),
               width = 1, linetype = "solid") +
  scale_fill_manual(values=condCols) +
  ylab("Age of Sample (Years)") +
  facet_grid(cols=vars(Site), space="free", scale="free") +
  theme_classic(base_size=32) + theme(legend.position="none")

# Save
ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate",
                  "sampleDateDistribution_allSites_discoveryData_trainTest.png"),
       p,
       units="px", width=8000, height=3500)

# Add to list
newList <- list(p)
names(newList) <- c("sampleDateDistribution")
plotList <- append(plotList,
                   newList)



# Z-scale age
trainMetaData["Z_Age"] <- scale(trainMetaData$Age)
trainMetaData_noVbilt["Z_Age"] <- scale(trainMetaData_noVbilt$Age)
trainMetaData_incVblt["Z_Age"] <- scale(trainMetaData_incVblt$Age)
trainMetaData_dresden["Z_Age"] <- scale(trainMetaData_dresden$Age)
trainMetaData_noLivpl["Z_Age"] <- scale(trainMetaData_noLivpl$Age)
trainMetaData_liverpl["Z_Age"] <- scale(trainMetaData_liverpl$Age)
trainMetaData_lvplH["Z_Age"] <- scale(trainMetaData_lvplH$Age)
trainMetaData_lvVbH["Z_Age"] <- scale(trainMetaData_lvVbH$Age)

# Z-scale only date
trainAddMetaData["Z_Date"] <- scale(trainAddMetaData$DateOfSample_asNumber)
trainAddMetaData_drsdn["Z_Date"] <- scale(trainAddMetaData_drsdn$DateOfSample_asNumber)
trainAddMetaData_drsVb["Z_Date"] <- scale(trainAddMetaData_drsVb$DateOfSample_asNumber)
trainAddMetaData_drsLp["Z_Date"] <- scale(trainAddMetaData_drsLp$DateOfSample_asNumber)
trainAddMetaData_lvplH["Z_Date"] <- scale(trainAddMetaData_lvplH$DateOfSample_asNumber)
trainAddMetaData_lvVbH["Z_Date"] <- scale(trainAddMetaData_lvVbH$DateOfSample_asNumber)
trainAddMetaData["Z_Age"] <- scale(trainAddMetaData$Age)
trainAddMetaData_drsdn["Z_Age"] <- scale(trainAddMetaData_drsdn$Age)
trainAddMetaData_drsVb["Z_Age"] <- scale(trainAddMetaData_drsVb$Age)
trainAddMetaData_drsLp["Z_Age"] <- scale(trainAddMetaData_drsLp$Age)
trainAddMetaData_lvplH["Z_Age"] <- scale(trainAddMetaData_lvplH$Age)
trainAddMetaData_lvVbH["Z_Age"] <- scale(trainAddMetaData_lvVbH$Age)



# Create output directories
dir.create(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom"), showWarnings=FALSE, recursive=TRUE)
dir.create(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate"), showWarnings=FALSE, recursive=TRUE)

# Save data
write.table(trainMetaData, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
            paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_allSites.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainMetaData_noVbilt, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
            paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_noVanderbilt.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainMetaData_incVblt, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
            paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_justVanderbilt.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainMetaData_dresden, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
            paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_justDresden.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainMetaData_noLivpl, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
            paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_noLiverpool.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainMetaData_liverpl, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
            paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_justLiverpool.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainMetaData_lvplH, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
            paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_justDresdenVanderbilt_andLiverpoolHealthy.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainAddMetaData_lvVbH, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
            paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_justDreden_andVanderbiltLiverpoolHealthy.tsv")), sep="\t", quote=FALSE, row.names=TRUE)

write.table(trainAddMetaData, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", 
                       paste0("trainData_justChildSamples_samplesWithDate_splitProp",trainSplitProp,"_allSites.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainAddMetaData_drsdn, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", 
                       paste0("trainData_justChildSamples_samplesWithDate_splitProp",trainSplitProp,"_justDresden.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainAddMetaData_drsVb, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", 
                       paste0("trainData_justChildSamples_samplesWithDate_splitProp",trainSplitProp,"_justDresdenVanderbilt.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainAddMetaData_drsLp, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", 
                       paste0("trainData_justChildSamples_samplesWithDate_splitProp",trainSplitProp,"_justDredenLiverpool.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainAddMetaData_lvplH, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", 
                       paste0("trainData_justChildSamples_samplesWithDate_splitProp",trainSplitProp,"_justDresdenVanderbilt_andLiverpoolHealthy.tsv")), sep="\t", quote=FALSE, row.names=TRUE)
write.table(trainAddMetaData_lvVbH, 
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", 
                       paste0("trainData_justChildSamples_samplesWithDate_splitProp",trainSplitProp,"_justDreden_andVanderbiltLiverpoolHealthy.tsv")), sep="\t", quote=FALSE, row.names=TRUE)



# Plot
ageBySite_plt <- ggplot(trainMetaData, aes(x=Age, fill=Site)) + 
  geom_density(alpha=0.65) + 
  scale_fill_manual(values=siteCols) +
  facet_grid(cols=vars(Condition)) +
  xlab("Age") + ylab("Density") + ggtitle("Patient Age by Site") +
  theme_classic(base_size=24) +
           theme(axis.text    = element_text(size=42),
                 axis.title   = element_text(size=48),
                 plot.title   = element_text(size=60),
                 strip.text   = element_text(size=32),
                 legend.text  = element_text(size=32),
                 legend.title = element_text(size=42))

# Save
ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
                  paste0("ageDensityPlt_colBySite_facetByCondition_trainSubsetData.png")), 
       ageBySite_plt, width=10000, height=4000, units="px")


# Enumerate sex by site by condition
sexFreq <- data.frame(table(trainMetaData[,c("Sex","Condition","Site")]))

# Plot
ageBySite_plt <- ggplot(sexFreq, aes(x=Sex, y=Freq, fill=Site)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=siteCols) +
  facet_grid(cols=vars(Condition)) +
  xlab("Sex") + ylab("Density") + ggtitle("Patient Sex by Site") +
  theme_classic(base_size=24) +
           theme(axis.text    = element_text(size=42),
                 axis.title   = element_text(size=48),
                 plot.title   = element_text(size=60),
                 strip.text   = element_text(size=20),
                 legend.text  = element_text(size=32),
                 legend.title = element_text(size=42))

# Save
ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
                  paste0("sexFreqPlot_colBySite_facetByCondition_trainSubsetData.png")), 
       ageBySite_plt, width=6000, height=3000, units="px")

# Plot
ageBySite_plt <- ggplot(sexFreq, aes(x=Sex, y=Freq, fill=Condition)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=condCols) +
  facet_grid(cols=vars(Site)) +
  xlab("Sex") + ylab("Density") + ggtitle("Patient Sex by Site") +
  theme_classic(base_size=24) +
           theme(axis.text    = element_text(size=42),
                 axis.title   = element_text(size=48),
                 plot.title   = element_text(size=60),
                 strip.text   = element_text(size=20),
                 legend.text  = element_text(size=32),
                 legend.title = element_text(size=42))

# Save
ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
                  paste0("sexFreqPlot_colByCondition_facetBySite_trainSubsetData.png")), 
       ageBySite_plt, width=6000, height=3000, units="px")

# Plot
ethnFreq <- data.frame(table(trainMetaData[,c("Ethnicity","Condition","Site")]))
condByEthn_plt <- ggplot(ethnFreq, aes(x=Ethnicity, y=Freq, fill=Condition)) + 
  geom_bar(stat="identity") + 
  scale_fill_manual(values=condCols) +
  facet_grid(cols=vars(Site)) +
  xlab("Sex") + ylab("Density") + ggtitle("Patient Sex by Site") +
  theme_classic(base_size=24) +
           theme(axis.text    = element_text(size=42),
                 axis.text.x  = element_text(angle=90,
                                             hjust=1,
                                             vjust=0.5,
                                             size=22),
                 axis.title   = element_text(size=48),
                 plot.title   = element_text(size=60),
                 strip.text   = element_text(size=20),
                 legend.text  = element_text(size=32),
                 legend.title = element_text(size=42))

# Save
ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", 
                  paste0("ethnicityFreqPlot_colByCondition_facetBySite_trainSubsetData.png")), 
       condByEthn_plt, width=6000, height=3000, units="px")





# Define alternate model fitting
varPart_altControl <- lme4::lmerControl(optimizer = 'Nelder_Mead', calc.derivs=FALSE)



# If run variance analysis
if (runVarAnalysis) {

  # Loop over data subsets
  dataSubsets <- list(trainMetaData, trainMetaData_noVbilt, trainMetaData_incVblt, trainMetaData_dresden, trainMetaData_noLivpl, trainMetaData_liverpl, trainMetaData_lvplH, trainMetaData_lvVbH)
  names(dataSubsets) <- c("All Sites", "Liverpool And Dresden", "Only Vanderbilt", "Only Dresden", "Dresden And Vanderbilt", "Just Liverpool", "Dresden, Vanderbilt and Healthy Liverpool", "Dresden, Vanderbilt Healthy And Liverpool Healthy")
  for (subSet in names(dataSubsets)) {
    
    print(subSet)
    
    # Get subsetted samples
    meta <- dataSubsets[[subSet]]
    
    # If more than one site
    if (length(unique(meta$Site)) > 1) {
      
      # Include site terms
      crossSectional_formula <- ~ Z_Age +
                                  (1|PLATE) + (1|Sex) + (1|Condition) + (1|Site) +
                                  (1|Condition:Site)
        
      # Define variables to test - with site
      factorsForPVCA <- c("Plate","Sex","Condition","Site","Ethnicity")
    
    # Else exclude site
    } else {
      
      # Include site terms
      crossSectional_formula <- ~ Z_Age +
                                  (1|PLATE) + (1|Sex) + (1|Condition)
      
      # Define variables to test
      factorsForPVCA <- c("Plate","Sex","Condition")
      
    }
    
    # Randomly sample orders
    meta <- meta[!is.na(meta$Age),]
    
    # Subset samples in quant table
    subQuant <- trainQuantTable[,meta$Sample_ID]
    
    # Remove proteins with NAs
    subQuant <- subQuant[!is.na(rowSums(subQuant)),]
    
    # PVCA - try multiple times to ensure consistent results
    varPartRes <- do.call("rbind",lapply(1:10, function(i) {
      
      # Sample metadata 
      meta <- meta[sample(1:nrow(meta), nrow(meta)),]
    
      # Fit variance partition
      vpInteraction <- fitExtractVarPartModel(subQuant[,meta$Sample_ID], crossSectional_formula, meta,
                                              control=varPart_altControl)
      
      # Tidy column names
      colnames(vpInteraction) <- mgsub(colnames(vpInteraction), c("Z_","PLATE"), c("","Plate"))
      
      # Melt
      vpTall <- data.frame(Sample=i,
                           suppressWarnings(reshape2::melt(data.frame(ID=rownames(vpInteraction), vpInteraction), id.vars="ID")))
      colnames(vpTall) <- c("Sample", "ID", "Term", "Variance")
      
      # Return
      return(vpTall)
      
    }))
    
    # Define colours
    plotColors <- hue_pal()(length(unique(varPartRes$Term)))
    names(plotColors) <- unique(varPartRes$Term)
    plotColors["Residuals"] <- "lightgrey"
    
    # Sort term orders
    varPartRes <- varPartRes[order(varPartRes$Term, decreasing=TRUE),]
    varPartRes <- rbind(varPartRes[varPartRes$Term != "Residuals",],
                        varPartRes[varPartRes$Term == "Residuals",])
    
    # Plot
    plt <- ggplot(varPartRes,
                  aes(x=Sample, y=Variance, fill=Term, group=Sample)) +
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_violin(scale = "width") +
                      geom_boxplot(fill="white", width=0.1) +
                      xlab("Iterations") + ylab("Proportion of Variance") +
                      facet_grid(cols=vars(Term), scale="free_x", space="free_x") +
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      theme_classic(base_size=24) +
                               theme(axis.text.x  = element_blank(),
                                     axis.text.y  = element_text(size=36),
                                     axis.title   = element_text(size=48),
                                     legend.text  = element_text(size=24),
                                     legend.title = element_blank(),
                                     plot.title   = element_text(size=60))
    agg_plt <- ggplot(aggregate(. ~ ID + Term, data=varPartRes, FUN=median),
                  aes(x=Term, y=Variance, fill=Term, group=Term)) +
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_violin(scale = "width") +
                      geom_boxplot(fill="white", width=0.1) +
                      xlab("Terms") + ylab("Proportion of Variance") +
                      facet_grid(cols=vars(Term), scale="free_x", space="free_x") +
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      theme_classic(base_size=24) +
                               theme(axis.text.y  = element_text(size=36),
                                     axis.title   = element_text(size=48),
                                     legend.text  = element_text(size=24),
                                     legend.title = element_blank(),
                                     plot.title   = element_text(size=60))
    
    
    # Save
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "varPartition", 
                      paste0("varPartPlt_allIterations_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           plt, width=10000, height=4000, units="px")
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "varPartition", 
                      paste0("varPartPlt_medianAcrossIterations_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           agg_plt, width=10000, height=4000, units="px")
    
    # Fix wonky names
    colnames(meta) <- gsub("PLATE","Plate",colnames(meta))
    
    # PVCA - try multiple times to ensure consistent results
    pvcaRes <- do.call("rbind",lapply(1:10, function(i) {
      
      # Randomly shuffle sample IDs
      meta <- meta[sample(1:nrow(meta), nrow(meta)),]
      
      # Run PVCA
      out <- data.frame(Sample=i, suppressWarnings(performPVCA(m=subQuant[,meta$Sample_ID], meta=meta, batchFactors=factorsForPVCA)))
      return(out)
      
    }))
    
    # Set consistent order
    pvcaRes <- pvcaRes[order(pvcaRes$Effect, decreasing=FALSE),]
    
    # Put residuals at the end
    pvcaRes <- rbind(pvcaRes[pvcaRes$Effect != "resid",],
                     pvcaRes[pvcaRes$Effect == "resid",])
    pvcaRes$Effect <- factor(pvcaRes$Effect, levels=unique(pvcaRes$Effect))
    
    # Define colours
    plotColors <- hue_pal()(length(unique(pvcaRes$Effect)))
    names(plotColors) <- unique(pvcaRes$Effect)
    plotColors["resid"] <- "lightgrey"
    
    # Plot
    box_plt <- ggplot(pvcaRes, 
                      aes(x=Effect, y=PercVar, fill=Effect)) + 
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_boxplot() + 
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      xlab("") + ylab("% Variance\nExplained") + ylim(0,100) +
                           theme_classic(base_size=24) +
                                    theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=32),
                                          axis.text.y  = element_text(size=36),
                                          axis.title   = element_text(size=48),
                                          legend.text  = element_text(size=24),
                                          legend.title = element_blank(),
                                          plot.title   = element_text(size=42))
    bar_plt <- ggplot(aggregate(. ~ Effect + Method, data=pvcaRes, FUN=median),
                      aes(x=Effect, y=PercVar, fill=Effect, label=paste0(round(PercVar, 2), "%"))) + 
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_bar(stat="identity", color="black") + 
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      xlab("") + ylab("% Variance\nExplained") + ylim(0,100) +
                      geom_text(nudge_y=3, size=6) +
                           theme_classic(base_size=24) +
                                    theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=32),
                                          axis.text.y  = element_text(size=36),
                                          axis.title   = element_text(size=48),
                                          legend.text  = element_text(size=24),
                                          legend.title = element_blank(),
                                          plot.title   = element_text(size=42))
    
    # Save
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "pvca", 
                      paste0("pvcaBoxPlot_allIterations_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           box_plt, width=6000, height=4000, units="px")
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "pvca", 
                      paste0("pvcaBarPlot_medianAcrossIterations_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           bar_plt, width=6000, height=4000, units="px")
    
  }

}



# If run variance analysis
if (runVarAnalysis) {

  # Create outputs
  dir.create(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", "varPartition"), showWarnings=FALSE, recursive=TRUE)
  dir.create(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", "pvca"), showWarnings=FALSE, recursive=TRUE)
  
  # Loop over data subsets
  dataSubsets <- list(trainAddMetaData, trainAddMetaData_drsdn, trainAddMetaData_drsVb, trainAddMetaData_drsLp, trainAddMetaData_lvplH, trainAddMetaData_lvVbH)
  names(dataSubsets) <- c("All Sites", "Only Dresden", "Only Dresden And Vanderbilt", "Only Dresden And Liverpool", "Dresden, Vanderbilt and Healthy Liverpool", "Dresden, Vanderbilt Healthy And Liverpool Healthy")
  for (subSet in names(dataSubsets)) {
    
    print(subSet)
    
    # Get subsetted samples
    meta <- dataSubsets[[subSet]]
    
    # If more than one site
    if (length(unique(meta$Site)) > 1) {
      
      # Include site terms
      crossSectional_formula <- ~ Z_Age +
                                  (1|PLATE) + (1|Sex) + (1|Condition) + (1|Site) +
                                  (1|Site:Condition)
      
      # Define variables to test - with site
      factorsForPVCA <- c("Plate","Sex","Condition","Site","Ethnicity")
    
    # Else exclude site
    } else {
      
      # Include site terms
      crossSectional_formula <- ~ Z_Age +
                                  (1|PLATE) + (1|Sex) + (1|Condition)
      
      # Define variables to test
      factorsForPVCA <- c("Plate","Sex","Condition")
      
    }
    
    # Randomly sample orders
    meta <- meta[!is.na(meta$Age) & !is.na(meta$DateOfSample_asNumber),]
    
    # Subset samples in quant table
    subQuant <- trainQuantTable[,meta$Sample_ID]
    
    # Remove proteins with NAs
    subQuant <- subQuant[!is.na(rowSums(subQuant)),]
    
    # PVCA - try multiple times to ensure consistent results
    varPartRes <- do.call("rbind",lapply(1:10, function(i) {
      
      # Sample metadata 
      meta <- meta[sample(1:nrow(meta), nrow(meta)),]
    
      # Fit variance partition
      vpInteraction <- suppressMessages(suppressWarnings(fitExtractVarPartModel(subQuant[,meta$Sample_ID], crossSectional_formula, meta,
                                                                                control=varPart_altControl)))
      
      # Tidy column names
      colnames(vpInteraction) <- mgsub(colnames(vpInteraction), c("Z_","PLATE"), c("","Plate"))
      
      # Melt
      vpTall <- data.frame(Sample=i,
                           suppressWarnings(reshape2::melt(data.frame(ID=rownames(vpInteraction), vpInteraction), id.vars="ID")))
      colnames(vpTall) <- c("Sample", "ID", "Term", "Variance")
      
      # Return
      return(vpTall)
      
    }))
    
    # Define colours
    plotColors <- hue_pal()(length(unique(varPartRes$Term)))
    names(plotColors) <- unique(varPartRes$Term)
    plotColors["Residuals"] <- "lightgrey"
    
    # Sort term orders
    varPartRes <- varPartRes[order(varPartRes$Term, decreasing=TRUE),]
    varPartRes <- rbind(varPartRes[varPartRes$Term != "Residuals",],
                        varPartRes[varPartRes$Term == "Residuals",])
    
    # Plot
    plt <- ggplot(varPartRes,
                  aes(x=Sample, y=Variance, fill=Term, group=Sample)) +
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_violin(scale = "width") +
                      geom_boxplot(fill="white", width=0.1) +
                      xlab("Iterations") + ylab("Proportion of Variance") +
                      facet_grid(cols=vars(Term), scale="free_x", space="free_x") +
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      theme_classic(base_size=24) +
                               theme(axis.text.x  = element_blank(),
                                     axis.text.y  = element_text(size=36),
                                     axis.title   = element_text(size=48),
                                     legend.text  = element_text(size=24),
                                     legend.title = element_blank(),
                                     plot.title   = element_text(size=42))
    agg_plt <- ggplot(aggregate(. ~ ID + Term, data=varPartRes, FUN=median),
                  aes(x=Term, y=Variance, fill=Term, group=Term)) +
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_violin(scale = "width") +
                      geom_boxplot(fill="white", width=0.1) +
                      xlab("Terms") + ylab("Proportion of Variance") +
                      facet_grid(cols=vars(Term), scale="free_x", space="free_x") +
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      theme_classic(base_size=24) +
                               theme(axis.text.y  = element_text(size=36),
                                     axis.title   = element_text(size=48),
                                     legend.text  = element_text(size=24),
                                     legend.title = element_blank(),
                                     plot.title   = element_text(size=42))
    
    
    # Save
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), 
                      "sampleRandom", "samplesWithDate", "varPartition", paste0("varPartPlt_allIterations_excludeDate_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           plt, width=10000, height=4000, units="px")
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), 
                      "sampleRandom", "samplesWithDate", "varPartition", paste0("varPartPlt_medianAcrossIterations_excludeDate_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           agg_plt, width=10000, height=4000, units="px")
    
    # Fix wonky names
    colnames(meta) <- gsub("PLATE","Plate",colnames(meta))
    
    # PVCA - try multiple times to ensure consistent results
    pvcaRes <- do.call("rbind",lapply(1:10, function(i) {
      
      # Randomly shuffle sample IDs
      meta <- meta[sample(1:nrow(meta), nrow(meta)),]
      
      # Run PVCA
      out <- data.frame(Sample=i, suppressMessages(suppressWarnings(performPVCA(m=subQuant[,meta$Sample_ID], meta=meta, batchFactors=factorsForPVCA))))
      return(out)
      
    }))
    
    # Set consistent order
    pvcaRes <- pvcaRes[order(pvcaRes$Effect, decreasing=FALSE),]
    
    # Put residuals at the end
    pvcaRes <- rbind(pvcaRes[pvcaRes$Effect != "resid",],
                     pvcaRes[pvcaRes$Effect == "resid",])
    pvcaRes$Effect <- factor(pvcaRes$Effect, levels=unique(pvcaRes$Effect))
    
    # Define colours
    plotColors <- hue_pal()(length(unique(pvcaRes$Effect)))
    names(plotColors) <- unique(pvcaRes$Effect)
    plotColors["resid"] <- "lightgrey"
    
    # Plot
    box_plt <- ggplot(pvcaRes, 
                      aes(x=Effect, y=PercVar, fill=Effect)) + 
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_boxplot() + 
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      xlab("") + ylab("% Variance\nExplained") + ylim(0,100) +
                           theme_classic(base_size=24) +
                                    theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=32),
                                          axis.text.y  = element_text(size=36),
                                          axis.title   = element_text(size=48),
                                          legend.text  = element_text(size=16),
                                          legend.title = element_blank(),
                                          plot.title   = element_text(size=48)) +
                      guides(fill=guide_legend(ncol=1))
    bar_plt <- ggplot(aggregate(. ~ Effect + Method, data=pvcaRes, FUN=median),
                      aes(x=Effect, y=PercVar, fill=Effect, label=paste0(round(PercVar, 2), "%"))) + 
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_bar(stat="identity", color="black") + 
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      xlab("") + ylab("% Variance\nExplained") + ylim(0,100) +
                      geom_text(nudge_y=3, size=4) +
                           theme_classic(base_size=24) +
                                    theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=32),
                                          axis.text.y  = element_text(size=36),
                                          axis.title   = element_text(size=48),
                                          legend.text  = element_text(size=16),
                                          legend.title = element_blank(),
                                          plot.title   = element_text(size=48)) +
                      guides(fill=guide_legend(ncol=1))
    
    # Save
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), 
                      "sampleRandom", "samplesWithDate", "pvca", paste0("pvcaBoxPlot_allIterations_excludeDate_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           box_plt, width=6000, height=4000, units="px")
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), 
                      "sampleRandom", "samplesWithDate", "pvca", paste0("pvcaBarPlot_medianAcrossIterations_excludeDate_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           bar_plt, width=6000, height=4000, units="px")
    
  }

}



# If run variance analysis
if (runVarAnalysis) {
  
  # Create outputs
  dir.create(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", "varPartition"), showWarnings=FALSE, recursive=TRUE)
  dir.create(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", "pvca"), showWarnings=FALSE, recursive=TRUE)
  
  # Loop over data subsets
  dataSubsets <- list(trainAddMetaData, trainAddMetaData_drsdn, trainAddMetaData_drsVb, trainAddMetaData_drsLp, trainAddMetaData_lvplH, trainAddMetaData_lvVbH)
  names(dataSubsets) <- c("All Sites", "Only Dresden", "Only Dresden And Vanderbilt", "Only Dresden And Liverpool", "Dresden, Vanderbilt and Healthy Liverpool", "Dresden, Vanderbilt Healthy And Liverpool Healthy")
  for (subSet in names(dataSubsets)[-5]) {
    
    print(subSet)
    
    # Get subsetted samples
    meta <- dataSubsets[[subSet]]
      
    # Rename date as categocial
    colnames(meta)[grepl("Date_asCategorical",colnames(meta))] <- "Date"
    
    # If more than one site
    if (length(unique(meta$Site)) > 1) {
      
      # Include site terms
      crossSectional_formula <- ~ Z_Age + Z_Date +
                                  (1|PLATE) + (1|Sex) + (1|Condition) + (1|Site) +
                                  (1|Site:Condition) + (1|Condition:Z_Date)
      
      # Define variables to test - with site
      factorsForPVCA <- c("Plate","Sex","Condition","Site","Ethnicity","Date")
    
    # Else exclude site
    } else {
      
      # Include site terms
      crossSectional_formula <- ~ Z_Age + Z_Date +
                                  (1|PLATE) + (1|Sex) + (1|Condition) + 
                                  (1|Condition:Z_Date)
      
      # Define variables to test
      factorsForPVCA <- c("Plate","Sex","Condition","Date")
      
    }
    
    # Randomly sample orders
    meta <- meta[!is.na(meta$Age) & !is.na(meta$Date),]
    
    # Subset samples in quant table
    subQuant <- trainQuantTable[,meta$Sample_ID]
    
    # Remove proteins with NAs
    subQuant <- subQuant[!is.na(rowSums(subQuant)),]
    
    # PVCA - try multiple times to ensure consistent results
    varPartRes <- do.call("rbind",lapply(1:10, function(i) {
      
      # Sample metadata 
      meta <- meta[sample(1:nrow(meta), nrow(meta)),]
    
      # Fit variance partition
      vpInteraction <- suppressMessages(suppressWarnings(fitExtractVarPartModel(subQuant[,meta$Sample_ID], crossSectional_formula, meta,
                                                                                control=varPart_altControl)))
      
      # Tidy column names
      colnames(vpInteraction) <- mgsub(colnames(vpInteraction), c("Z_","PLATE"), c("","Plate"))
      
      # Melt
      vpTall <- data.frame(Sample=i,
                           suppressWarnings(reshape2::melt(data.frame(ID=rownames(vpInteraction), vpInteraction), id.vars="ID")))
      colnames(vpTall) <- c("Sample", "ID", "Term", "Variance")
      
      # Return
      return(vpTall)
      
    }))
    
    # Define colours
    plotColors <- hue_pal()(length(unique(varPartRes$Term)))
    names(plotColors) <- unique(varPartRes$Term)
    plotColors["Residuals"] <- "lightgrey"
    
    # Sort term orders
    varPartRes <- varPartRes[order(varPartRes$Term, decreasing=TRUE),]
    varPartRes <- rbind(varPartRes[varPartRes$Term != "Residuals",],
                        varPartRes[varPartRes$Term == "Residuals",])
    
    # Plot
    plt <- ggplot(varPartRes,
                  aes(x=Sample, y=Variance, fill=Term, group=Sample)) +
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_violin(scale = "width") +
                      geom_boxplot(fill="white", width=0.1) +
                      xlab("Iterations") + ylab("Proportion of Variance") +
                      facet_grid(cols=vars(Term), scale="free_x", space="free_x") +
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      theme_classic(base_size=24) +
                               theme(axis.text.x  = element_blank(),
                                     axis.text.y  = element_text(size=36),
                                     axis.title   = element_text(size=48),
                                     legend.text  = element_text(size=24),
                                     legend.title = element_blank(),
                                     plot.title   = element_text(size=60))
    agg_plt <- ggplot(aggregate(. ~ ID + Term, data=varPartRes, FUN=median),
                  aes(x=Term, y=Variance, fill=Term, group=Term)) +
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_violin(scale = "width") +
                      geom_boxplot(fill="white", width=0.1) +
                      xlab("Terms") + ylab("Proportion of Variance") +
                      facet_grid(cols=vars(Term), scale="free_x", space="free_x") +
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      theme_classic(base_size=24) +
                               theme(axis.text.y  = element_text(size=36),
                                     axis.title   = element_text(size=48),
                                     legend.text  = element_text(size=24),
                                     legend.title = element_blank(),
                                     plot.title   = element_text(size=60))
    
    
    # Save
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", "varPartition", 
                      paste0("varPartPlt_allIterations_withDate_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           plt, width=10000, height=4000, units="px")
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", "varPartition", 
                      paste0("varPartPlt_medianAcrossIterations_withDate_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           agg_plt, width=10000, height=4000, units="px")
    
    # Fix wonky names
    colnames(meta) <- gsub("PLATE","Plate",colnames(meta))
    
    # PVCA - try multiple times to ensure consistent results
    pvcaRes <- do.call("rbind",lapply(1:10, function(i) {
      
      # Randomly shuffle sample IDs
      meta <- meta[sample(1:nrow(meta), nrow(meta)),]
      
      # Run PVCA
      out <- data.frame(Sample=i, suppressMessages(suppressWarnings(performPVCA(m=subQuant[,meta$Sample_ID], meta=meta, batchFactors=factorsForPVCA))))
      return(out)
      
    }))
    
    # Set consistent order
    pvcaRes <- pvcaRes[order(pvcaRes$Effect, decreasing=FALSE),]
    
    # Put residuals at the end
    pvcaRes <- rbind(pvcaRes[pvcaRes$Effect != "resid",],
                     pvcaRes[pvcaRes$Effect == "resid",])
    pvcaRes$Effect <- factor(pvcaRes$Effect, levels=unique(pvcaRes$Effect))
    
    # Define colours
    plotColors <- hue_pal()(length(unique(pvcaRes$Effect)))
    names(plotColors) <- unique(pvcaRes$Effect)
    plotColors["resid"] <- "lightgrey"
    
    # Plot
    box_plt <- ggplot(pvcaRes, 
                      aes(x=Effect, y=PercVar, fill=Effect)) + 
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_boxplot() + 
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      xlab("") + ylab("% Variance\nExplained") + ylim(0,100) +
                           theme_classic(base_size=24) +
                                    theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=32),
                                          axis.text.y  = element_text(size=36),
                                          axis.title   = element_text(size=48),
                                          legend.text  = element_text(size=16),
                                          legend.title = element_blank(),
                                          plot.title   = element_text(size=48)) +
                      guides(fill=guide_legend(ncol=1))
    bar_plt <- ggplot(aggregate(. ~ Effect + Method, data=pvcaRes, FUN=median),
                      aes(x=Effect, y=PercVar, fill=Effect, label=paste0(round(PercVar, 2), "%"))) + 
                      geom_hline(yintercept=0,  color="darkgrey") +
                      geom_bar(stat="identity", color="black") + 
                      ggtitle(paste0(subSet, " Child Samples")) + 
                      scale_fill_manual(values=plotColors) +
                      xlab("") + ylab("% Variance\nExplained") + ylim(0,100) +
                      geom_text(nudge_y=3, size=4) +
                           theme_classic(base_size=24) +
                                    theme(axis.text.x  = element_text(angle=90, hjust=1, vjust=0.5, size=32),
                                          axis.text.y  = element_text(size=36),
                                          axis.title   = element_text(size=48),
                                          legend.text  = element_text(size=16),
                                          legend.title = element_blank(),
                                          plot.title   = element_text(size=48)) +
                      guides(fill=guide_legend(ncol=1))
    
    # Save
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", "pvca", 
                      paste0("pvcaBoxPlot_allIterations_withDate_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           box_plt, width=6000, height=4000, units="px")
    ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "samplesWithDate", "pvca", 
                      paste0("pvcaBarPlot_medianAcrossIterations_withDate_",gsub(" |[,]","",subSet), "_trainSubsetData.png")), 
           bar_plt, width=6000, height=4000, units="px")
    
  }

}
  


# Check all condiitions
checkAllConditions <- all(c("CNO", "Healthy", "Crohn", "Oncology") %in% unique(trainMetaData$Condition))

# If not just healthy
if (checkAllConditions) {
  
  # Fix Site and condition
  trainMetaData_dresden$Condition <- as.factor(as.character(trainMetaData_dresden$Condition))
  
  # Test model AICs
  modelAICs1 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainMetaData_dresden$Sample_ID],
                              meta=trainMetaData_dresden,
                              toDrop=c("Sex","Condition:Z_Age","Condition:PLATE","Condition:Sex","PLATE","Z_Age"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age","Condition:PLATE","Condition:Sex"))
  
  # Test model AICs - drop condition:plate
  modelAICs2 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainMetaData_dresden$Sample_ID],
                              meta=trainMetaData_dresden,
                              toDrop=c("Sex","Condition:Z_Age","Condition:Sex","PLATE","Z_Age"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age","Condition:Sex"))
  
  # Test model AICs - drop condition:sex
  modelAICs3 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainMetaData_dresden$Sample_ID],
                              meta=trainMetaData_dresden,
                              toDrop=c("Sex","Condition:Z_Age","PLATE","Z_Age"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age"))
  
  # Test model AICs - drop condition:age
  modelAICs4 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainMetaData_dresden$Sample_ID],
                              meta=trainMetaData_dresden,
                              toDrop=c("Sex","PLATE","Z_Age"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex"))
  
  # Drop full for each
  modelAICs1$IC <- modelAICs1$IC[,-1]
  modelAICs2$IC <- modelAICs2$IC[,-1]
  modelAICs3$IC <- modelAICs3$IC[,-1]
  modelAICs4$IC <- modelAICs4$IC[,-1]
  
  # Combine to big list 
  allAICs <- list.append(list(modelAICs1), modelAICs2, modelAICs3, modelAICs4) #, , modelAICs5, modelAICs6)
  
  # Combine
  dropped <- c("None","Condition:Plate","Condition:Sex","Condition:Age", "Condition:Date", "Date")
  all <- lapply(seq_len(length(allAICs)), function(i) data.frame(ModelSet=dropped[[i]], Protein=rownames(allAICs[[i]]$IC), allAICs[[i]]$IC - allAICs[[1]]$IC[,1]))
  allTall <- do.call("rbind",lapply(all, function(model) reshape2::melt(model, id.vars=c("ModelSet","Protein"))))
  
  # Colour
  allTall$variable <- as.character(allTall$variable)
  allTall["DroppedTerm"] <- "No"
  allTall[allTall$variable == "Condition.PLATE",]["DroppedTerm"] <- "Yes"
  allTall[allTall$variable == "Condition.Sex" & allTall$ModelSet == "Condition:Plate",]["DroppedTerm"] <- "Yes"
  allTall[(allTall$variable == "Condition.Z_Age" & allTall$ModelSet == "Condition:Sex") | (allTall$variable == "Condition.Z_Age" & allTall$ModelSet == "Condition:Sex"),]["DroppedTerm"] <- "Yes"
  
  # Calculate median
  medAll <- aggregate(. ~ ModelSet + variable + DroppedTerm, data=allTall[,!grepl("Protein",colnames(allTall))], FUN=median)
  
  # Replace "None" with all
  allTall[allTall$ModelSet == "None",]["ModelSet"] <- "All Terms"
  medAll[medAll$ModelSet == "None",]["ModelSet"] <- "All Terms"
  
  # Plot
  allTall$ModelSet <- factor(allTall$ModelSet, levels=c("All Terms","Condition:Plate","Condition:Sex","Condition:Age", "Condition:Date", "Date"))
  p1 <- ggplot(allTall, aes(x=gsub("Z_","",variable), y=value, fill=DroppedTerm)) + 
    geom_hline(yintercept=0,  color="darkgrey") +
    geom_point(size=8, alpha=0.45, shape=21) + 
    stat_summary(fun.y = median, geom = "errorbar", 
                 aes(ymax = ..y.., ymin = ..y.., group = variable),
                 width = 1, linetype = "solid") +
    xlab("Variable being Dropped") + 
    ylab("Difference in AIC") +
    scale_fill_manual(values=c("darkgrey","#99E6FF")) +
    facet_grid(cols=vars(ModelSet), scale="free", space="free") +
    theme_classic(base_size=48) + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          axis.title=element_text(size=64),
          strip.text=element_text(size=36),
          legend.position="top",
          legend.direction="horizontal") +
    guides(fill=guide_legend(override.aes=list(size=14,
                                               alpha=0.85)))
  p2 <- ggplot(allTall, aes(x=variable, y=value, fill=DroppedTerm)) + 
    geom_hline(yintercept=0,  color="darkgrey") +
    geom_boxplot() +
    xlab("Variable being Dropped") +
    ylab("Difference in AIC") +
    scale_fill_manual(values=c("darkgrey","#99E6FF")) +
    facet_grid(cols=vars(ModelSet), scale="free", space="free") +
    theme_classic(base_size=48) + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          axis.title=element_text(size=64),
          strip.text=element_text(size=36),
          legend.position="top",
          legend.direction="horizontal")
  print(p1)
  print(p2)
  
  # Save
  ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom",
                    "modelOptimisation_termDrop_differenceInAIC.png"),
         p1,
         units="px", width=8000, height=4000)
  
}



# If a subset
if (!checkAllConditions) {
  
  # Fix Site and condition
  trainMetaData$Condition <- as.factor(as.character(trainMetaData$Condition))
  
  # Test model AICs
  modelAICs1 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),
                                               trainMetaData$Sample_ID],
                              meta=trainMetaData,
                              toDrop=c("Sex","Condition:Z_Age","Condition:PLATE","Condition:Sex","PLATE","Z_Age"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age","Condition:PLATE","Condition:Sex"))
    
  # Plot
  plotAIC(modelAICs1,
          dropTerm="Condition:Sex")
  
  # Test model AICs
  modelAICs2 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),
                                               trainMetaData$Sample_ID],
                              meta=trainMetaData,
                              toDrop=c("Sex","Condition:Z_Age","Condition:PLATE","PLATE","Z_Age"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age","Condition:PLATE"))
    
  # Plot
  plotAIC(modelAICs2,
          dropTerm="Condition:PLATE")
  
  # Test model AICs
  modelAICs3 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),
                                               trainMetaData$Sample_ID],
                              meta=trainMetaData,
                              toDrop=c("Sex","Condition:Z_Age","PLATE","Z_Age"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age"))
    
  # Plot
  plotAIC(modelAICs3)
  
}



# Add to list
newList <- list(p1)
names(newList) <- c("modelAICOpt")
plotList <- append(plotList,
                   newList)



# If not just healthy
if (checkAllConditions) {
  
  # Fix Site and condition
  trainAddMetaData_drsdn$Condition <- as.factor(as.character(trainAddMetaData_drsdn$Condition))
  
  # Test model AICs
  modelAICs1 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainAddMetaData_drsdn$Sample_ID],
                              meta=trainAddMetaData_drsdn,
                              toDrop=c("Sex","Condition:Z_Age","Condition:PLATE","Condition:Sex","PLATE","Z_Age","Z_Date","Condition:Z_Date"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Z_Date","Condition:Z_Age","Condition:PLATE","Condition:Sex","Condition:Z_Date"))
  
  # Test model AICs - drop condition:plate
  modelAICs2 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainAddMetaData_drsdn$Sample_ID],
                              meta=trainAddMetaData_drsdn,
                              toDrop=c("Sex","Condition:Z_Age","Condition:Sex","PLATE","Z_Age","Z_Date","Condition:Z_Date"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Z_Date","Condition:Z_Age","Condition:Sex","Condition:Z_Date"))
  
  # Test model AICs - drop condition:sex
  modelAICs3 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainAddMetaData_drsdn$Sample_ID],
                              meta=trainAddMetaData_drsdn,
                              toDrop=c("Sex","Condition:Z_Age","PLATE","Z_Age","Condition:Z_Date"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age","Z_Date","Condition:Z_Date"))
  
  # Test model AICs - drop condition:age
  modelAICs4 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainAddMetaData_drsdn$Sample_ID],
                              meta=trainAddMetaData_drsdn,
                              toDrop=c("Sex","PLATE","Z_Age","Z_Date","Condition:Z_Date"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Z_Date","Condition:Z_Date"))
  
  # Test model AICs - drop condition:date
  modelAICs5 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainAddMetaData_drsdn$Sample_ID],
                              meta=trainAddMetaData_drsdn,
                              toDrop=c("Sex","PLATE","Z_Age","Z_Date"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Z_Date"))
  
  # Test model AICs - drop date
  modelAICs6 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),trainAddMetaData_drsdn$Sample_ID],
                              meta=trainAddMetaData_drsdn,
                              toDrop=c("Sex","PLATE","Z_Age"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex"))
  
  # Drop full for each
  modelAICs1$IC <- modelAICs1$IC[,-1]
  modelAICs2$IC <- modelAICs2$IC[,-1]
  modelAICs3$IC <- modelAICs3$IC[,-1]
  modelAICs4$IC <- modelAICs4$IC[,-1]
  modelAICs5$IC <- modelAICs5$IC[,-1]
  modelAICs6$IC <- modelAICs6$IC[,-1]
  
  # Combine to big list
  allAICs <- list.append(list(modelAICs1), modelAICs2, modelAICs3, modelAICs4, modelAICs5, modelAICs6)
  
  # Combine
  dropped <- c("None","Condition:Plate","Condition:Sex","Condition:Age", "Condition:Date", "Date")
  all <- lapply(seq_len(length(allAICs)), function(i) data.frame(ModelSet=dropped[[i]], Protein=rownames(allAICs[[i]]$IC), allAICs[[i]]$IC - allAICs[[1]]$IC[,1]))
  allTall <- do.call("rbind",lapply(all, function(model) reshape2::melt(model, id.vars=c("ModelSet","Protein"))))
  
  # Colour
  allTall$variable <- as.character(allTall$variable)
  allTall["DroppedTerm"] <- "No"
  allTall[allTall$variable == "Condition.PLATE",]["DroppedTerm"] <- "Yes"
  allTall[allTall$variable == "Condition.Sex" & allTall$ModelSet == "Condition:Plate",]["DroppedTerm"] <- "Yes"
  allTall[(allTall$variable == "Condition.Z_Age" & allTall$ModelSet == "Condition:Sex") | (allTall$variable == "Condition.Z_Age" & allTall$ModelSet == "Condition:Sex"),]["DroppedTerm"] <- "Yes"
  allTall[(allTall$variable == "Condition.Z_Date" & allTall$ModelSet == "Condition:Age") | (allTall$variable == "Condition.Z_Date" & allTall$ModelSet == "Condition:Age"),]["DroppedTerm"] <- "Yes"
  allTall[(allTall$variable == "Z_Date" & allTall$ModelSet == "Condition:Date") | (allTall$variable == "Z_Date" & allTall$ModelSet == "Condition:Date"),]["DroppedTerm"] <- "Yes"
  
  # Calculate median
  medAll <- aggregate(. ~ ModelSet + variable + DroppedTerm, data=allTall[,!grepl("Protein",colnames(allTall))], FUN=median)
  
  # Replace "None" with all
  allTall[allTall$ModelSet == "None",]["ModelSet"] <- "All Terms"
  medAll[medAll$ModelSet == "None",]["ModelSet"] <- "All Terms"
  
  # Plot
  allTall$ModelSet <- factor(allTall$ModelSet, levels=c("All Terms","Condition:Plate","Condition:Sex","Condition:Age", "Condition:Date", "Date"))
  p1 <- ggplot(allTall, aes(x=gsub("Z_","",variable), y=value, fill=DroppedTerm)) + 
    geom_hline(yintercept=0,  color="darkgrey") +
    geom_point(size=8, alpha=0.15, shape=21) + 
    stat_summary(fun.y = median, geom = "errorbar", 
                 aes(ymax = ..y.., ymin = ..y.., group = variable),
                 width = 1, linetype = "solid") +
    xlab("Variable being Dropped") + 
    ylab("Difference in AIC") +
    scale_fill_manual(values=c("darkgrey","#99E6FF")) +
    facet_grid(cols=vars(ModelSet), scale="free", space="free") +
    theme_classic(base_size=48) + 
    theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
          axis.title=element_text(size=64),
          strip.text=element_text(size=28),
          legend.position="top",
          legend.direction="horizontal") +
    guides(fill=guide_legend(override.aes=list(size=14,
                                               alpha=0.85)))
  p2 <- ggplot(allTall, aes(x=variable, y=value, fill=DroppedTerm)) + 
    geom_hline(yintercept=0,  color="darkgrey") +
    geom_boxplot() +
    xlab("Variable being Dropped") + ylab("AIC") +
    scale_fill_manual(values=c("darkgrey","#99E6FF")) +
    facet_grid(cols=vars(ModelSet), scale="free", space="free") +
    theme_classic(base_size=16) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  print(p1)
  print(p2)
  
  # Save
  ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom",
                    "modelOptimisation_termDrop_differenceInAIC_samplesWithDate.png"),
         p1,
         units="px", width=8000, height=4000)

  # Define formulas for models
  childCross_modelFormulas <- list(childCrossFormulas_condition=" ~ Condition + Z_Age + PLATE + Sex",
                                   childCrossFormulas_site=" ~ Site + Z_Age + PLATE + Sex",
                                   childCrossFormula_condition_noIntercept=" ~ 0 + Condition + Z_Age + Sex + PLATE",
                                   childCrossFormulas_withDate=" ~ Condition + Z_Age + PLATE + Sex + Z_Date", #  + Condition:Z_Date
                                   childCrossFormulas_withDate_noIntercept=" ~ 0 + Condition + Z_Age + PLATE + Sex + Z_Date") #  + Condition:Z_Date
  
}



# # If a subset
if (!checkAllConditions) {

  # Fix Site and condition
  trainAddMetaData$Condition <- as.factor(as.character(trainAddMetaData$Condition))

  # Test model AICs
  modelAICs1 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),
                                               trainAddMetaData$Sample_ID],
                              meta=trainAddMetaData,
                              toDrop=c("Sex","Condition:Z_Age","Condition:PLATE","Condition:Sex","PLATE","Z_Age","Z_Date","Condition:Z_Date"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age","Condition:PLATE","Condition:Sex","Z_Date","Condition:Z_Date"))

  # Plot
  plotAIC(modelAICs1,
          dropTerm="Condition:Sex")

  # Test model AICs
  modelAICs2 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),
                                               trainAddMetaData$Sample_ID],
                              meta=trainAddMetaData,
                              toDrop=c("Sex","Condition:Z_Age","Condition:PLATE","PLATE","Z_Age","Z_Date","Condition:Z_Date"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age","Condition:PLATE","Z_Date","Condition:Z_Date"))

  # Plot
  plotAIC(modelAICs2,
          dropTerm="Condition:PLATE")

  # Test model AICs
  modelAICs3 <- TestModelAICs(m=normQuantTable[!is.na(rowSums(normQuantTable)),
                                               trainAddMetaData$Sample_ID],
                              meta=trainAddMetaData,
                              toDrop=c("Sex","Condition:Z_Age","PLATE","Z_Age","Z_Date","Condition:Z_Date"),
                              fullModel=c("Condition","Z_Age","PLATE","Sex","Condition:Z_Age","Z_Date","Condition:Z_Date"))

  # Plot
  plotAIC(modelAICs3)

  # Define formulas for models
  childCross_modelFormulas <- list(childCrossFormulas_condition=" ~ Condition + Z_Age + PLATE + Sex",
                                   childCrossFormulas_site=" ~ Site + Z_Age + PLATE + Sex",
                                   childCrossFormula_condition_noIntercept=" ~ 0 + Condition  + Z_Age + Sex + PLATE",
                                   childCrossFormulas_withDate=" ~ Condition + Z_Age + PLATE + Sex + Z_Date", #  + Condition:Z_Date
                                   childCrossFormulas_withDate_noIntercept=" ~ 0 + Condition + Z_Age + PLATE + Sex + Z_Date") #  + Condition:Z_Date

}
  


# Add to list
newList <- list(p1)
names(newList) <- c("modelAICOpt_incDate")
plotList <- append(plotList,
                   newList)



# Make Healthy the reference
## Not considering date
trainMetaData$Condition <- as.factor(as.character(trainMetaData$Condition))
trainMetaData_noVbilt$Condition <- as.factor(as.character(trainMetaData_noVbilt$Condition))
trainMetaData_incVblt$Condition <- as.factor(as.character(trainMetaData_incVblt$Condition))
trainMetaData_dresden$Condition <- as.factor(as.character(trainMetaData_dresden$Condition))
trainMetaData_noLivpl$Condition <- as.factor(as.character(trainMetaData_noLivpl$Condition))
trainMetaData_liverpl$Condition <- as.factor(as.character(trainMetaData_liverpl$Condition))
trainMetaData_lvplH$Condition <- as.factor(as.character(trainMetaData_lvplH$Condition))
trainMetaData_lvVbH$Condition <- as.factor(as.character(trainMetaData_lvVbH$Condition))
if (length(excludeConditions) == 0 | is.null(length(excludeConditions))) { trainMetaData$Condition <- relevel(trainMetaData$Condition, "JIA") }
trainMetaData_noVbilt$Condition <- relevel(trainMetaData_noVbilt$Condition, "Healthy")
trainMetaData_incVblt$Condition <- relevel(trainMetaData_incVblt$Condition, "Healthy")
if (length(excludeConditions) == 0 | is.null(length(excludeConditions))) { trainMetaData_dresden$Condition <- relevel(trainMetaData_dresden$Condition, "JIA") }
trainMetaData_noLivpl$Condition <- relevel(trainMetaData_noLivpl$Condition, "Healthy")
trainMetaData_liverpl$Condition <- relevel(trainMetaData_liverpl$Condition, "Healthy")
trainMetaData_lvplH$Condition <- relevel(trainMetaData_lvplH$Condition, "Healthy")
trainMetaData_lvVbH$Condition <- relevel(trainMetaData_lvVbH$Condition, "Healthy")

## With date
trainAddMetaData$Condition <- relevel(as.factor(as.character(trainAddMetaData$Condition)), "Healthy")
if (length(excludeConditions) == 0 | is.null(length(excludeConditions))) { trainAddMetaData_drsdn$Condition <- relevel(as.factor(as.character(trainAddMetaData_drsdn$Condition)), "JIA") }
trainAddMetaData_drsVb$Condition <- relevel(as.factor(as.character(trainAddMetaData_drsVb$Condition)), "Healthy")
trainAddMetaData_drsLp$Condition <- relevel(as.factor(as.character(trainAddMetaData_drsLp$Condition)), "Healthy")
trainAddMetaData_lvplH$Condition <- relevel(as.factor(as.character(trainAddMetaData_lvplH$Condition)), "Healthy")
trainAddMetaData_lvVbH$Condition <- relevel(as.factor(as.character(trainAddMetaData_lvVbH$Condition)), "Healthy")



# Get number of conditions per site
freqCondPerSite <- melt(data.frame(table(trainMetaData[trainMetaData$Age_Group == "Child",][,c("Condition","Site")])))

# Get proportion of conditions per site
freqCondPerSite["Prop"] <- NA
sites <- unique(freqCondPerSite$Site)
for (i in seq_len(length(sites))) { freqCondPerSite[freqCondPerSite$Site == sites[[i]],]["Prop"] <- freqCondPerSite[freqCondPerSite$Site == sites[[i]],]$value / sum(freqCondPerSite[freqCondPerSite$Site == sites[[i]],]$value) }

# Plot freq
plt1 <- ggplot(freqCondPerSite,
       aes(x=Site, y=value, fill=Condition)) + 
  geom_bar(stat="identity", color="black", position="stack") +
  scale_fill_manual(values=condCols) +
  ylab("Number of\nSamples") +
  theme_classic(base_size=36) + theme(axis.title=element_text(size=64),
                                 axis.text=element_text(size=42))

# Plot freq
plt2 <- ggplot(freqCondPerSite,
       aes(x=Site, y=Prop, fill=Condition)) + 
  geom_bar(stat="identity", color="black", position="stack") +
  scale_fill_manual(values=condCols) +
  ylab("Proportion of\nSamples") +
  theme_classic(base_size=36) + theme(axis.title=element_text(size=64),
                                 axis.text=element_text(size=42))

plt1
plt2

# save
ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom",
                  paste0("sampleGroupNumbers_bySite_afterFiltration.png")),
       plt1, width=5000, height=3500, units="px")
ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom",
                  paste0("sampleGroupProportions_bySite_afterFiltration.png")),
       plt2, width=5000, height=3500, units="px")



# Add to list
newList <- list(plt1, plt2)
names(newList) <- c("conditionGroupNumbers_trainData",
                    "conditionGroupProportions_trainData")
plotList <- append(plotList,
                   newList)



# Get number of conditions per site
freqCondPerSite <- melt(data.frame(table(trainAddMetaData[trainAddMetaData$Age_Group == "Child",][,c("Condition","Site")])))

# Get proportion of conditions per site
freqCondPerSite["Prop"] <- NA
sites <- unique(freqCondPerSite$Site)
for (i in seq_len(length(sites))) { freqCondPerSite[freqCondPerSite$Site == sites[[i]],]["Prop"] <- freqCondPerSite[freqCondPerSite$Site == sites[[i]],]$value / sum(freqCondPerSite[freqCondPerSite$Site == sites[[i]],]$value) }

# Plot freq
plt1 <- ggplot(freqCondPerSite,
       aes(x=Site, y=value, fill=Condition)) + 
  geom_bar(stat="identity", color="black", position="stack") +
  scale_fill_manual(values=condCols) +
  ylab("Number of\nSamples") +
  theme_classic(base_size=36) + theme(axis.title=element_text(size=64),
                                 axis.text=element_text(size=42))

# Plot freq
plt2 <- ggplot(freqCondPerSite,
       aes(x=Site, y=Prop, fill=Condition)) + 
  geom_bar(stat="identity", color="black", position="stack") +
  scale_fill_manual(values=condCols) +
  ylab("Proportion of\nSamples") +
  theme_classic(base_size=36) + theme(axis.title=element_text(size=64),
                                 axis.text=element_text(size=42))

plt1
plt2

# save
ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom",
                  paste0("sampleGroupNumbers_bySite_afterFiltration_justSamplesWithDate.png")),
       plt1, width=5000, height=3500, units="px")
ggsave(here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom",
                  paste0("sampleGroupProportions_bySite_afterFiltration_justSamplesWithDate.png")),
       plt2, width=5000, height=3500, units="px")

# Add to list
newList <- list(plt1, plt2)
names(newList) <- c("conditionGroupNumbers_trainData_justSamplesWithDate",
                    "conditionGroupProportions_trainData_justSamplesWithDate")
plotList <- append(plotList,
                   newList)



# Make designs
childCrossDesign         <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_condition), 
                                         data=trainMetaData)
childCrossDesignSite     <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_site), 
                                         data=trainMetaData)
childCrossDesign_noVbilt <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_condition), 
                                         data=trainMetaData_noVbilt)

# If not just CNO/Healthy
if (checkAllConditions) {

  childCrossDesign_incVblt <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_condition), 
                                           data=trainMetaData_incVblt)
  childCrossDesign_dresden <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_condition), 
                                           data=trainMetaData_dresden)
  is.fullrank(childCrossDesign_incVblt)
  is.fullrank(childCrossDesign_dresden)

}

childCrossDesign_noLivpl <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_condition), 
                                         data=trainMetaData_noLivpl)
childCrossDesign_lvplH   <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_condition), 
                                         data=trainMetaData_lvplH)
childCrossDesign_lvVbH   <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_condition), 
                                         data=trainMetaData_lvVbH)

childAddDesgin       <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate), 
                                                data=trainAddMetaData)

# If not just CNO/Healthy
if (checkAllConditions) {
  
  childAddDesgin_drsdn <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate), 
                                                  data=trainAddMetaData_drsdn)
  is.fullrank(childAddDesgin_drsdn)
  
}

childAddDesgin_drsVb <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate), 
                                                data=trainAddMetaData_drsVb)
childAddDesgin_drsLp <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate), 
                                                data=trainAddMetaData_drsLp)
childAddDesgin_lvplH <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate), 
                                                data=trainAddMetaData_lvplH)
childAddDesgin_lvVbH <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate), 
                                                data=trainAddMetaData_lvVbH)

# Check full rank
is.fullrank(childCrossDesign)
is.fullrank(childCrossDesignSite)
is.fullrank(childCrossDesign_noVbilt)
is.fullrank(childCrossDesign_noLivpl)
is.fullrank(childCrossDesign_lvplH)
is.fullrank(childCrossDesign_lvVbH)

is.fullrank(childAddDesgin)
is.fullrank(childAddDesgin_drsVb)
is.fullrank(childAddDesgin_drsLp)
is.fullrank(childAddDesgin_lvplH)



# Define levels to test
levels2Test <- paste0("Condition",unique(trainMetaData$Condition))
levels2Test_site <- paste0("Site",unique(trainMetaData$Site))

# Run Limma models
childLimmaRes_ANOVASite              <- runLimma(trainQuantTable[,rownames(childCrossDesignSite)], designM=childCrossDesignSite,
                                                 levels2Test=levels2Test_site[which(levels2Test_site %in% colnames(childCrossDesignSite))])
                                                 

childLimmaRes_ANOVACondition         <- runLimma(trainQuantTable[,rownames(childCrossDesign)], designM=childCrossDesign,
                                                 levels2Test=levels2Test[which(levels2Test %in% colnames(childCrossDesign))])
                                                 

childLimmaRes_ANOVACondition_noVbilt <- runLimma(trainQuantTable[,rownames(childCrossDesign_noVbilt)], designM=childCrossDesign_noVbilt, 
                                                 levels2Test=levels2Test[which(levels2Test %in% colnames(childCrossDesign_noVbilt))])

# If not just CNO/Healthy
if (checkAllConditions) {
  
  childLimmaRes_ANOVACondition_incVblt <- runLimma(trainQuantTable[,rownames(childCrossDesign_incVblt)], designM=childCrossDesign_incVblt, 
                                                   levels2Test=levels2Test[which(levels2Test %in% colnames(childCrossDesign_incVblt))])

  childLimmaRes_ANOVACondition_dresden <- runLimma(trainQuantTable[,rownames(childCrossDesign_dresden)], designM=childCrossDesign_dresden, 
                                                   levels2Test=levels2Test[which(levels2Test %in% colnames(childCrossDesign_dresden))])
  
}
                                                 

childLimmaRes_ANOVACondition_noLivpl <- runLimma(trainQuantTable[,rownames(childCrossDesign_noLivpl)], designM=childCrossDesign_noLivpl, 
                                                 levels2Test=levels2Test[which(levels2Test %in% colnames(childCrossDesign_noLivpl))])
                                                 

childLimmaRes_ANOVACondition_lvplH <- runLimma(trainQuantTable[,rownames(childCrossDesign_lvplH)], designM=childCrossDesign_lvplH, 
                                                 levels2Test=levels2Test[which(levels2Test %in% colnames(childCrossDesign_lvplH))])
                                                 

childLimmaRes_ANOVACondition_lvVbH <- runLimma(trainQuantTable[,rownames(childCrossDesign_lvVbH)], designM=childCrossDesign_lvVbH, 
                                                 levels2Test=levels2Test[which(levels2Test %in% colnames(childCrossDesign_lvVbH))])



childLimmaRes_ANOVACondition_addMeta       <- runLimma(trainQuantTable[,rownames(childAddDesgin)], designM=childAddDesgin, 
                                                       levels2Test=levels2Test[which(levels2Test %in% colnames(childAddDesgin))])
                                                 
# If not just CNO/Healthy
if (checkAllConditions) {
  
  childLimmaRes_ANOVACondition_addMeta_drsdn <- runLimma(trainQuantTable[,rownames(childAddDesgin_drsdn)], designM=childAddDesgin_drsdn, 
                                                         levels2Test=levels2Test[which(levels2Test %in% colnames(childAddDesgin_drsdn))])
  
}
                                                 

childLimmaRes_ANOVACondition_addMeta_drsVb <- runLimma(trainQuantTable[,rownames(childAddDesgin_drsVb)], designM=childAddDesgin_drsVb, 
                                                       levels2Test=levels2Test[which(levels2Test %in% colnames(childAddDesgin_drsVb))])
                                                 

childLimmaRes_ANOVACondition_addMeta_drsLp <- runLimma(trainQuantTable[,rownames(childAddDesgin_drsLp)], designM=childAddDesgin_drsLp, 
                                                       levels2Test=levels2Test[which(levels2Test %in% colnames(childAddDesgin_drsLp))])
                                                 

childLimmaRes_ANOVACondition_addMeta_lvplH <- runLimma(trainQuantTable[,rownames(childAddDesgin_lvplH)], designM=childAddDesgin_lvplH, 
                                                       levels2Test=levels2Test[which(levels2Test %in% colnames(childAddDesgin_lvplH))])
                                                 

childLimmaRes_ANOVACondition_addMeta_lvVbH <- runLimma(trainQuantTable[,rownames(childAddDesgin_lvVbH)], designM=childAddDesgin_lvVbH, 
                                                       levels2Test=levels2Test[which(levels2Test %in% colnames(childAddDesgin_lvVbH))])
                                                 


# Check number of sig
childLimmaRes_ANOVACondition[childLimmaRes_ANOVACondition$adj.P.Val < 0.05,]
childLimmaRes_ANOVASite[childLimmaRes_ANOVASite$adj.P.Val < 0.05,]
childLimmaRes_ANOVACondition_noVbilt[childLimmaRes_ANOVACondition_noVbilt$adj.P.Val < 0.05,]
# If not just CNO/Healthy
if (checkAllConditions) {
  childLimmaRes_ANOVACondition_incVblt[childLimmaRes_ANOVACondition_incVblt$adj.P.Val < 0.05,]
  childLimmaRes_ANOVACondition_dresden[childLimmaRes_ANOVACondition_dresden$adj.P.Val < 0.05,]
}

childLimmaRes_ANOVACondition_noLivpl[childLimmaRes_ANOVACondition_noLivpl$adj.P.Val < 0.05,]
childLimmaRes_ANOVACondition_lvplH[childLimmaRes_ANOVACondition_lvplH$adj.P.Val < 0.05,]
childLimmaRes_ANOVACondition_lvVbH[childLimmaRes_ANOVACondition_lvVbH$adj.P.Val < 0.05,]
childLimmaRes_ANOVACondition_addMeta[childLimmaRes_ANOVACondition_addMeta$adj.P.Val < 0.05,]
# If not just CNO/Healthy
if (checkAllConditions) {
  childLimmaRes_ANOVACondition_addMeta_drsdn[childLimmaRes_ANOVACondition_addMeta_drsdn$adj.P.Val < 0.05,]
} 
childLimmaRes_ANOVACondition_addMeta_drsVb[childLimmaRes_ANOVACondition_addMeta_drsVb$adj.P.Val < 0.05,]
childLimmaRes_ANOVACondition_addMeta_drsLp[childLimmaRes_ANOVACondition_addMeta_drsLp$adj.P.Val < 0.05,]
childLimmaRes_ANOVACondition_addMeta_lvplH[childLimmaRes_ANOVACondition_addMeta_lvplH$adj.P.Val < 0.05,]
childLimmaRes_ANOVACondition_addMeta_lvVbH[childLimmaRes_ANOVACondition_addMeta_lvVbH$adj.P.Val < 0.05,]

# Define proteins associated with site
siteProts <- childLimmaRes_ANOVASite[childLimmaRes_ANOVASite$adj.P.Val < 0.05,]$PG_ID

# Remove "Condition"
colnames(childLimmaRes_ANOVACondition) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition))
colnames(childLimmaRes_ANOVASite) <- gsub("Condition","",colnames(childLimmaRes_ANOVASite))
colnames(childLimmaRes_ANOVACondition_noVbilt) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_noVbilt))
# If not just CNO/Healthy
if (checkAllConditions) {
  colnames(childLimmaRes_ANOVACondition_incVblt) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_incVblt))
  colnames(childLimmaRes_ANOVACondition_dresden) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_dresden))
}
colnames(childLimmaRes_ANOVACondition_lvplH) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_lvplH))
colnames(childLimmaRes_ANOVACondition_lvVbH) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_lvVbH))

colnames(childLimmaRes_ANOVACondition_addMeta) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_addMeta))
# If not just CNO/Healthy
if (checkAllConditions) {
  colnames(childLimmaRes_ANOVACondition_addMeta_drsdn) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_addMeta_drsdn))
}
colnames(childLimmaRes_ANOVACondition_addMeta_drsVb) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_addMeta_drsVb))
colnames(childLimmaRes_ANOVACondition_addMeta_drsLp) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_addMeta_drsLp))
colnames(childLimmaRes_ANOVACondition_addMeta_lvplH) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_addMeta_lvplH))
colnames(childLimmaRes_ANOVACondition_addMeta_lvVbH) <- gsub("Condition","",colnames(childLimmaRes_ANOVACondition_addMeta_lvVbH))

# If not just CNO/Healthy
if (checkAllConditions) {
  
  # Intersect within-Dresden and within-Liverpool
  intersect(childLimmaRes_ANOVACondition_incVblt[childLimmaRes_ANOVACondition_incVblt$adj.P.Val < 0.05,]$PG_ID, 
            childLimmaRes_ANOVACondition_dresden[childLimmaRes_ANOVACondition_dresden$adj.P.Val < 0.05,]$PG_ID)
  
}

# Save
write.table(childLimmaRes_ANOVACondition,
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                       "limmaResults_ANOVAwrtCondition_crossSectionalData.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(childLimmaRes_ANOVASite,
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                       "limmaResults_ANOVAwrtSite_crossSectionalData.tsv"), sep="\t", quote=FALSE, row.names=FALSE)


# If not just CNO/Healthy
if (checkAllConditions) {
    
  write.table(childLimmaRes_ANOVACondition_dresden,
              here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                         "limmaResults_ANOVAwrtCondition_crossSectionalData_justDresdenSamples.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  
}



# Make Design_noInts
childCrossDesign_noInt         <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormula_condition_noIntercept), 
                                         data=trainMetaData)
childCrossDesign_noInt_noVbilt <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormula_condition_noIntercept), 
                                         data=trainMetaData_noVbilt)
# If not just CNO/Healthy
if (checkAllConditions) {
  childCrossDesign_noInt_dresden <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormula_condition_noIntercept), 
                                           data=trainMetaData_dresden)
}
childCrossDesign_noInt_noLivpl <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormula_condition_noIntercept), 
                                         data=trainMetaData_noLivpl)
childCrossDesign_noInt_lvplH <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormula_condition_noIntercept), 
                                         data=trainMetaData_lvplH)
childCrossDesign_noInt_lvVbH <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormula_condition_noIntercept), 
                                         data=trainMetaData_lvVbH)

childAddDesign_noInt       <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate_noIntercept), 
                                           data=trainAddMetaData)
# If not just CNO/Healthy
if (checkAllConditions) {
  childAddDesign_noInt_drsdn <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate_noIntercept), 
                                             data=trainAddMetaData_drsdn)
}
childAddDesign_noInt_drsVb <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate_noIntercept), 
                                           data=trainAddMetaData_drsVb)
childAddDesign_noInt_drsLp <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate_noIntercept), 
                                           data=trainAddMetaData_drsLp)
childAddDesign_noInt_lvplH <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate_noIntercept), 
                                           data=trainAddMetaData_lvplH)
childAddDesign_noInt_lvVbH <- model.matrix(as.formula(childCross_modelFormulas$childCrossFormulas_withDate_noIntercept), 
                                           data=trainAddMetaData_lvVbH)

# Fix column names
colnames(childCrossDesign_noInt) <- gsub("[:]","_",colnames(childCrossDesign_noInt))
colnames(childCrossDesign_noInt_noVbilt) <- gsub("[:]","_",colnames(childCrossDesign_noInt_noVbilt))
# If not just CNO/Healthy
if (checkAllConditions) {
  colnames(childCrossDesign_noInt_dresden) <- gsub("[:]","_",colnames(childCrossDesign_noInt_dresden))
}
colnames(childCrossDesign_noInt_noLivpl) <- gsub("[:]","_",colnames(childCrossDesign_noInt_noLivpl))
colnames(childCrossDesign_noInt_lvplH) <- gsub("[:]","_",colnames(childCrossDesign_noInt_lvplH))
colnames(childCrossDesign_noInt_lvVbH) <- gsub("[:]","_",colnames(childCrossDesign_noInt_lvVbH))

colnames(childAddDesign_noInt) <- gsub("[:]","_",colnames(childAddDesign_noInt))
# If not just CNO/Healthy
if (checkAllConditions) {
  colnames(childAddDesign_noInt_drsdn) <- gsub("[:]","_",colnames(childAddDesign_noInt_drsdn))
}
colnames(childAddDesign_noInt_drsVb) <- gsub("[:]","_",colnames(childAddDesign_noInt_drsVb))
colnames(childAddDesign_noInt_drsLp) <- gsub("[:]","_",colnames(childAddDesign_noInt_drsLp))
colnames(childAddDesign_noInt_lvplH) <- gsub("[:]","_",colnames(childAddDesign_noInt_lvplH))
colnames(childAddDesign_noInt_lvVbH) <- gsub("[:]","_",colnames(childAddDesign_noInt_lvVbH))

# Check full rank
is.fullrank(childCrossDesign_noInt)
is.fullrank(childCrossDesign_noInt_noVbilt)
# If not just CNO/Healthy
if (checkAllConditions) { is.fullrank(childCrossDesign_noInt_dresden) }
is.fullrank(childCrossDesign_noInt_noLivpl)
is.fullrank(childCrossDesign_noInt_lvplH)
is.fullrank(childCrossDesign_noInt_lvVbH)

is.fullrank(childAddDesign_noInt)
# If not just CNO/Healthy
if (checkAllConditions) { is.fullrank(childAddDesign_noInt_drsdn) }
is.fullrank(childAddDesign_noInt_drsVb)
is.fullrank(childAddDesign_noInt_drsLp)
is.fullrank(childAddDesign_noInt_lvplH)
is.fullrank(childAddDesign_noInt_lvplH)
is.fullrank(childAddDesign_noInt_lvVbH)



# If not just CNO/Healthy
if (checkAllConditions) {
  
  # Make constrasts
  childCross_contrastM <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                        CNOvsCrohn   = ConditionCNO-ConditionCrohn,
                                        CNOvsOnc      = ConditionCNO-ConditionOncology,
                                        CNOvsInfOsteo = ConditionCNO-ConditionInfOsteitis,
                                        levels=childCrossDesign_noInt)
  childCross_contrastM_noVbilt <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                                CNOvsCrohn   = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc      = ConditionCNO-ConditionOncology,
                                                levels=childCrossDesign_noInt_noVbilt)
  childCross_contrastM_dresden <- makeContrasts(CNOvsCrohn  = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc     = ConditionCNO-ConditionOncology,
                                                levels=childCrossDesign_noInt_dresden)
  childCross_contrastM_noLivpl <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                                CNOvsCrohn   = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc      = ConditionCNO-ConditionOncology,
                                                CNOvsInfOsteo = ConditionCNO-ConditionInfOsteitis,
                                                levels=childCrossDesign_noInt_noLivpl)
  childCross_contrastM_lvplH <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                                CNOvsCrohn   = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc      = ConditionCNO-ConditionOncology,
                                                CNOvsInfOsteo = ConditionCNO-ConditionInfOsteitis,
                                                levels=childCrossDesign_noInt_lvplH)
  childCross_contrastM_lvVbH <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                                CNOvsCrohn   = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc      = ConditionCNO-ConditionOncology,
                                                levels=childCrossDesign_noInt_lvVbH)
  childCross_contrastM_addMeta <- makeContrasts(CNOvsCrohn  = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc     = ConditionCNO-ConditionOncology,
                                                CNOvsInfOsteo = ConditionCNO-ConditionInfOsteitis,
                                                CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt)
  childCross_contrastM_addMeta_drsdn <- makeContrasts(CNOvsCrohn  = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc     = ConditionCNO-ConditionOncology,
                                                levels=childAddDesign_noInt_drsdn)
  childCross_contrastM_addMeta_drsVb <- makeContrasts(CNOvsCrohn  = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc     = ConditionCNO-ConditionOncology,
                                                CNOvsInfOsteo = ConditionCNO-ConditionInfOsteitis,
                                                CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt_drsVb)
  childCross_contrastM_addMeta_drsLp <- makeContrasts(CNOvsCrohn  = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc     = ConditionCNO-ConditionOncology,
                                                CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt_drsLp)
  childCross_contrastM_addMeta_lvplH <- makeContrasts(CNOvsCrohn  = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc     = ConditionCNO-ConditionOncology,
                                                CNOvsInfOsteo = ConditionCNO-ConditionInfOsteitis,
                                                CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt_lvplH)
  childCross_contrastM_addMeta_lvVbH <- makeContrasts(CNOvsCrohn  = ConditionCNO-ConditionCrohn,
                                                CNOvsOnc     = ConditionCNO-ConditionOncology,
                                                CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt_lvVbH)
  
  childCross_contrastM_CNOvsAll <- makeContrasts(CNOvsInfOsteo = ConditionCNO-ConditionInfOsteitis,
                                                 CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                 CNOvsCrohn = ConditionCNO-ConditionCrohn,
                                                 CNOvsOncology = ConditionCNO-ConditionOncology,
                                                 levels=childCrossDesign_noInt)
  childCross_contrastM_CNOvsIOH <- makeContrasts(CNOvsInfOsteo = ConditionCNO-ConditionInfOsteitis,
                                                 CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                 levels=childCrossDesign_noInt)

} else {
  
  # Make constrasts
  childCross_contrastM <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                        levels=childCrossDesign_noInt)
  childCross_contrastM_noVbilt <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                                levels=childCrossDesign_noInt_noVbilt)
  childCross_contrastM_noLivpl <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                                levels=childCrossDesign_noInt_noLivpl)
  childCross_contrastM_lvplH <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                                levels=childCrossDesign_noInt_lvplH)
  childCross_contrastM_lvVbH <- makeContrasts(CNOvsHealthy  = ConditionCNO-ConditionHealthy,
                                                levels=childCrossDesign_noInt_lvVbH)
  childCross_contrastM_addMeta <- makeContrasts(CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt)
  childCross_contrastM_addMeta_drsVb <- makeContrasts(CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt_drsVb)
  childCross_contrastM_addMeta_drsLp <- makeContrasts(CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt_drsLp)
  childCross_contrastM_addMeta_lvplH <- makeContrasts(CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt_lvplH)
  childCross_contrastM_addMeta_lvVbH <- makeContrasts(CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                levels=childAddDesign_noInt_lvVbH)
  
  childCross_contrastM_CNOvsAll <- makeContrasts(CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                 levels=childCrossDesign_noInt)
  childCross_contrastM_CNOvsIOH <- makeContrasts(CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                 levels=childCrossDesign_noInt)
  
}  



# Run Limma models
childLimmaRes_vsCNO_ANOVA_Condition         <- runLimma(trainQuantTable[,rownames(childCrossDesign_noInt)], 
                                                        designM=childCrossDesign_noInt,
                                                        levels2Test=colnames(childCross_contrastM),
                                                        contrasts2Fit=childCross_contrastM)
                                                 
childLimmaRes_vsCNO_ANOVA_Condition_noVbilt <- runLimma(trainQuantTable[,rownames(childCrossDesign_noInt_noVbilt)], 
                                                        designM=childCrossDesign_noInt_noVbilt, 
                                                        levels2Test=colnames(childCross_contrastM_noVbilt),
                                                        contrasts2Fit=childCross_contrastM_noVbilt)
# If not just CNO/Healthy
if (checkAllConditions) {
  childLimmaRes_vsCNO_ANOVA_Condition_dresden <- runLimma(trainQuantTable[,rownames(childCrossDesign_noInt_dresden)], 
                                                          designM=childCrossDesign_noInt_dresden, 
                                                          levels2Test=colnames(childCross_contrastM_dresden),
                                                          contrasts2Fit=childCross_contrastM_dresden)
}
childLimmaRes_vsCNO_ANOVA_Condition_noLivpl <- runLimma(trainQuantTable[,rownames(childCrossDesign_noInt_noLivpl)], 
                                                        designM=childCrossDesign_noInt_noLivpl,  
                                                        levels2Test=colnames(childCross_contrastM_noLivpl),
                                                        contrasts2Fit=childCross_contrastM_noLivpl)
childLimmaRes_vsCNO_ANOVA_Condition_lvplH <- runLimma(trainQuantTable[,rownames(childCrossDesign_noInt_lvplH)], 
                                                        designM=childCrossDesign_noInt_lvplH,  
                                                        levels2Test=colnames(childCross_contrastM_lvplH),
                                                        contrasts2Fit=childCross_contrastM_lvplH)
                                                 
childLimmaRes_vsCNO_ANOVA_Condition_lvVbH <- runLimma(trainQuantTable[,rownames(childCrossDesign_noInt_lvVbH)], 
                                                        designM=childCrossDesign_noInt_lvVbH,  
                                                        levels2Test=colnames(childCross_contrastM_lvVbH),
                                                        contrasts2Fit=childCross_contrastM_lvVbH)
                                                 


childLimmaRes_vsCNO_ANOVA_Condition_addMeta <- runLimma(normQuantTable[,rownames(childAddDesign_noInt)], 
                                                        designM=childAddDesign_noInt, 
                                                        levels2Test=colnames(childCross_contrastM_addMeta),
                                                        contrasts2Fit=childCross_contrastM_addMeta)
# If not just CNO/Healthy
if (checkAllConditions) {
  childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn <- runLimma(normQuantTable[,rownames(childAddDesign_noInt_drsdn)], 
                                                          designM=childAddDesign_noInt_drsdn, 
                                                          levels2Test=colnames(childCross_contrastM_addMeta_drsdn),
                                                          contrasts2Fit=childCross_contrastM_addMeta_drsdn)
}
childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsVb <- runLimma(normQuantTable[,rownames(childAddDesign_noInt_drsVb)], 
                                                        designM=childAddDesign_noInt_drsVb, 
                                                        levels2Test=colnames(childCross_contrastM_addMeta_drsVb),
                                                        contrasts2Fit=childCross_contrastM_addMeta_drsVb)
                                                 
childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsLp <- runLimma(normQuantTable[,rownames(childAddDesign_noInt_drsLp)], 
                                                        designM=childAddDesign_noInt_drsLp, 
                                                        levels2Test=colnames(childCross_contrastM_addMeta_drsLp),
                                                        contrasts2Fit=childCross_contrastM_addMeta_drsLp)
                                                 
childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvplH <- runLimma(normQuantTable[,rownames(childAddDesign_noInt_lvplH)], 
                                                        designM=childAddDesign_noInt_lvplH, 
                                                        levels2Test=colnames(childCross_contrastM_addMeta_lvplH),
                                                        contrasts2Fit=childCross_contrastM_addMeta_lvplH)
                                                 
childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvVbH <- runLimma(normQuantTable[,rownames(childAddDesign_noInt_lvVbH)], 
                                                        designM=childAddDesign_noInt_lvVbH, 
                                                        levels2Test=colnames(childCross_contrastM_addMeta_lvVbH),
                                                        contrasts2Fit=childCross_contrastM_addMeta_lvVbH)
                                                 

# Check number of sig
childLimmaRes_vsCNO_ANOVA_Condition[childLimmaRes_vsCNO_ANOVA_Condition$adj.P.Val < 0.05,]
childLimmaRes_vsCNO_ANOVA_Condition_noVbilt[childLimmaRes_vsCNO_ANOVA_Condition_noVbilt$adj.P.Val < 0.05,]
if (checkAllConditions) { childLimmaRes_vsCNO_ANOVA_Condition_dresden[childLimmaRes_vsCNO_ANOVA_Condition_dresden$adj.P.Val < 0.05,] }
childLimmaRes_vsCNO_ANOVA_Condition_noLivpl[childLimmaRes_vsCNO_ANOVA_Condition_noLivpl$adj.P.Val < 0.05,]
childLimmaRes_vsCNO_ANOVA_Condition_lvplH[childLimmaRes_vsCNO_ANOVA_Condition_lvplH$adj.P.Val < 0.05,]
childLimmaRes_vsCNO_ANOVA_Condition_lvVbH[childLimmaRes_vsCNO_ANOVA_Condition_lvVbH$adj.P.Val < 0.05,]

childLimmaRes_vsCNO_ANOVA_Condition_addMeta[childLimmaRes_vsCNO_ANOVA_Condition_addMeta$adj.P.Val < 0.05,]
if (checkAllConditions) { childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn[childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn$adj.P.Val < 0.05,] }
childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsVb[childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsVb$adj.P.Val < 0.05,]
childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsLp[childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsLp$adj.P.Val < 0.05,]
childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvplH[childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvplH$adj.P.Val < 0.05,]
childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvVbH[childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvVbH$adj.P.Val < 0.05,]

# Remove "Condition"
colnames(childLimmaRes_vsCNO_ANOVA_Condition) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition))
colnames(childLimmaRes_vsCNO_ANOVA_Condition_noVbilt) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_noVbilt))
if (checkAllConditions) { colnames(childLimmaRes_vsCNO_ANOVA_Condition_dresden) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_dresden)) }
colnames(childLimmaRes_vsCNO_ANOVA_Condition_noLivpl) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_noLivpl))
colnames(childLimmaRes_vsCNO_ANOVA_Condition_lvplH) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_lvplH))
colnames(childLimmaRes_vsCNO_ANOVA_Condition_lvVbH) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_lvVbH))

colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta))
if (checkAllConditions) { colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn)) }
colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsVb) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsVb))
colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsLp) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsLp))
colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvplH) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvplH))
colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvVbH) <- gsub("Condition","",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_lvVbH))
  
# Define all conditions/sites proteins
allProts_padj <- childLimmaRes_vsCNO_ANOVA_Condition[childLimmaRes_vsCNO_ANOVA_Condition$adj.P.Val < candidatePThr &
                                           apply(childLimmaRes_vsCNO_ANOVA_Condition[,grepl("^CNO|logFC",colnames(childLimmaRes_vsCNO_ANOVA_Condition)),drop=FALSE], 1, function(row) any(abs(row) > log2ThrVal)),]$PG_ID
allProts_pval <- childLimmaRes_vsCNO_ANOVA_Condition[childLimmaRes_vsCNO_ANOVA_Condition$P.Value < candidatePThr &
                                           apply(childLimmaRes_vsCNO_ANOVA_Condition[,grepl("^CNO|logFC",colnames(childLimmaRes_vsCNO_ANOVA_Condition)),drop=FALSE], 1, function(row) any(abs(row) > log2ThrVal)),]$PG_ID
allSiteProts_padj <- childLimmaRes_ANOVASite[childLimmaRes_ANOVASite$adj.P.Val < candidatePThr &
                                             apply(childLimmaRes_ANOVASite[,grepl("^Site|logFC",colnames(childLimmaRes_ANOVASite)),drop=FALSE], 1, function(row) any(abs(row) > log2ThrVal)),]$PG_ID
allSiteProts_pval <- childLimmaRes_ANOVASite[childLimmaRes_ANOVASite$P.Value < candidatePThr &
                                             apply(childLimmaRes_ANOVASite[,grepl("^Site|logFC",colnames(childLimmaRes_ANOVASite)),drop=FALSE], 1, function(row) any(abs(row) > log2ThrVal)),]$PG_ID


# If we want to check all possible conditions (FALSE if focussing on just subsets of conditions)
if (checkAllConditions) { 
  
  # Define dresden proteins
  dresdenProts_padj <- childLimmaRes_vsCNO_ANOVA_Condition_dresden[childLimmaRes_vsCNO_ANOVA_Condition_dresden$adj.P.Val < candidatePThr &
                                                apply(childLimmaRes_vsCNO_ANOVA_Condition_dresden[,grepl("^CNO",colnames(childLimmaRes_vsCNO_ANOVA_Condition_dresden))], 1, function(row) any(abs(row) > log2ThrVal)),]$PG_ID
  
  dresdenDateProts_padj <- childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn[childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn$adj.P.Val < candidatePThr &
                                                apply(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn[,grepl("^CNO",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn))], 1, function(row) any(abs(row) > log2ThrVal)),]$PG_ID
  
  dresdenProts_pval <- childLimmaRes_vsCNO_ANOVA_Condition_dresden[childLimmaRes_vsCNO_ANOVA_Condition_dresden$P.Value < candidatePThr &
                                                apply(childLimmaRes_vsCNO_ANOVA_Condition_dresden[,grepl("^CNO",colnames(childLimmaRes_vsCNO_ANOVA_Condition_dresden))], 1, function(row) any(abs(row) > log2ThrVal)),]$PG_ID
  
  dresdenDateProts_pval <- childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn[childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn$P.Value < candidatePThr &
                                                apply(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn[,grepl("^CNO",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn))], 1, function(row) abs(any(abs(row) > log2ThrVal))),]$PG_ID

}

length(Reduce(intersect, list(allSiteProts_padj, allProts_padj, dresdenProts_padj, dresdenDateProts_padj)))

# Save
write.table(childLimmaRes_vsCNO_ANOVA_Condition,
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                       "limmaResults_ANOVA_wrtCNOvContrasts_crossSectionalData.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
# 
if (checkAllConditions) { 

write.table(childLimmaRes_vsCNO_ANOVA_Condition_dresden,
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                       "limmaResults_ANOVA_wrtCNOvContrasts_crossSectionalData_justDresdenSamples.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

}



# Loop over possible contratss
allPW_limmaRes <- data.frame()
for (contrastName in colnames(childCross_contrastM_CNOvsAll)) {

  # Test
  limmaRes <- runLimma(trainQuantTable[,rownames(childCrossDesign_noInt)], 
                       designM=childCrossDesign_noInt,  
                       levels2Test=contrastName,
                       contrasts2Fit=childCross_contrastM_CNOvsAll)
  limmaRes["CONTRAST"] <- contrastName
  allPW_limmaRes <- rbind(allPW_limmaRes, limmaRes)
  
}
iohPW_limmaRes <- data.frame()
for (contrastName in colnames(childCross_contrastM_CNOvsIOH)) {

  # Test
  limmaRes <- runLimma(trainQuantTable[,rownames(childCrossDesign_noInt)], 
                       designM=childCrossDesign_noInt,  
                       levels2Test=contrastName,
                       contrasts2Fit=childCross_contrastM_CNOvsIOH)
  limmaRes["CONTRAST"] <- contrastName
  iohPW_limmaRes <- rbind(iohPW_limmaRes, limmaRes)
  
}

# Enumerate hits
## Padj
nConsistentProts_allPW_padj <- c(table(allPW_limmaRes[allPW_limmaRes$adj.P.Val < candidatePThr & abs(allPW_limmaRes$logFC) > log2ThrVal,]$PG_ID))
nConsistentProts_iohPW_padj <- c(table(iohPW_limmaRes[iohPW_limmaRes$adj.P.Val < candidatePThr & abs(iohPW_limmaRes$logFC) > log2ThrVal,]$PG_ID))
consitentProts_allPW_padj <- names(nConsistentProts_allPW_padj[nConsistentProts_allPW_padj == ncol(childCross_contrastM_CNOvsAll)])
consitentProts_iohPW_padj <- names(nConsistentProts_iohPW_padj[nConsistentProts_iohPW_padj == ncol(childCross_contrastM_CNOvsIOH)])

## PVal
nConsistentProts_allPW_pval <- c(table(allPW_limmaRes[allPW_limmaRes$P.Value < candidatePThr & abs(allPW_limmaRes$logFC) > log2ThrVal,]$PG_ID))
nConsistentProts_iohPW_pval <- c(table(iohPW_limmaRes[iohPW_limmaRes$P.Value < candidatePThr & abs(iohPW_limmaRes$logFC) > log2ThrVal,]$PG_ID))
consitentProts_allPW_pval <- names(nConsistentProts_allPW_pval[nConsistentProts_allPW_pval == ncol(childCross_contrastM_CNOvsAll)])
consitentProts_iohPW_pval <- names(nConsistentProts_iohPW_pval[nConsistentProts_iohPW_pval == ncol(childCross_contrastM_CNOvsIOH)])

# Save
write.table(allPW_limmaRes,
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                       "limmaResults_tTestCNOvsConditions_pairwiseContrasts_crossSectionalData.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(iohPW_limmaRes,
            here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                       "limmaResults_tTestCNOvsHealthyInfOsteo_pairwiseContrasts_crossSectionalData.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

# Overlap dresden sets with pairiwse ioh sets
if (checkAllConditions) { 
    
  length(intersect(dresdenProts_padj, consitentProts_iohPW_padj))
  length(intersect(dresdenDateProts_padj, consitentProts_iohPW_padj))

}



# Get filtered
allSiteSig <- childLimmaRes_vsCNO_ANOVA_Condition[childLimmaRes_vsCNO_ANOVA_Condition$adj.P.Val < 0.05 &
                                                  apply(childLimmaRes_vsCNO_ANOVA_Condition[,grepl("^CNO",colnames(childLimmaRes_vsCNO_ANOVA_Condition))], 1, function(row) any(abs(row) > log2ThrVal)),]
dresdenSig <- childLimmaRes_vsCNO_ANOVA_Condition_dresden[childLimmaRes_vsCNO_ANOVA_Condition_dresden$adj.P.Val < 0.05 &
                                                  apply(childLimmaRes_vsCNO_ANOVA_Condition_dresden[,grepl("^CNO",colnames(childLimmaRes_vsCNO_ANOVA_Condition_dresden))], 1, function(row) any(abs(row) > log2ThrVal)),]
allSiteDateSig <- childLimmaRes_vsCNO_ANOVA_Condition_addMeta[childLimmaRes_vsCNO_ANOVA_Condition_addMeta$adj.P.Val < 0.05 &
                                                  apply(childLimmaRes_vsCNO_ANOVA_Condition_addMeta[,grepl("^CNO",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta))], 1, function(row) any(abs(row) > log2ThrVal)),]
dresdenDateSig <- childLimmaRes_vsCNO_ANOVA_Condition[childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn$adj.P.Val < 0.05 &
                                                  apply(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn[,grepl("^CNO",colnames(childLimmaRes_vsCNO_ANOVA_Condition_addMeta_drsdn))], 1, function(row) any(abs(row) > log2ThrVal)),]



# Fit probabilistic PCA
set.seed(35353)
resPPCA <- pcaMethods::pca(t(trainQuantTable[ metaCrossSectional_child$Sample_ID]), 
                           method="ppca", center=TRUE, scale="pareto", cv="q2", nPcs=10)

# Combine with metadata
resPPCA_data <- merge(metaCrossSectional_child, pcaMethods::scores(resPPCA), by.x="Sample_ID", "row.names")

# Get variance explained
varExpl <- list(variance.percent=(sDev(resPPCA) / sum(sDev(resPPCA)))  * 100)



# Plot
site_plt <- Plot2DPCA(meta=resPPCA_data, var=varExpl, shape_by="Age_Group", color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=7) + 
    scale_shape_manual(values=c(21, 24)) + 
    scale_fill_manual(values=siteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Site") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
cond_plt <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=7) + 
    scale_fill_manual(values=tidyCondCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
plate_plt <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=7) + 
    scale_fill_manual(values=plateCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72))
condSite_plt1 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=7) + 
      scale_fill_manual(values=condSiteCols) +
      theme_classic(base_size=36) +
      theme(legend.position="top", 
            legend.direction="horizontal",
            legend.title.position="left",
            axis.text=element_text(size=48),
            legend.title=element_text(size=88),
            legend.text=element_text(size=64),
            axis.title=element_text(size=104),
            plot.title=element_text(size=116))
condSite_plt2 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=7) + 
      scale_fill_manual(values=condSiteCols) +
      theme_classic(base_size=36) +
      theme(legend.position="top", 
            legend.direction="horizontal",
            legend.title.position="left",
            axis.text=element_text(size=48),
            legend.title=element_text(size=88),
            legend.text=element_text(size=64),
            axis.title=element_text(size=104),
            plot.title=element_text(size=116))
condSite_plt3 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=7) + 
      scale_fill_manual(values=condSiteCols) +
      theme_classic(base_size=36) +
      theme(legend.position="top", 
            legend.direction="horizontal",
            legend.title.position="left",
            axis.text=element_text(size=48),
            legend.title=element_text(size=88),
            legend.text=element_text(size=64),
            axis.title=element_text(size=104),
            plot.title=element_text(size=116))
condSite_plt4 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=7) + 
      scale_fill_manual(values=condSiteCols) +
      theme_classic(base_size=36) +
      theme(legend.position="top", 
            legend.direction="horizontal",
            legend.title.position="left",
            axis.text=element_text(size=48),
            legend.title=element_text(size=88),
            legend.text=element_text(size=64),
            axis.title=element_text(size=104),
            plot.title=element_text(size=116))
    labs(fill="CondSite") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt5 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=7) + 
      scale_fill_manual(values=condSiteCols) +
      theme_classic(base_size=36) +
      theme(legend.position="top", 
            legend.direction="horizontal",
            legend.title.position="left",
            axis.text=element_text(size=48),
            legend.title=element_text(size=88),
            legend.text=element_text(size=64),
            axis.title=element_text(size=104),
            plot.title=element_text(size=116))
condSite_plt6 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=7) + 
      scale_fill_manual(values=condSiteCols) +
      theme_classic(base_size=36) +
      theme(legend.position="top", 
            legend.direction="horizontal",
            legend.title.position="left",
            axis.text=element_text(size=48),
            legend.title=element_text(size=88),
            legend.text=element_text(size=64),
            axis.title=element_text(size=104),
            plot.title=element_text(size=116))

# Add to list
newList <- list(site_plt,
                cond_plt,
                plate_plt,
                condSite_plt1,
                condSite_plt2,
                condSite_plt3,
                condSite_plt4,
                condSite_plt5,
                condSite_plt6)
names(newList) <- c("allProtsPPCA_colBySite",
                    "allProtsPPCA_colByCondition",
                    "allProtsPPCA_colByPlate",
                    "allProtsPPCA_colByCondSitePC12",
                    "allProtsPPCA_colByCondSitePC23",
                    "allProtsPPCA_colByCondSitePC34",
                    "allProtsPPCA_colByCondSitePC45",
                    "allProtsPPCA_colByCondSitePC67",
                    "allProtsPPCA_colByCondSitePC78")
plotList <- append(plotList,
                   newList)



# Fit probabilistic PCA
set.seed(35353)
resPPCA <- pcaMethods::pca(t(trainQuantTable[dresdenSig$PG_ID,
                                             metaCrossSectional_child$Sample_ID]), 
                           method="ppca", center=TRUE, scale="pareto", cv="q2", nPcs=10)

# Combine with metadata
resPPCA_data <- merge(metaCrossSectional_child, pcaMethods::scores(resPPCA), by.x="Sample_ID", "row.names")

# Get variance explained
varExpl <- list(variance.percent=(sDev(resPPCA) / sum(sDev(resPPCA)))  * 100)



# Plot
site_plt <- Plot2DPCA(meta=resPPCA_data, var=varExpl, shape_by="Age_Group", color_by="Site", label_id="Sample_ID", pcs=c(1, 2), pch_size=9) + 
    scale_shape_manual(values=c(21, 24)) + 
    scale_fill_manual(values=siteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Site") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
cond_plt <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="Condition", label_id="Sample_ID", pcs=c(1,2), pch_size=9) + 
    scale_fill_manual(values=tidyCondCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
plate_plt <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="PLATE", label_id="Sample_ID", pcs=c(1, 2), pch_size=9) + 
    scale_fill_manual(values=plateCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt1 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt2 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(2, 3), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt3 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(3, 4), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt4 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(4, 5), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt5 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(5, 6), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt6 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(6, 7), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))

# Add to list
newList <- list(site_plt,
                cond_plt,
                plate_plt,
                condSite_plt1,
                condSite_plt2,
                condSite_plt3,
                condSite_plt4,
                condSite_plt5,
                condSite_plt6)
names(newList) <- c("sigDresdenProtsPPCA_colBySite",
                    "sigDresdenProtsPPCA_colByCondition",
                    "sigDresdenProtsPPCA_colByPlate",
                    "sigDresdenProtsPPCA_colByCondSitePC12",
                    "sigDresdenProtsPPCA_colByCondSitePC23",
                    "sigDresdenProtsPPCA_colByCondSitePC34",
                    "sigDresdenProtsPPCA_colByCondSitePC45",
                    "sigDresdenProtsPPCA_colByCondSitePC67",
                    "sigDresdenProtsPPCA_colByCondSitePC78")
plotList <- append(plotList,
                   newList)



# Fit probabilistic PCA
set.seed(35353)
resPPCA <- pcaMethods::pca(t(trainQuantTable[allSiteSig$PG_ID,
                                             metaCrossSectional_child$Sample_ID]), 
                           method="ppca", center=TRUE, scale="pareto", cv="q2", nPcs=10)

# Combine with metadata
resPPCA_data <- merge(metaCrossSectional_child, pcaMethods::scores(resPPCA), by.x="Sample_ID", "row.names")

# Get variance explained
varExpl <- list(variance.percent=(sDev(resPPCA) / sum(sDev(resPPCA)))  * 100)



# Plot
site_plt <- Plot2DPCA(meta=resPPCA_data, var=varExpl, shape_by="Age_Group", color_by="Site", label_id="Sample_ID", pcs=c(1, 2), pch_size=9) + 
    scale_shape_manual(values=c(21, 24)) + 
    scale_fill_manual(values=siteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Site") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
cond_plt <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="Condition", label_id="Sample_ID", pcs=c(1,2), pch_size=9) + 
    scale_fill_manual(values=tidyCondCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
plate_plt <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="PLATE", label_id="Sample_ID", pcs=c(1, 2), pch_size=9) + 
    scale_fill_manual(values=plateCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt1 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(1, 2), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt2 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(2, 3), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt3 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(3, 4), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt4 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(4, 5), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt5 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(5, 6), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
condSite_plt6 <- Plot2DPCA(meta=resPPCA_data, var=varExpl, color_by="CondSite", label_id="Sample_ID", pcs=c(6, 7), pch_size=9) + 
    scale_fill_manual(values=condSiteCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))

# Add to list
newList <- list(site_plt,
                cond_plt,
                plate_plt,
                condSite_plt1,
                condSite_plt2,
                condSite_plt3,
                condSite_plt4,
                condSite_plt5,
                condSite_plt6)
names(newList) <- c("sigProtsPPCA_colBySite",
                    "sigProtsPPCA_colByCondition",
                    "sigProtsPPCA_colByPlate",
                    "sigProtsPPCA_colByCondSitePC12",
                    "sigProtsPPCA_colByCondSitePC23",
                    "sigProtsPPCA_colByCondSitePC34",
                    "sigProtsPPCA_colByCondSitePC45",
                    "sigProtsPPCA_colByCondSitePC67",
                    "sigProtsPPCA_colByCondSitePC78")
plotList <- append(plotList,
                   newList)



# Plot venn diagram
venn.diagram(
  x = list(allSiteSig$PG_ID, 
           dresdenSig$PG_ID), #, 
          #  dresdenDateSig$PG_ID),
  category.names = c("All Sites",
                     "Dresden"),
                     # "Dresden\nwith Date"),
  filename = here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                        "sigProteinSet_overlapVennDiagram.png"),
  output=TRUE,
  
        # Output features
        imagetype="png" ,
        height = 10000 , 
        width = 10000 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 20,
        fill = siteCols[c("Liverpool", "Dresden")],
        col = c("white", "white"),
        
        # Numbers
        cex = 16,
        fontface = "bold",
        fontfamily = "sans",
  
        # Category labels
        cat.cex = 12,
        cat.dist = c(0.1, 0.1),
        cat.default.pos = "outer",
        cat.fontfamily = "sans",
  
        ext.line.lwd=0,
        ext.dist=-0.275,
        ext.length=0,
      
        # Margins
        margin = 0.1)# ,




# Plot venn diagram
venn.diagram(
  x = list(allSiteSig$PG_ID, 
           dresdenSig$PG_ID, 
           dresdenDateSig$PG_ID), #, 
          #  dresdenDateSig$PG_ID),
  category.names = c("All Sites",
                     "Dresden",
                     "Dresden\nDate-Adjusted"),
  filename = here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                        "sigProteinSet_overlapVennDiagram_withDateAdjustment.png"),
  output=TRUE,
  
        # Output features
        imagetype="png" ,
        height = 10000 , 
        width = 10000 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 20,
        fill = siteCols[c("Liverpool", "Dresden", "Naproxen")],
        col = c("white", "white", "white"),
        
        # Numbers
        cex = 14,
        fontface = "bold",
        fontfamily = "sans",
  
        # Category labels
        cat.cex = 10,
        cat.dist = c(0.125, 0.125, 0.3),
        cat.default.pos = "outer",
        cat.fontfamily = "sans",
      
        # Margins
        margin = 0.08)# ,



# Plot venn diagram
venn.diagram(
  x = list(allSiteProts_padj,
           allSiteSig$PG_ID, 
           dresdenSig$PG_ID, 
           allSiteDateSig$PG_ID),
           # dresdenDateSig$PG_ID),
  category.names = c("Site",
                     "All Sites",
                     "Dresden",
                     "All Sites\nwith Date"),
                     # "Dresden\nwith Date"),
  filename = here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance",
                        "sigProteinSet_overlapVennDiagram_incSite.png"),
  output=TRUE,
  
        # Output features
        imagetype="png" ,
        height = 10000 , 
        width = 10000 , 
        resolution = 300,
        compression = "lzw",
        
        # Circles
        lwd = 20,,
        fill = siteCols[c("Liverpool", "Leiden", "Dresden", "Naproxen")],
        col = c("white", "white", "white", "white"),
        
        # Numbers
        cex = 16,
        fontface = "bold",
        fontfamily = "sans",
  
        # Category labels
        cat.cex = 12,
        cat.dist = c(0.225, 0.225, 0.115, 0.15),
        cat.default.pos = "outer",
        cat.fontfamily = "sans",
      
        # Margins
        margin = 0.05)


# Design matrix
designMatrix <- model.matrix( ~ 0 + Condition + Z_Age + PLATE + Sex,
                             data=trainMetaData)
is.fullrank(designMatrix)

# Calculate duplicate correlation 
# Calculating this across all proteins for better assessment of patient correlation?
dupCorr <- duplicateCorrelation(trainQuantTable[,trainMetaData$Sample_ID],
                                design=designMatrix, 
                                block=trainMetaData$Site)

# Fit
fit <- lmFit(trainQuantTable[,trainMetaData$Sample_ID], 
             design=designMatrix[trainMetaData$Sample_ID,], 
             correlation=dupCorr$consensus.correlation, 
             block=trainMetaData$Site)

# Fit contrasts
fit <- contrasts.fit(fit,
                     contrast=childCross_contrastM)

# Ebayes
fit <- eBayes(fit, 
              trend=TRUE, 
              robust=TRUE)

# Get top hits
siteCorrRes <- topTable(fit, coef=colnames(childCross_contrastM), n=Inf)

commonProts <- intersect(rownames(siteCorrRes[siteCorrRes$adj.P.Val < 0.05,]),
          rownames(childLimmaRes_vsCNO_ANOVA_Condition[childLimmaRes_vsCNO_ANOVA_Condition$adj.P.Val < 0.05,]))



# If not just CNO/Healthy
write.table(trainMetaData,
          here::here(metaDir, paste0("trainData_justChildSamples_splitProp",trainSplitProp,"_allSites.tsv")), sep="\t", quote=FALSE, row.names=FALSE)



# Get sig CNO vs healthy
cnoVsHealthy_sigData <- iohPW_limmaRes[iohPW_limmaRes$adj.P.Val < candidatePThr & abs(iohPW_limmaRes$logFC) > log2ThrVal & iohPW_limmaRes$CONTRAST == "CNOvsHealthy",]

# Create contrast
childCross_contrastM_CNOvsIOH_withDate <- makeContrasts(CNOvsHealthy = ConditionCNO-ConditionHealthy,
                                                        levels=childAddDesign_noInt)
  
# Re-run with date
limmaRes_cnoHealthy_wDate <- runLimma(trainQuantTable[,rownames(childAddDesign_noInt)], 
                                      designM=childAddDesign_noInt,  
                                      levels2Test="CNOvsHealthy",
                                      contrasts2Fit=childCross_contrastM_CNOvsIOH_withDate)

# Get sig]
cnoHealthyDateProts <- limmaRes_cnoHealthy_wDate[limmaRes_cnoHealthy_wDate$adj.P.Val < 0.05 & abs(limmaRes_cnoHealthy_wDate$logFC) > log2ThrVal,]$PG_ID




# If not just CNO/Healthy
if (checkAllConditions) {
  
  print(paste0("chosenProteins_featureSelection_pThr",as.character(candidatePThr),"_treatLog2Thr",as.character(log2ThrVal),".rds"))
  
  # Generate list of outputs
  allProtSets <- list(limma=list(
                                 padj=list(AllSig=allProts_padj,
                                           AllSite=allSiteProts_padj,
                                           AllSite_withDate=allSiteDateSig$PG_ID,
                                           DresdenOnly=dresdenProts_padj,
                                           DresdenOnly_withDate=dresdenDateProts_padj,
                                           Dresden_allPW=intersect(dresdenProts_padj, consitentProts_allPW_padj),
                                           Dresden_infOsteoHealthyPW=intersect(dresdenProts_padj, consitentProts_iohPW_padj),
                                           Dresden_withDate_allPW=intersect(dresdenDateProts_padj, consitentProts_allPW_padj),
                                           Dresden_withDate_infOsteoHealthyPW=intersect(dresdenDateProts_padj, consitentProts_iohPW_padj)),
                                 pval=list(AllSig=allProts_pval,
                                           AllSite=allSiteProts_pval,
                                           DresdenOnly=dresdenProts_pval,
                                           DresdenOnly_withDate=dresdenDateProts_pval,
                                           Dresden_allPW=intersect(dresdenProts_pval, consitentProts_allPW_pval),
                                           Dresden_infOsteoHealthyPW=intersect(dresdenProts_pval, consitentProts_iohPW_pval),
                                           Dresden_withDate_allPW=intersect(dresdenDateProts_pval, consitentProts_allPW_pval),
                                           Dresden_withDate_infOsteoHealthyPW=intersect(dresdenDateProts_pval, consitentProts_iohPW_pval)),
                                 All=rownames(trainQuantTable))
          )
  
  # Save
  print(paste0("// Saving to: ",here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance", 
                     paste0("chosenProteins_featureSelection_pThr",as.character(candidatePThr),"_treatLog2Thr",as.character(log2ThrVal),".rds"))))
  saveRDS(allProtSets,
          here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance", 
                     paste0("chosenProteins_featureSelection_pThr",as.character(candidatePThr),"_treatLog2Thr",as.character(log2ThrVal),".rds")))
 
} else {
  
  # Generate list of outputs
  allProtSets <- list(limma=list(
                                 padj=list(AllSig=allProts_padj),
                                 pval=list(AllSig=allProts_pval),
                                 All=rownames(trainQuantTable))
          )
  
  # Save
  saveRDS(allProtSets,
          here::here(resCrossDiscoDir, paste0("trainSplitProp",trainSplitProp), "sampleRandom", "differentialAbundance", 
                     paste0("chosenProteins_featureSelection_pThr",as.character(candidatePThr),"_treatLog2Thr",as.character(log2ThrVal),".rds")))
  
}
  


# # Save as RDS
# saveRDS(plotList,
#         here::here(figDir, 
#                    "06_differentialAbundance_figs.rds"))

