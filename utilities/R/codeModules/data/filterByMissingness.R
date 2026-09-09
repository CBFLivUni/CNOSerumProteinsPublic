

rownames(quantTable_crossSectional) <- quantTable_crossSectional$PG.ProteinGroups
dim(quantTable_crossSectional)

# Drop proteins with high frequencies of missing values
crossPropProts <- EnumerateNAs(quantTable_crossSectional[,metaCrossSectional$Sample_ID], idx=1)

quantTable_crossSectional <- quantTable_crossSectional[crossPropProts[crossPropProts$PROP_MISSING < maxPropMissing,]$ID,]

# Define samples with high frequencies of missing values
crossPropSmpls <- EnumerateNAs(quantTable_crossSectional[,metaCrossSectional$Sample_ID], idx=2)
crossPropSmpls <- merge(crossPropSmpls, metaCrossSectional, by.x="ID", by.y="Sample_ID")

quantTable_crossSectional <- data.frame(PG.ProteinAccessions=rownames(quantTable_crossSectional), 
                                        quantTable_crossSectional[,c(crossPropSmpls[crossPropSmpls$PROP_MISSING < maxPropMissing,]$ID)])

# Subset metadata
metaCrossSectional <- metaCrossSectional[which(metaCrossSectional$Sample_ID %in% colnames(quantTable_crossSectional)),]

# Filter proteins by abundance
# protMeans <- rowMeans(log2(quantTable[,metaCrossSectional$Sample_ID]), na.rm=TRUE)
# quantTable <- quantTable[protMeans > 5,]

# Check max missingness
if (max(EnumerateNAs(quantTable_crossSectional[,metaCrossSectional$Sample_ID], idx=1)$PROP_MISSING) > maxPropMissing) { warning("Protein Missingness for Cross-Sectional not controlled at desired level") } 
if (max(EnumerateNAs(quantTable_crossSectional[,metaCrossSectional$Sample_ID], idx=2)$PROP_MISSING) > maxPropMissing) { warning("Sample Missingness for Cross-Sectional not controlled at desired level") } 

# Filter proteins again to make sure proteins are < 20% after removing samples
crossPropProts <- EnumerateNAs(quantTable_crossSectional[,metaCrossSectional$Sample_ID], idx=1)
quantTable_crossSectional <- data.frame(PG.ProteinAccessions=rownames(quantTable_crossSectional), 
                                        quantTable_crossSectional[,c(crossPropSmpls[crossPropSmpls$PROP_MISSING < maxPropMissing,]$ID)])
dim(quantTable_crossSectional)
