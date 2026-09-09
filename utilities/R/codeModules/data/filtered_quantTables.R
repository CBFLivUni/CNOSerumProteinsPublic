
# Save cross-Sectional
quantTable_crossSectional <- read.table(here::here(diaDir, "filtered", paste0("crossSectional_quantTable.filtBy_serumQCProts_maxMissing",maxPropMissing,".tsv")), sep="\t", header=TRUE, row.names=1)
