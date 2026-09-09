
# Read in pivot table
## Read in protein abundance table of actual samples
quantTable_crossSectional <- read.table(here::here(diaDir, "processed", 
                                                   "crossSectional_pivotTable_PGQuantity.sampleIDs_noSerumQCs.tsv"), sep="\t", quote="\"", header=TRUE)