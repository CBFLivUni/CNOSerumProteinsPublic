

## Read in protein abundance table of QC tables
### Serum QC
quantTable_crossSectional_serumQC <- read.table(here::here(diaDir, "processed", "QCs", "crossSectional_pivotTable_PGQuantity.sampleIDs_justSerumQCs.tsv"),
                                                sep="\t", quote="\"", header=TRUE)

### Instrument QC
quantTable_crossSectional_instrQC <- read.table(here::here(diaDir, "processed", "QCs", "crossSectional_pivotTable_PGQuantity.justInstrumentQCs.tsv"),
                                                sep="\t", quote="\"", header=TRUE)

# Set rownames
rownames(quantTable_crossSectional_serumQC) <- quantTable_crossSectional_serumQC$PG.ProteinGroups
rownames(quantTable_crossSectional_instrQC) <- quantTable_crossSectional_instrQC$PG.ProteinGroups
