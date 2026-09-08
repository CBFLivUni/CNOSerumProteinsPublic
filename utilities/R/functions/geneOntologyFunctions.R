


runEnrichmentAnalysis <- function(protList, backgroundProts, orgDbObj, term2Gene, ...) {
  
  require(clusterProfiler)
  require(org.Hs.eg.db)
  
  # Convert UniProt to Entrez ID using bitr function
  entzProts  <- bitr(protList, 
                     fromType = "UNIPROT", 
                     toType = "ENTREZID", 
                     OrgDb = orgDbObj)
  
  # Get background universe
  backGround <- bitr(backgroundProts, 
                     fromType = "UNIPROT", 
                     toType = "ENTREZID", 
                     OrgDb = orgDbObj)
  
  # Make term2gene character
  term2Gene$term <- as.character(term2Gene$term)
  term2Gene$gene <- as.character(term2Gene$gene)
  
  # Enrichment analysis
  enrichRes <- enricher(as.character(entzProts$ENTREZID),
                        pvalueCutoff=1,
                        pAdjustMethod="BH",
                        universe=as.character(backGround$ENTREZID),
                        TERM2GENE=term2Gene,
                        ...)
  
  # Get results
  return(enrichRes)
  
}

