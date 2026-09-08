
performPVCA <- function(m, meta, batchFactors, method_key="", pcaThr=0.7) {
  
  # Expressionset data
  expr <- ExpressionSet(assayData = as.matrix(m[,rownames(meta)]),
                        phenoData = new("AnnotatedDataFrame", data = meta))
  
  # Run PVCA
  res <- pvcaBatchAssess(abatch = expr, batch.factors = batchFactors, threshold = pcaThr)
  
  # Get as df
  df <- data.frame(as.data.frame(res$label), t(as.data.frame(res$dat)), Method=method_key)
  colnames(df) <- c("Effect", "Variance", "Method")
  
  # Add percentage variance
  df["PercVar"] <- round(df$Variance * 100, 2)
  df["VarThrPCA"] <- pcaThr
  
  # Return
  return(df)
  
}