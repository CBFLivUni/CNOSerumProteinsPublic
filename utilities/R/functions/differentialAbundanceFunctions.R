


runLimma <- function(qntTable, designM, levels2Test, contrasts2Fit=NULL, returnWhat=NULL, decideTestsPVal=0.05, log2Thr=NULL) {
  
  # Dependencies
  require(limma)
  
  # Fit
  fit <- lmFit(qntTable, designM)
  
  # If specified contrasts
  if (!is.null(contrasts2Fit)) { fit <- contrasts.fit(fit, contrasts=contrasts2Fit) }
  
  # Run eBayes
  if (is.null(log2Thr)) { fit <- eBayes(fit) }
  if (!is.null(log2Thr)) { fit <- treat(fit, lfc=log2Thr) }
  
  # If we want to return the top table
  if (is.null(returnWhat)) {
    
    # Get results and tidy
    resTable <- topTable(fit, coef=levels2Test, number=Inf)
    resTable <- data.frame(PG_ID=rownames(resTable), resTable)
    
    # Return
    return(resTable)
    
    # Else if return output of "decideTests()"
  } else if (returnWhat == "decideTests") {
    
    # Decide tests
    sigHits <- decideTests(fit, p.value=decideTestsPVal)
    
    # Return
    return(sigHits)
    
  }
  
}


getDiffSign_vsChosenLevel <- function(sig, refLevel, altLevels) {
  
  # Generate list
  prots <- c()
  
  # Loop over proteins
  for (prot in rownames(sig)) {
    
    # Get sign of CNO and other levels
    cnoSign <- sign(sig[prot,][[refLevel]])
    otherSign <- unlist(sign(sig[prot,][,altLevels]))
    
    # Test
    if (all(cnoSign != otherSign)) { prots <- c(prots, prot) }
    
  }
  
  # Return
  return(prots)
  
}


TestModelAICs <- function(m, meta, fullModel, toDrop) {
  
  # Dependencies
  require(limma)
  
  # Drop all models
  testDesigns <- list(model.matrix(as.formula(paste0(" ~ ",paste(fullModel, collapse=" + "))), data=meta))
  
  # Loop over variables
  for (i in seq_len(length(toDrop))) {
    
    # Get interaction terms
    interactTerms <- fullModel[grepl("[:]",fullModel)]
    
    # If dropping simple term, remove equivalent interaction term
    if (any(grepl(toDrop[[i]], interactTerms))) { 
      
      # Add to design
      subFormula <- as.formula(paste0(" ~ ",paste(setdiff(fullModel, c(toDrop[[i]], interactTerms[grepl(toDrop[[i]],interactTerms)])), collapse=" + ")))
      testDesigns[i+1] <- list(model.matrix(subFormula, data=meta))
      
      # Else just drop interaction term
    } else {
      
      # Add to design
      subFormula <- as.formula(paste0(" ~ ",paste(setdiff(fullModel, toDrop[[i]]), collapse=" + ")))
      testDesigns[i+1] <- list(model.matrix(subFormula, data=meta))
      
    }
    
  }
  names(testDesigns) <- c("Full",toDrop)
  
  # Choose models
  modelAICs <- selectModel(m, testDesigns)
  print(table(modelAICs$pref))
  
  # Get difference between model and base model
  modelAICs$IC <- as.matrix(modelAICs$IC)
  modelAICs <- rlist::list.append(modelAICs, deltaIC=as.matrix(modelAICs$IC))
  modelAICs$deltaIC <- modelAICs$IC[,1] - modelAICs$IC
  
  # Return
  return(modelAICs)
  
}

runProDA <- function(quantTable, designFormula, meta, refLevel, outDir, extraFileName, lrtFTest=FALSE, contrasts2Test=NULL, onlyReadInFile=FALSE, reGenFiles=TRUE) {
  
  # Dependencies
  require(proDA)
  
  # If we want to only read in the file
  if (onlyReadInFile & !reGenerateFiles) {
    
    proDAFit <- NULL
  
  # Else set proDA fit to null
  } else {
    
    # Fit model
    print("///// RUNNING ProDA")
    proDAFit <- proDA(as.matrix(quantTable), design = as.formula(designFormula), 
                      col_data = meta, reference_level = refLevel)
    
  }
  
  # If contrasts are NULL, set to all results names
  if (is.null(contrasts2Test)) { contrasts2Test <- result_names(proDAFit) }
    
  # If we want to do an LRT F test
  if (lrtFTest) {
    
    # Perform LRT test on chosen factor
    print(paste0("// Performing LRT test for :- ", paste(contrasts2Test, collapse=", ")))
    proDA_resTable <- performLRTTests_forProDA(proDAFit, lrtContrasts=contrasts2Test, fullFormula=as.formula(designFormula),
                                               reGenFile=reGenFiles, outDir=outDir, fileExtra=extraFileName)
    
  } else {
    
    # Get results per contrast
    proDA_resTable <- data.frame()
    for (contr in contrasts2Test) {
      
      # Perform test
      print(paste0("// Running individual Wald tests for :- ",contr))
      proDA_res <- test_diff(proDAFit, contrast=contr)
      proDA_resTable  <- rbind(proDA_resTable, 
                               data.frame(CONTRAST=contr, proDA_res))
      
    }
    
  }
  
  # Return
  return(proDA_resTable)
  
}



performLRTTests_forProDA <- function(fitObj, lrtContrasts, fullFormula, outDir="./", fileExtra="", reGenFile=TRUE) {
  
  # Fit LRT models
  proDA_lrtResTable <- data.frame()
  for (contr in lrtContrasts) {
    
    # Report
    print(paste0("///// Working on :- ", contr))
    
    # Save
    filename <- here::here(outDir, gsub("[.][.]","",paste0("proDAResults_lrtTest",gsub("[:]","-",contr),".",fileExtra,".tsv")))
    if (file.exists(filename) & !reGenFile) { 
      
      # Read in from file
      print(paste0("// Reading in from file: ",filename))
      proDA_testRes <- read.table(filename, sep="\t", header=TRUE, row.names=NULL)
      proDA_lrtResTable <- rbind(proDA_lrtResTable, 
                                 data.frame(CONTRAST=contr, proDA_testRes)) 
      
    } else {
      
      # Generate reduced formula
      reducedFormula <- gsub(paste0(contr, " "), " ", as.character(fullFormula))
      reducedFormula <- gsub("[+]\\s*[+]", "+", paste(reducedFormula, collapse=" "))
      reducedFormula <- gsub("[~]\\s*[+]","~",reducedFormula)
      print(paste0("// Dropping :- ", contr))
      print(paste0("// Testing with reduced formula: ",reducedFormula))
      
      # Get sig 
      proDA_testRes <- test_diff(fitObj, 
                                 reduced_model = as.formula(reducedFormula))
      proDA_lrtResTable <- rbind(proDA_lrtResTable, 
                                 data.frame(CONTRAST=contr, proDA_testRes)) 
      
      # Save
      print(paste0("// Saving to file: ",filename))
      write.table(proDA_testRes, filename, sep="\t", quote=FALSE, row.names=FALSE)
      
    }
    
  }
  
  # Return
  return(proDA_lrtResTable)
  
}


plotAIC <- function(aicData,
                    dropTerm=NULL) {
  
  # Get tall
  tallAIC <- reshape2::melt(aicData$IC[,-1], id.vars=c("ModelSet","Protein"))
  
  # Get
  tallAIC["Chosen"] <- "No"
  if (!is.null(dropTerm)) { tallAIC[which(tallAIC$Models %in% dropTerm),]["Chosen"] <- "Dropped" }
  
  # Plot
  p1 <- ggplot(tallAIC, 
               aes(x=Models, 
                   y=value,
                   fill=Chosen)) + 
    geom_hline(yintercept=0, linetype="dashed", color="darkgrey") +
    geom_point(size=2, alpha=0.15, shape=21) + 
    stat_summary(fun.y = median, geom = "errorbar", 
                 aes(ymax = ..y.., ymin = ..y.., group = Models),
                 width = 1, linetype = "solid") +
    xlab("Variable being Dropped") + ylab("AIC") +
    scale_fill_manual(values=c("No"="darkgrey","Dropped"="skyblue")) +
    theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  p2 <- ggplot(tallAIC, 
               aes(x=Models, 
                   y=value,
                   fill=Chosen)) +
    geom_hline(yintercept=0, linetype="dashed", color="darkgrey") +
    geom_boxplot() +
    xlab("Variable being Dropped") + ylab("AIC") +
    scale_fill_manual(values=c("No"="darkgrey","Dropped"="skyblue")) +
    theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  
  # Return
  return(p1 + p2)
  
}


GetLmerSeqContrasts <- function(fit, contr_m) {
  
  # Loop over rows of contrast matrix
  out_fits <- lapply(colnames(contr_m), function(contr_name) {
    
    # Get results
    print(paste0("##### Fitting contrast :- ",contr_name))
    contr <- lmerSeq.contrast(fit, contrast=rbind(t(contr_m)[contr_name,]), include_singular=T, sort_results = T)
    
    # Return
    return(contr)
    
  }); names(out_fits) <- colnames(contr_m)
  
  # Return
  return(out_fits)
  
}

GetLmerSeqCoef <- function(fit, coefs) {
  
  # Loop over rows of contrast matrix
  out_fits <- lapply(coefs, function(coef_name) {
    
    # Get results
    print(paste0("##### Extracting coefficient for :- ",coef_name))
    summary <- lmerSeq.summary(fit, coefficient=coef_name, include_singular=T, sort_results = T)
    
    # Return
    return(summary)
    
  }); names(out_fits) <- coefs
  
  # Return
  return(out_fits)
  
}

