

plotDetectDist <- function (se, extraLayers=list(geom_blank())) {
  
  # Function adapted from plot_detect() in the DEP package 
  
  # Dependencies
  require(radiant.data)
  
  assertthat::assert_that(inherits(se, "SummarizedExperiment"))
  se_assay <- assay(se)
  if (!any(is.na(se_assay))) {
    stop("No missing values in '", deparse(substitute(se)), 
         "'", call. = FALSE)
  }
  df <- se_assay %>% data.frame() %>% rownames_to_column() %>% 
    gather(ID, val, -rowname)
  stat <- df %>% group_by(rowname) %>% summarize(mean = mean(val, 
                                                             na.rm = TRUE), missval = any(is.na(val)))
  cumsum <- stat %>% group_by(missval) %>% arrange(mean) %>% 
    mutate(num = 1, cs = cumsum(num), cs_frac = cs/n())
  p1 <- ggplot(stat, aes(mean, col = missval)) + geom_density(na.rm = TRUE) + 
    labs(x = expression(log[2] ~ "Intensity"), y = "Density") + 
    guides(col = guide_legend(title = "Missing values")) + 
    theme_DEP1()
  p2 <- ggplot(cumsum, aes(mean, cs_frac, col = missval)) + 
    geom_line() + labs(x = expression(log[2] ~ "Intensity"), 
                       y = "Cumulative fraction") + guides(col = guide_legend(title = "Missing values")) + 
    theme_DEP1()
  
  # Set extra layers
  for (layer in extraLayers) { p1 <- p1 + layer }
  for (layer in extraLayers) { p2 <- p2 + layer }
  
  # Return
  return(list(p1, p2))
  
}


PlotAllDistributions <- function(qntTable, filts, filtValList=c(-5, 0, 0.5, 1, 1.5, 2, 2.5, 3), extraLayers=list(geom_blank())) {
  
  # This function depends on the plotDetectDist() function in the
  # plottingFunctions.R script file in utilities/
  
  # Loop over mean-filtration thresholds
  allFiltValPlts <- lapply(filtValList, function(filtVal) {
    
    print(filtVal)
    
    # Loop over different missing filtration measures
    allPlts <- lapply(filts, 
                      function(filtTable) {
                        
                        # Remove samples with all NAs
                        qntFilt <- qntTable[apply(assay(qntTable), 1, function(row) !all(is.na(row))),]
                        
                        # Filter by missingness
                        qntFilt <- qntFilt[which(rowData(qntFilt)$PG.ProteinAccessions %in% filtTable$ID),]
                        
                        # Filter by log abundnace
                        protMeans <- rowMeans(assay(qntFilt), na.rm=TRUE)
                        
                        # Filter by log abundance
                        qntFilt <- qntFilt[protMeans > filtVal,]
                        
                        # Plot intensities between samples with/without missing values
                        plt <- plotDetectDist(qntFilt, extraLayers=extraLayers)
                        
                        # Add Number of proteins
                        plt[[1]] <- plt[[1]] + ggtitle(paste0("N = ",nrow(assay(qntFilt))))
                        plt[[2]] <- plt[[2]] + ggtitle(paste0("N = ",nrow(assay(qntFilt))))
                        
                        # Sort axis
                        plt[[2]] <- plt[[2]] + ylab("Cumulative")
                        
                        # Return
                        return(plt)
                        
                      })
    
    # Set all plots except first to have no axis labesl
    for (i in 2:length(allPlts)) {
      
      allPlts[[i]][[1]] <- allPlts[[i]][[1]] + theme(axis.text.y=element_blank())
      allPlts[[i]][[2]] <- allPlts[[i]][[2]] + theme(axis.text.y=element_blank())
      
    }
    
    # Combine
    allPlt_empDist <- wrap_plots(lapply(allPlts, function(plts) plts[[1]]), ncol=6) +
      plot_layout(axis_titles = "collect")
    allPlt_cumDist <- wrap_plots(lapply(allPlts, function(plts) plts[[2]]), ncol=6) +
      plot_layout(axis_titles = "collect")
    
    # Add row number
    allPlt_empDist <- wrap_plots(wrap_elements(panel = textGrob(filtVal, gp=gpar(fontsize = 25))), allPlt_empDist) + plot_layout(widths = c(0.05, 0.95))
    allPlt_cumDist <- wrap_plots(wrap_elements(panel = textGrob(filtVal, gp=gpar(fontsize = 25))), allPlt_cumDist) + plot_layout(widths = c(0.05, 0.95))
    
    # Return
    return(list(allPlt_empDist, allPlt_cumDist))
    
  })
  
  # Tidy final plot
  empPlt <- wrap_plots(lapply(allFiltValPlts, function(plts) plts[[1]]), nrow=length(allFiltValPlts)) +   
    plot_layout(guides = "collect", axis_titles = "collect") & theme(legend.position = "bottom")
  cumPlt <- wrap_plots(lapply(allFiltValPlts, function(plts) plts[[2]]), nrow=length(allFiltValPlts)) +   
    plot_layout(guides = "collect", axis_titles = "collect") & theme(legend.position = "bottom")
  
  # Return
  return(list(EMPIRICAL_DIST=empPlt, CUMULATIVE_DIST=cumPlt))
  
}


Plot2DPCA <- function(meta, pcs=c(1,2), var=NULL, pch_size=3, color_by=NA, shape_by=NA, alpha_by=NA, label_id="SAMPLE_ID", default_color="black", reorder_colors=NA) {
  
  # Set colour factor levels
  if (all(!is.na(reorder_colors))) { meta[[color_by]] <- factor(meta[[color_by]], levels=reorder_colors)} 
  
  # If adding a separate shape
  if (is.na(shape_by)) {
    
    if (is.na(alpha_by)) {
      
      # Plot PCA plot, coloured by color-by
      pca_plt <- ggplot(meta,
                        aes(x=.data[[paste0("PC",pcs[[1]])]],
                            y=.data[[paste0("PC",pcs[[2]])]], )) + 
        geom_point(size=pch_size, alpha=0.75, color=default_color, shape=21,
                   aes(
                     fill=.data[[color_by]],
                     color=.data[[color_by]])) + 
        theme_classic(base_size=16)
      
    } else {
      
      # Plot PCA plot, coloured by color-by
      pca_plt <- ggplot(meta,
                        aes(x=.data[[paste0("PC",pcs[[1]])]],
                            y=.data[[paste0("PC",pcs[[2]])]], )) + 
        geom_point(size=pch_size, color=default_color,
                   aes(
                     alpha=.data[[alpha_by]],
                     fill=.data[[color_by]],
                     color=.data[[color_by]])) + 
        theme_classic(base_size=16)
      
    }
    
  } else {
    
    if (is.na(alpha_by)) {
      
      # Plot PCA plot, coloured by color-by
      pca_plt <- ggplot(meta,
                        aes(x=.data[[paste0("PC",pcs[[1]])]],
                            y=.data[[paste0("PC",pcs[[2]])]], )) + 
        geom_point(size=pch_size, alpha=0.75, color=default_color, shape=21,
                   aes(
                     shape=.data[[shape_by]],
                     group=.data[[shape_by]],
                     fill=.data[[color_by]],
                     color=.data[[color_by]])) + 
        theme_classic(base_size=16)
      
    } else {
      
      # Plot PCA plot, coloured by color-by
      pca_plt <- ggplot(meta,
                        aes(x=.data[[paste0("PC",pcs[[1]])]],
                            y=.data[[paste0("PC",pcs[[2]])]], )) + 
        geom_point(size=pch_size, color=default_color,
                   aes(
                     alpha=.data[[alpha_by]],
                     shape=.data[[shape_by]],
                     group=.data[[shape_by]],
                     fill=.data[[color_by]],
                     color=.data[[color_by]])) +
        theme_classic(base_size=16)
      
    }
    
  }
  
  # If we want to add variance explained to PCs
  if (!is.null(var)) { 

      # Set cell lines and axis labels
      xlabel <- paste("PC",pcs[[1]]," ~ ",round(var$variance.percent[pcs[[1]]],2),"%",sep="")
      ylabel <- paste("PC",pcs[[2]]," ~ ",round(var$variance.percent[pcs[[2]]],2),"%",sep="")
    
    # Add labels with variance explained
    pca_plt <- pca_plt + xlab(xlabel) + ylab(ylabel) 
    
  } else {
    
    # Just add "PC"
    pca_plt <- pca_plt + xlab(paste0("PC",pcs[[1]])) + ylab(paste0("PC",pcs[[2]]))
    
  }
  
  # Return
  return(pca_plt)

}



plotProtBoxplot <- function(qntTable, meta, varOfInterest, chosenCols, varOfInterest_2=NULL, refLevel="", nRow=2) {
  
  # If second var of interest unset
  if (is.null(varOfInterest_2)) { varOfInterest_2 <- varOfInterest }
  
  # Plot
  toPlot <- reshape2::melt(data.frame(FACTOR=rownames(qntTable), as.matrix(qntTable)), id.vars="FACTOR") # base::log2(
  toPlot <- merge(toPlot, meta, by.x="variable", by.y="Sample_ID")
  
  # Plot
  plt <- ggplot(toPlot[toPlot[[varOfInterest]] != refLevel,], 
                aes(x=.data[[varOfInterest]], y=value, fill=.data[[varOfInterest_2]], group=.data[[varOfInterest_2]])) +
    geom_boxplot() + 
    facet_wrap(~FACTOR, nrow=nRow, scale="free") +
    scale_fill_manual(values=chosenCols) +
    xlab(varOfInterest) + ylab("Log2 Abundance") +
    theme_classic(base_size=16) + theme(axis.title=element_text(size=28),
                                   axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
  plt
  
  # Return
  return(plt)
  
}



plotProtScatterplot <- function(qntTable, meta, varOfInterest, chosenCols, varIsDate=FALSE, varOfInterest_2=NULL, varOfInterest_3=NULL, refLevel="", nRow=2) {
  
  # Dependencies
  require(reshape2)
  require(lubridate)
  require(ggplot2)
  
  # If second var of interest unset
  if (is.null(varOfInterest_2)) { varOfInterest_2 <- varOfInterest }
  if (is.null(varOfInterest_3)) { varOfInterest_3 <- varOfInterest }
  
  # Plot
  toPlot <- reshape2::melt(data.frame(FACTOR=rownames(qntTable), as.matrix(qntTable)), id.vars="FACTOR")
  toPlot <- merge(toPlot, meta, by.x="variable", by.y="Sample_ID")
  
  # Drop NA measurements
  toPlot <- toPlot[!is.na(toPlot[[varOfInterest]]),]
  
  # Calculate correlations
  for (varLevel in unique(toPlot$FACTOR)) { 
    
    # If using date
      if (!varIsDate) {
    
        # Add correlation value to facet plots  
        toPlot[toPlot$FACTOR == varLevel,]["FACTOR"] <- paste0(toPlot[toPlot$FACTOR == varLevel,]$FACTOR, ", ρ = ",
                                                               round(cor(toPlot[toPlot$FACTOR == varLevel,][[varOfInterest]], 
                                                                         toPlot[toPlot$FACTOR == varLevel,]$value,
                                                                         use="complete.obs", method="spearman"), 2))
    
      } else {
        
        # Add correlation value to facet plots  
        toPlot[toPlot$FACTOR == varLevel,]["FACTOR"] <- paste0(toPlot[toPlot$FACTOR == varLevel,]$FACTOR, ", ρ = ",
                                                               round(cor(as.numeric(ymd(toPlot[toPlot$FACTOR == varLevel,][[varOfInterest]])), 
                                                                         toPlot[toPlot$FACTOR == varLevel,]$value,
                                                                         use="complete.obs", method="spearman"), 2))
        
      }
  }
  
  # Update object- BROKEN
  # if (varIsDate) { toPlot[[varOfInterest]] <- as.Date(dmy(ymd(toPlot[[varOfInterest]]) }
  
  # Plot
  plt <- ggplot(toPlot[toPlot[[varOfInterest]] != refLevel,], 
                aes(x=.data[[varOfInterest]], 
                    y=value, 
                    fill=.data[[varOfInterest_2]], 
                    group=.data[[varOfInterest_2]], 
                    shape=.data[[varOfInterest_3]])) +
    geom_point(size=3, alpha=0.8) + 
    facet_wrap(~FACTOR, nrow=nRow, scale="free") +
    xlab(varOfInterest) + ylab("Log2 Abundance") +
    theme_classic(base_size=16) + theme(axis.title=element_text(size=48),
                                   axis.text=element_text(size=20),
                                   strip.text=element_text(size=18))
  
  # If variable is date, update labels
  # if (varIsDate) { plt <- plt + scale_x_date(labels = date_format("%b")) }
  
  # Return
  return(plt)
  
}




plotMedianAbundances_forSelectedProts <- function(qnt, metaData, selectProts, labelProts) {
  
  # Loop over 
  plts <- lapply(unique(metaData$Condition), function(grp) {
    
    # Get group data
    grpData <- cbind(Group=grp,
                     reshape2::melt(cbind(Protein=rownames(qnt), 
                                          qnt[,metaData[metaData$Condition == grp,]$Sample_ID]), id.vars="Protein"))
    
    # Order by median abudnance
    medianAbGrp <- aggregate(value ~ Group + Protein, data=grpData[,!grepl("variable",colnames(grpData))], FUN=median)
    medianAbGrp <- medianAbGrp[order(medianAbGrp$value, decreasing=TRUE),]
    grpData$Protein <- factor(grpData$Protein, levels=unique(medianAbGrp$Protein))
    
    # Add if select proteins
    grpData["SelectedProtein"] <- "No"
    grpData[which(grpData$Protein %in% selectProts),]["SelectedProtein"] <- "Yes"
    
    # Add if label protein
    grpData["LabelledProtein"] <- "No"
    grpData[which(grpData$Protein %in% names(labelProts)),]["LabelledProtein"] <- "Yes"
    grpData["ProteinLabel"] <- NA
    grpData[grpData$LabelledProtein == "Yes",]["ProteinLabel"] <- labelProts
    
    # Plot
    plt <- ggplot(grpData[grpData$SelectedProtein == "Yes",], 
                  aes(x=Protein, 
                      y=value, 
                      fill=factor(Protein, levels=selectProts))) + 
      geom_point(data=grpData[grpData$SelectedProtein == "No",], 
                 mapping=aes(x=Protein, y=value), 
                 alpha=0.45, color="darkgrey", size=1) +
      geom_boxplot(data=grpData[grpData$SelectedProtein == "Yes",],
                   outlier.shape=NA, 
                   fill="white",
                   width=length(unique(grpData$Protein)) / 20) + 
      geom_beeswarm(color="black", 
                    size=3.5,
                    alpha=0.85) +
      # geom_text_repel(data=aggregate(value ~ Group + Protein + ProteinLabel, 
      #                          data=grpData[grpData$LabelledProtein == "Yes",!grepl("variable",colnames(grpData))], 
      #                          FUN=max),
      #           mapping=aes(x=Protein, 
      #                       y=value*1.1, 
      #                       label=ProteinLabel), 
      #           size=20, nudge_x=5) +
      scale_x_discrete(breaks = levels(grpData$Protein)[seq(0, length(unique(grpData$Protein)), by = 125)[-1]],
                       labels = seq(0, length(unique(grpData$Protein)), by = 125)[-1],
                       drop=FALSE) +
      xlab("Rank of Protein Median Abundance") + 
      ylab("Log2 Relative\nAbundance") + 
      labs(fill="Protein") + 
      ggtitle(paste0(grp)) +
      theme_classic(base_size=22) +  theme(axis.text.x=element_text(size=40),
                                           axis.text.y=element_text(size=48),
                                           axis.title.x=element_text(size=72),
                                           axis.title.y=element_text(size=72, angle=90),
                                           plot.title=element_text(size=88),
                                           legend.position="none")
    
    # Return
    return(plt)
    
  })
  
  # Return
  names(plts) <- unique(metaData$Condition)
  return(plts)
  
}


plotPerGroupVariableImportance <- function(modelSets,
                                           protData,
                                           metaData) {
  
  # Get protein IDs
  prots <- protData$Protein
  
  # Get variable importances per protein per split
  varImpData <- do.call("rbind",lapply(seq_along(modelSets), function(i) {
    
    # Split to model set
    splitSet <- modelSets[[i]]
    
    # Get local importance
    localImp <- splitSet$rfmodelRes$Model$variable.importance.local
    
    # melt
    localImpTall <- reshape2::melt(localImp)
    
    # Add groups
    localImpTall <- merge(localImpTall,
                          metaData[,c("Sample_ID", "Condition")],
                          by.x="Var1",
                          by.y="Sample_ID")
    
    # Aggregate
    localImpTall <- aggregate(value ~ Var2 + Condition,
                              data=localImpTall[,colnames(localImpTall) != "Var1"],
                              FUN=mean)
    
    # Add split
    localImpTall["Split"] <- i
    
    # Return
    return(localImpTall)
    
  }))
  
  # Add columns
  colnames(varImpData) <- c("Variable", "Condition", "VariableImportance", "Split")
  
  # Add if protein or covaraite
  varImpData["VarType"] <- "Protein"
  checkExternalVariables <- !varImpData$Variable %in% prots
  if (sum(checkExternalVariables) > 0) { varImpData[checkExternalVariables,]["VarType"] <- "Covariate" }
  
  # Drop non-protein covariates
  varImpData <- varImpData[varImpData$VarType == "Protein",]
  
  # Add 
  varImpData <- merge(varImpData,
                      protData[,c("Protein", "SYMBOL")],
                      by.x="Variable",
                      by.y="Protein")
  
  # Aggregate
  aggVarImpData <- aggregate(VariableImportance ~ SYMBOL + Condition + VarType, 
                             data=varImpData[,!grepl("Split|Variable$",colnames(varImpData))],
                             FUN=mean)
  
  # Order by CNO
  cnoVarImp <- aggVarImpData[aggVarImpData$Condition == "CNO",]
  cnoVarImp <- cnoVarImp[order(cnoVarImp$VariableImportance, decreasing=FALSE),]
  varImpData$SYMBOL <- factor(varImpData$SYMBOL,
                              cnoVarImp$SYMBOL)
  
  # Plot
  plt <- ggplot(varImpData,
                aes(y=SYMBOL,
                    x=VariableImportance,
                    color=Condition)) + 
    geom_vline(xintercept=0,
               
               color="darkgrey") +
    geom_boxplot(outlier.shape=NA,
                 color="black") +
    geom_jitter(alpha=0.65,
                size=4,
                height=0.25) +
    facet_grid(cols=vars(Condition)) +
    ylab("Variable") +
    theme_classic(base_size=30) +
    theme(legend.position="none",
          axis.title=element_text(size=36),
          strip.text=element_text(size=36))
  
  # Return
  return(list(Plt=plt,
              PermImpData=varImpData,
              MeanPermImpData=aggVarImpData))
  
}



plotOverallVariableImportance <- function(modelSets,
                                          protData) {
  
  # Get protein IDs
  prots <- protData$Protein
  
  # Get variable importances per protein per split
  varImpM <- sapply(modelSets, function(splitSet) splitSet$rfmodelRes$Model$variable.importance)
  
  # Make tall
  varImpData <- reshape2::melt(varImpM)
  colnames(varImpData) <- c("Variable", "Split", "VariableImportance")
  
  # Add if protein or covaraite
  varImpData["VarType"] <- "Protein"
  checkExternalVariables <- !varImpData$Variable %in% prots
  if (sum(checkExternalVariables) > 0) { varImpData[checkExternalVariables,]["VarType"] <- "Covariate" }
  
  # Drop non-protein covariates
  varImpData <- varImpData[varImpData$VarType == "Protein",]
  
  # Add 
  varImpData <- merge(varImpData,
                      protData[,c("Protein", "SYMBOL")],
                      by.x="Variable",
                      by.y="Protein")
  
  # Plot
  plt <- ggplot(varImpData,
                aes(y=reorder(SYMBOL,
                              VariableImportance),
                    x=VariableImportance)) + 
    geom_vline(xintercept=0,
               
               color="darkgrey") +
    geom_boxplot(outlier.shape=NA,
                 color="black") +
    geom_jitter(alpha=0.65,
                size=2,
                height=0.25,
                color="grey") +
    # facet_grid(rows=vars(VarType),
    #            scale="free",
    #            space="free") +
    xlab("Impurity Corrected") +
    ylab("Variable") +
    theme_classic(base_size=24) +
    theme(legend.position="none",
          axis.title=element_text(size=36),
          strip.text=element_text(size=36))
  
  # Return
  return(plt)
  
}



plotModelMetrics <- function(plotData) {
  
  # Get plots
  plt <- ggplot(plotData,
                aes(x=Group,
                    y=value,
                    fill=Group)) + 
    geom_hline(yintercept=0,
               
               color="darkgrey") +
    geom_boxplot() +
    scale_fill_manual(values=condCols) +
    facet_grid(cols=vars(variable)) +
    ylab("Score") + 
    xlab("") +
    scale_y_continuous(breaks=seq(0, 1, 0.2)) +
    theme_classic(base_size=16) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  # Return
  return(plt)
  
}


plotModelMetricsBar <- function(plotData,
                                cols) {
  
  # Get plots
  plt <- ggplot(plotData,
                aes(x=Group,
                    y=value,
                    fill=Group)) +
    geom_col() + 
    geom_hline(yintercept=0,
               
               color="darkgrey") +
    facet_grid(cols=vars(variable)) +
    scale_fill_manual(values=cols) +
    ylab("Score") + 
    xlab("") +
    scale_y_continuous(breaks=seq(0, 1, 0.2)) +
    theme_classic(base_size=16) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  # Return
  return(plt)
  
}


plotMisMatchSampleAbundance <- function(m, 
                                        misMatch,
                                        proteins) {
  
  # Subset matrix
  subM <- m[proteins$Protein,misMatch$ID]
  rownames(subM) <- proteins$SYMBOL
  
  # Make tall
  tallM <- reshape2::melt(cbind(Protein=rownames(subM),
                                subM),
                          id.vars="Protein")
  
  # Combine
  tallM <- merge(misMatch, 
                 tallM,
                 by.x="ID",
                 by.y="variable")
  
  # Get mismatched samples
  misMatchIDs <- misMatch[misMatch$Mismatch == "Yes",]$ID
  
  # PLots
  plts <- vector(mode="list",
                 length=length(misMatchIDs))
  
  # Re-level
  tallM$Protein <- factor(tallM$Protein,
                          levels=proteins[order(proteins$MultiModalVarImp_Validation, 
                                                decreasing=TRUE),]$SYMBOL)
  
  # Loop over samples
  for (i in seq_along(misMatchIDs)) {
    
    # Get ID
    id <- misMatchIDs[[i]]
    
    # Plot
    plt <- ggplot(tallM,
                  aes(x=Protein,
                      y=value,
                      color=Predicted)) +
      geom_boxplot(outlier.shape=NA,
                   color="black") +
      geom_jitter(data=tallM[tallM$ID != id,],
                  mapping=aes(x=Protein,
                              y=value),
                  size=2.5,
                  color="black") +
      geom_point(data=tallM[tallM$ID == id,],
                 mapping=aes(x=Protein,
                             y=value,
                             size=Mismatch),
                 size=4) +
      scale_color_manual(values=condCols) +
      labs(x="Candidate Protein",
           y="Log2 Abundance",
           title=paste0("Sample: ", id, " (Predicted as ",misMatch[misMatch$ID == id,]$Predicted,")")) +
      theme_classic(base_size=32) +
      theme(axis.text.x=element_text(angle=90,
                                     vjust=0.5,
                                     hjust=1,
                                     size=32),
            legend.text=element_text(size=36),
            legend.title=element_text(size=42),
            axis.title=element_text(size=48),
            plot.title=element_text(size=52))
    plts[[i]] <- plt
    
  }
  
  # COmbine
  bigPlt <- wrap_plots(plts,
                       ncol=1)
  
  # Return
  return(bigPlt)
  
}

