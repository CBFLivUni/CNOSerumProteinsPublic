

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

