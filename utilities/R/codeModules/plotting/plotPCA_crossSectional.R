
# Create output directories
dir.create(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples"), recursive=TRUE, showWarnings=FALSE)

# Loop over PCs
## All together
for (i in 1:9) {
  
  # Plot
  # Technical
  plate_plt <- Plot2DPCA(meta=crossChildPCA_metaData, var=get_eigenvalue(crossChildPCA), color_by="PLATE", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
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
  row_plt <- Plot2DPCA(meta=crossChildPCA_metaData, var=get_eigenvalue(crossChildPCA), color_by="ROW", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Row") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  col_plt <- Plot2DPCA(meta=crossChildPCA_metaData, var=get_eigenvalue(crossChildPCA), color_by="COLUMN", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Column") +
    guides(fill=guide_legend(ncol=1,title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  site_plt <- Plot2DPCA(meta=crossChildPCA_metaData, var=get_eigenvalue(crossChildPCA), shape_by="Age_Group", color_by="Site", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
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
  
  ## Biological
  cond_plt <- Plot2DPCA(meta=crossChildPCA_metaData, var=get_eigenvalue(crossChildPCA), color_by="Condition", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
    scale_fill_manual(values=tidyCondCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(shape="Age Group") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  condSite_plt <- Plot2DPCA(meta=crossChildPCA_metaData, var=get_eigenvalue(crossChildPCA), shape_by="Site", color_by="Condition", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
    scale_shape_manual(values=21:25) + 
    scale_fill_manual(values=tidyCondCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(shape="Site") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  sex_plt <- Plot2DPCA(meta=crossChildPCA_metaData, var=get_eigenvalue(crossChildPCA), color_by="Sex", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
    scale_fill_manual(values=sexCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    guides(fill=guide_legend(ncol=1,title.position="top", override.aes=list(shape=21)), shape=guide_legend(ncol=1,title.position="top"))
  age_plt <- Plot2DPCA(meta=crossChildPCA_metaData, var=get_eigenvalue(crossChildPCA), color_by="Age", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
    scale_fill_continuous(high="darkred", low="white") +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    guides(fill=guide_colourbar(barheight=10, title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  
  # Save
  ggsave(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples", paste0("PCA_biplot.pc",i,"to",i+1,".colByPlate.onlyNonMissingProts_crossSectional_childSamples.png")), plate_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples", paste0("PCA_biplot.pc",i,"to",i+1,".colByPlateRow.onlyNonMissingProts_crossSectional_childSamples.png")), row_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples", paste0("PCA_biplot.pc",i,"to",i+1,".colByPlateColumn.onlyNonMissingProts_crossSectional_childSamples.png")), col_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples", paste0("PCA_biplot.pc",i,"to",i+1,".colBySite.onlyNonMissingProts_crossSectional_childSamples.png")), site_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples", paste0("PCA_biplot.pc",i,"to",i+1,".colByCondition.onlyNonMissingProts_crossSectional_childSamples.png")), cond_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples", paste0("PCA_biplot.pc",i,"to",i+1,".colByCondition_shapeBySite.onlyNonMissingProts_crossSectional_childSamples.png")), 
         condSite_plt, width=6000, height=6000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples", paste0("PCA_biplot.pc",i,"to",i+1,".colByCondition_shapeByFacetBySite.onlyNonMissingProts_crossSectional_childSamples.png")), 
         condSite_plt + facet_wrap(~Site, ncol=2), width=6000, height=6000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples", paste0("PCA_biplot.pc",i,"to",i+1,".colBySex.onlyNonMissingProts_crossSectional_childSamples.png")), sex_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "PCA", "childSamples", paste0("PCA_biplot.pc",i,"to",i+1,".colByAge.onlyNonMissingProts_crossSectional_childSamples.png")), age_plt, width=6000, height=4000, units="px")
  
}

