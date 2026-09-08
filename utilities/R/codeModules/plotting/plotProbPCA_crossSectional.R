
# Create output directories
dir.create(here::here(resCrossDiscoDir, "eda", "probPCA", "childSamples"), recursive=TRUE, showWarnings=FALSE)

# Loop over PCs
## All together
for (i in 1:9) {
  
  # Plot
  # Technical
  plate_plt <- Plot2DPCA(meta=resPPCA_childData, var=varExpl_child, color_by="PLATE", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
    scale_fill_manual(values=plateCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Plate") +
    ggtitle("Child Samples") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  row_plt <- Plot2DPCA(meta=resPPCA_childData, var=varExpl_child, color_by="ROW", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Row") +
    ggtitle("Child Samples") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  col_plt <- Plot2DPCA(meta=resPPCA_childData, var=varExpl_child, color_by="COLUMN", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Column") +
    ggtitle("Child Samples") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  
  ## Biological
  cond_plt <- Plot2DPCA(meta=resPPCA_childData, var=varExpl_child, color_by="Condition", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
    scale_fill_manual(values=tidyCondCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Condition") +
    ggtitle("Child Samples") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  sex_plt <- Plot2DPCA(meta=resPPCA_childData, var=varExpl_child, color_by="Sex", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
    scale_fill_manual(values=sexCols) +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    labs(fill="Sex") +
    ggtitle("Child Samples") +
    guides(fill=guide_legend(ncol=1,override.aes=list(shape=21), title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  age_plt <- Plot2DPCA(meta=resPPCA_childData, var=varExpl_child, color_by="Age", label_id="Sample_ID", pcs=c(i, i+1), pch_size=7) + 
    scale_fill_continuous(high="darkred", low="white") +
    theme_classic(base_size=48) +
    theme(legend.position="right", 
          legend.direction="vertical",
          legend.title=element_text(size=64),
          legend.text=element_text(size=56),
          axis.title=element_text(size=64),
          plot.title=element_text(size=72)) +
    ggtitle("Child Samples") +
    guides(fill=guide_colourbar(barheight=10, title.position="top"), shape=guide_legend(ncol=1,title.position="top"))
  
  # Save
  ggsave(here::here(resCrossDiscoDir, "eda", "probPCA", "childSamples", paste0("probPCA_biplot.pc",i,"to",i+1,".colByPlate.crossSectional_childSamples.png")), plate_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "probPCA", "childSamples", paste0("probPCA_biplot.pc",i,"to",i+1,".colByPlateRow.crossSectional_childSamples.png")), row_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "probPCA", "childSamples", paste0("probPCA_biplot.pc",i,"to",i+1,".colByPlateColumn.crossSectional_childSamples.png")), col_plt, width=6000, height=4000, units="px")
  
  ggsave(here::here(resCrossDiscoDir, "eda", "probPCA", "childSamples", paste0("probPCA_biplot.pc",i,"to",i+1,".colByCondition.crossSectional_childSamples.png")), cond_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "probPCA", "childSamples", paste0("probPCA_biplot.pc",i,"to",i+1,".colBySex.crossSectional_childSamples.png")), sex_plt, width=6000, height=4000, units="px")
  ggsave(here::here(resCrossDiscoDir, "eda", "probPCA", "childSamples", paste0("probPCA_biplot.pc",i,"to",i+1,".colByAge.crossSectional_childSamples.png")), age_plt, width=6000, height=4000, units="px")
  
}

