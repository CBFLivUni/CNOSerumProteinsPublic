
# Heatmap
crossPlt <- Heatmap(quantTableCross_naZeros[,metaCrossSectional$Sample_ID],
                    
                    # Names
                    # name = "Child",
                    
                    # Colours
                    col = mainColFun,
                    
                    # Legend direction
                    heatmap_legend_param = list(title = "Log2 Abundance",
                                                direction = "horizontal",
                                                legend_width = unit(10, "cm"),
                                                grid_height = unit(1.25, "cm"),
                                                labels_gp = gpar(fontsize=38),
                                                legend_gp = gpar(fontsize=72)),
                    show_heatmap_legend =FALSE,
                    
                    # Clustering
                    clustering_distance_rows = "manhattan",
                    clustering_distance_columns = "manhattan",
                    
                    # Columns
                    column_title_rot = 65,
                    column_split = metaCrossSectional[colnames(quantTableCross_naZeros[,metaCrossSectional$Sample_ID]),]$PLATE,
                    top_annotation = HeatmapAnnotation(Condition = anno_simple(metaCrossSectional$Condition,
                                                                               simple_anno_size = unit(2, "cm"), col=condCols),
                                                       Age = anno_simple(metaCrossSectional$Age, col = ageColFun_cross,
                                                                         simple_anno_size = unit(2, "cm")),
                                                       Cohort = anno_simple(metaCrossSectional$Age_Group, col=c("Adult"="darkred","Child"="darkblue"),
                                                                            simple_anno_size = unit(2, "cm")),
                                                       Sex = anno_simple(metaCrossSectional$Sex, col=sexCols,
                                                                         simple_anno_size = unit(2, "cm")),
                                                       Plate = anno_simple(as.character(metaCrossSectional$PLATE), col=plateCols,
                                                                           simple_anno_size = unit(2, "cm")),
                                                       Site = anno_simple(as.character(metaCrossSectional$Site), col=siteCols,
                                                                          simple_anno_size = unit(2, "cm")),
                                                       annotation_name_gp= gpar(fontsize = 48)),
                    column_title_gp = gpar(fontsize = 42, fontface = "bold"),
                    
                    # Rows
                    show_row_names = FALSE,
                    border = TRUE
)

png(file=here::here(resCrossDiscoDir, "eda", "impCheck_naZeroHmap_splitByPlate_crossSectional.png"), units="px", width=2500, height=2500)
draw(crossPlt, ht_gap = unit(5, "cm"), heatmap_legend_side = "bottom", legend_title_gp = gpar(fontsize=72),
     column_title="Cross Sectional (Split by Plate)")
dev.off()

# Heatmap
crossPlt <- Heatmap(quantTableCross_naZeros[,metaCrossSectional$Sample_ID],
                    
                    # Names
                    # name = "Child",
                    
                    # Colours
                    col = mainColFun,
                    
                    # Legend direction
                    heatmap_legend_param = list(title = "Log2 Abundance",
                                                direction = "horizontal",
                                                legend_width = unit(10, "cm"),
                                                grid_height = unit(1.25, "cm"),
                                                labels_gp = gpar(fontsize=38),
                                                legend_gp = gpar(fontsize=72)),
                    show_heatmap_legend =FALSE,
                    
                    # Clustering
                    clustering_distance_rows = "manhattan",
                    clustering_distance_columns = "manhattan",
                    
                    # Columns
                    column_title_rot = 65,
                    column_split = metaCrossSectional[colnames(quantTableCross_naZeros[,metaCrossSectional$Sample_ID]),]$Condition,
                    top_annotation = HeatmapAnnotation(Condition = anno_simple(metaCrossSectional$Condition,
                                                                               simple_anno_size = unit(2, "cm"), col=condCols),
                                                       Age = anno_simple(metaCrossSectional$Age, col = ageColFun_cross,
                                                                         simple_anno_size = unit(2, "cm")),
                                                       Cohort = anno_simple(metaCrossSectional$Age_Group, col=c("Adult"="darkred","Child"="darkblue"),
                                                                            simple_anno_size = unit(2, "cm")),
                                                       Sex = anno_simple(metaCrossSectional$Sex, col=sexCols,
                                                                         simple_anno_size = unit(2, "cm")),
                                                       Plate = anno_simple(as.character(metaCrossSectional$PLATE), col=plateCols,
                                                                           simple_anno_size = unit(2, "cm")),
                                                       Site = anno_simple(as.character(metaCrossSectional$Site), col=siteCols,
                                                                          simple_anno_size = unit(2, "cm")),
                                                       annotation_name_gp= gpar(fontsize = 48)),
                    column_title_gp = gpar(fontsize = 42, fontface = "bold"),
                    
                    # Rows
                    show_row_names = FALSE,
                    border = TRUE
)

png(file=here::here(resCrossDiscoDir, "eda", "impCheck_naZeroHmap_splitByCondition_crossSectional.png"), units="px", width=2500, height=2500)
draw(crossPlt, ht_gap = unit(5, "cm"), heatmap_legend_side = "bottom", legend_title_gp = gpar(fontsize=72),
     column_title="Cross Sectional (Split by Condition)")
dev.off()

# Heatmap
crossPlt <- Heatmap(quantTableCross_naZeros[,metaCrossSectional$Sample_ID],
                    
                    # Names
                    # name = "Child",
                    
                    # Colours
                    col = mainColFun,
                    
                    # Legend direction
                    heatmap_legend_param = list(title = "Log2 Abundance",
                                                direction = "horizontal",
                                                legend_width = unit(10, "cm"),
                                                grid_height = unit(1.25, "cm"),
                                                labels_gp = gpar(fontsize=38),
                                                legend_gp = gpar(fontsize=72)),
                    show_heatmap_legend =FALSE,
                    
                    # Clustering
                    clustering_distance_rows = "manhattan",
                    clustering_distance_columns = "manhattan",
                    
                    # Columns
                    column_title_rot = 65,
                    top_annotation = HeatmapAnnotation(Condition = anno_simple(metaCrossSectional$Condition,
                                                                               simple_anno_size = unit(2, "cm"), col=condCols),
                                                       Age = anno_simple(metaCrossSectional$Age, col = ageColFun_cross,
                                                                         simple_anno_size = unit(2, "cm")),
                                                       Cohort = anno_simple(metaCrossSectional$Age_Group, col=c("Adult"="darkred","Child"="darkblue"),
                                                                            simple_anno_size = unit(2, "cm")),
                                                       Sex = anno_simple(metaCrossSectional$Sex, col=sexCols,
                                                                         simple_anno_size = unit(2, "cm")),
                                                       Plate = anno_simple(as.character(metaCrossSectional$PLATE), col=plateCols,
                                                                           simple_anno_size = unit(2, "cm")),
                                                       Site = anno_simple(as.character(metaCrossSectional$Site), col=siteCols,
                                                                          simple_anno_size = unit(2, "cm")),
                                                       annotation_name_gp= gpar(fontsize = 48)),
                    column_title_gp = gpar(fontsize = 42, fontface = "bold"),
                    
                    # Rows
                    show_row_names = FALSE,
                    border = TRUE
)

png(file=here::here(resCrossDiscoDir, "eda", "impCheck_naZeroHmap_crossSectional.png"), units="px", width=2500, height=2500)
draw(crossPlt, ht_gap = unit(5, "cm"), heatmap_legend_side = "bottom", legend_title_gp = gpar(fontsize=72),
     column_title="Cross Sectiona")
dev.off()
