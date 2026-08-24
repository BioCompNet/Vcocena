heatmap_top_var <- function(info_dataset, global_settings, 
                            working_directory, input, x, plot_HM, anno_groups = global_settings$voi){

  print(paste0("...creating heatmap of top variable entities..."))
  output <- list()
  
  filt_cutoff_data = input$corr_calc_out$correlation_df_filt %>% dplyr::filter(., rval > cutoff_vec[x])
  filt_cutoff_graph = igraph::graph_from_data_frame(filt_cutoff_data,directed=FALSE)
  filt_cutoff_counts = input$ds[row.names(input$ds) %in% names(V(filt_cutoff_graph)),]
  corresp_info = info_dataset[rownames(t(input$topvar))%in%rownames(info_dataset),]
  output[["corresp_info"]] <- corresp_info
  output[["filt_cutoff_graph"]] <- filt_cutoff_graph
  output[["filt_cutoff_data"]] <- filt_cutoff_data
  output[["used_cutoff"]] <- cutoff_vec[x]
  
  print(paste("After using the optimal cutoff of ",cutoff_vec[x], " the number of edges = ", 
              nrow(filt_cutoff_data), " and the number of nodes = ", nrow(filt_cutoff_counts)))
  #rownames(info_dataset) <- info_dataset$SampleID
  #print(colnames(filt_cutoff_counts))
  #print(rownames(info_dataset))
  #rownames(info_dataset) <- info_dataset$ID %>% as.character()
  all_anno_groups <- list()
  mat_colors <- list()
  if (!is.null(anno_groups) && length(anno_groups) > 0) {
    for(i in anno_groups){
      all_anno_groups[[i]] <- as.character(dplyr::pull(info_dataset, i))
      if(length(unique(all_anno_groups[[i]])) > 20){
        mat_colors[[i]] <- colorRampPalette(ggsci::pal_d3("category20")(20))(length(unique(all_anno_groups[[i]])))
      }else{
        mat_colors[[i]] <- ggsci::pal_d3("category20")(20)[1:length(unique(all_anno_groups[[i]]))]
      }
      names(mat_colors[[i]]) <- sort(unique(all_anno_groups[[i]]))
    }
  } else {
    anno_groups <- character(0)
  }
  # unique_anno_groups <- unique(all_anno_groups)
  
  # print(mat_colors)
  # print("unique:")
  # print(unique_anno_groups)
  
  # if(length(unique_anno_groups) > 7){
  #   
  #   mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(unique_anno_groups))
  #   # mat_colors <- list(group = rainbow(length(unique(dplyr::pull(info_dataset, global_settings$voi)))))
  #   mat_colors <- list(group = mycolors)
  #   names(mat_colors$group) <- sort(unique_anno_groups)
  #   names(mat_colors) <- all_anno_groups
  # }else{
  #   mat_colors <- list(group = ggsci::pal_jama("default")(length(unique_anno_groups)))
  #   
  #   names(mat_colors$group) <- sort(unique_anno_groups)
  #   names(mat_colors) <- all_anno_groups
  # }
  
  if(plot_HM){
    heatmap_filtered_counts <- pheatmap::pheatmap(
      mat = input$topvar,  # filt_cutoff_counts
      color = rev(RColorBrewer::brewer.pal(11, "RdBu")),
      scale = "row",
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      annotation_colors = if (length(anno_groups) > 0) mat_colors else NULL,
      annotation_col = if (length(anno_groups) > 0) info_dataset[anno_groups] else NULL,
      fontsize = 8,
      show_rownames = FALSE, 
      show_colnames = TRUE
    )
  }else{
    print("don't plot")
    pdf(file = NULL)
    heatmap_filtered_counts <- pheatmap::pheatmap(mat = filt_cutoff_counts ,
                                                  color=rev(RColorBrewer::brewer.pal(11, "RdBu")),
                                                  scale="row",
                                                  cluster_rows=T,
                                                  cluster_cols=T,
                                                  annotation_colors = mat_colors,
                                                  annotation_col=info_dataset[global_settings$voi],
                                                  fontsize = 8,
                                                  show_rownames = F, 
                                                  show_colnames = F)
    dev.off()
  }
  
  output[["heatmap"]] <- heatmap_filtered_counts
  ggsave(filename = paste0("Heatmap_topvar_genes",x,".pdf"), plot = heatmap_filtered_counts$gtable, device = cairo_pdf,
         path = paste0(working_directory,global_settings$save_folder), width = 7, height = 10, units = "in", limitsize = F)
  return(output)
}
