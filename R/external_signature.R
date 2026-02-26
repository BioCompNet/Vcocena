
make_rownames_unique <- function(counts_new){

  counts_new <- counts_new[!duplicated(counts_new$SYMBOL),] %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., "SYMBOL")
  # remove all non-numeric columns (description, gene id, etc.):
  for (x in colnames(counts_new)){
    if(!is.numeric(counts_new[[x]])){
      counts_new[[x]] <- NULL
    }
  }
  if(ncol(counts_new) == 0){
    print("All columns were deleted when removing non-numeric columns. 
          Please check the data type of your expression values.")
  }

  return(counts_new)
}



# Gene names will be used as rownames, please provide a column names "SYMBOL" !
read_expression_data <- function(file, rown = T, sep = "\t"){
  if(rown){
    expression_data <- read.table(file = file, row.names = 1,
                       stringsAsFactors = F, sep = sep, check.names = F)
  }else{
    expression_data <- read.table(file = file, 
                       stringsAsFactors = F, sep = sep, check.names = F)
  }
  
  expression_data <- make_rownames_unique(expression_data)
  return(expression_data)
}

# Sample IDs will be used as rownames, please provide a column names "SampleID" !
read_anno <- function(file, rown = T, sep = "\t"){
  if(rown){
    anno <- read.table(file = file, row.names = 1,
                                  stringsAsFactors = F, sep = sep, check.names = T) 
    
  }else{
    anno <- read.table(file = file, 
                                  stringsAsFactors = F, sep = sep, check.names = T)
  }
  
  rownames(anno) <- anno$SampleID
  return(anno)
}


calculate_GFCs <- function(expressions, anno, genes, grpvar, range_GFC){
  
  # filter expression data for genes of which the GFC shall be calculated
  expressions <- tibble::rownames_to_column(expressions, var = "SYMBOL")%>%
    dplyr::filter(., SYMBOL %in% genes)%>%
    tibble::column_to_rownames(., "SYMBOL")
  # calculate GFCs for groups
  if(!grpvar %in% colnames(anno)){
    print("The column you specified to contain the grouping variables does not exists in the annotation table.")
    return()
  }
  
  GFCs <- NULL
  for (g in genes){
    g_df <- expressions[rownames(expressions) == g,]
    colnames(g_df) <- colnames(expressions)
    
    # vector of group means:
    mean_vec <- NULL
    for (var in unique(anno[[grpvar]])){
      var_anno <- anno[anno[[grpvar]] == var,]
      var_exp <- g_df[, colnames(g_df) %in% var_anno$ID]
      if(nrow(var_anno) == 1){
        mean_vec <- c(mean_vec,mean(var_exp %>% as.numeric()))
      }else{
        mean_vec <- c(mean_vec,mean(var_exp[1,] %>% as.numeric()))
      }
      
      
    }
    # overall expression mean of current gene:
    g_mean <- mean(mean_vec)
    GFC_vec <- NULL
    for (x in 1:length(mean_vec)){
      GFC_vec <- c(GFC_vec, gtools::foldchange(mean_vec[x], g_mean) %>% 
                     ifelse(.>range_GFC, range_GFC,.) %>%
                     ifelse(.< (-range_GFC), -range_GFC,.))
    }
    # add this genes GFC vector to GFCs
    GFCs <- rbind( GFCs, GFC_vec)
  }
  GFCs <- as.data.frame(GFCs) %>%
    round(., digits = 3)
  # add column with gene names
  GFCs$Gene <- genes
  colnames(GFCs) <- c(unique(anno[[grpvar]]), "Gene") 
  rownames(GFCs) <- GFCs$"Gene"
  return(GFCs)
}



plot_cluster_heatmap <- function(cluster_info = integrated_output$cluster_calc$cluster_information, 
                                 GFCs = integrated_output$GFC_all_layers, 
                                 col_order = NULL, user_enrichment_1 = NULL, 
                                 user_enrichment_2 = NULL, row_order = NULL, cluster_columns, 
                                 cluster_rows, k = 0, column_anno_categorical = NULL, column_anno_numerical = NULL,
                                 return_HM = T){
  gc()
  # filter for included clusters (non-white)
  c_df <- dplyr::filter(cluster_info, cluster_included == "yes")
  mat_heatmap <- NULL
  
  
  if(!is.null(row_order)){
    
    rownames(GFCs) <- GFCs$cluster
    GFCs$cluster <- NULL
    c_GFC_means <- GFCs[row_order,]
    
    mat_heatmap <- rbind(mat_heatmap, c_GFC_means)
    
    rownames(mat_heatmap) <- row_order
    
  }else{
    rownames(GFCs) <- GFCs$cluster
    GFCs$cluster <- NULL
    c_GFC_means <- GFCs[c_df$color,]
    
    mat_heatmap <- rbind(mat_heatmap, c_GFC_means)
    
    rownames(mat_heatmap) <- c_df$color
  }
  
  
  colnames(mat_heatmap) <- colnames(GFCs)
  
  if(!is.null(col_order)){
    mat_heatmap <- mat_heatmap %>% as.data.frame()
    mat_heatmap <- mat_heatmap[, col_order] %>% as.matrix()
  }
  
  enrich_mat1 <- list()
  enrich_count1 <- list()
  enrich_mat2 <- list()
  enrich_count2 <- list()
  
  if(!is.null(row_order)){
    if(!is.null(user_enrichment_1)){
      for(x in row_order){
        enrich_mat1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x)%>%
          dplyr::pull(., count)
        enrich_count1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x)%>%
          dplyr::pull(., hits)%>%
          dplyr::first(.)
      }
    }
    if(!is.null(user_enrichment_2)){
      for(x in row_order){
        enrich_mat2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x)%>%
          dplyr::pull(., count)
        enrich_count2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x)%>%
          dplyr::pull(., hits)%>%
          dplyr::first(.)
      }
    }
    
  }else{
    for(x in unique(user_enrichment_1$cluster)){
      enrich_mat1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x)%>%
        dplyr::pull(., count)
      enrich_count1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x)%>%
        dplyr::pull(., hits)%>%
        dplyr::first(.)
    }
    for(x in unique(user_enrichment_2$cluster)){
      enrich_mat2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x)%>%
        dplyr::pull(., count)
      enrich_count2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x)%>%
        dplyr::pull(., hits)%>%
        dplyr::first(.)
    }
  }
  
  if(!length(enrich_mat1) == 0){
    enrich_mat1 <- matrix(unlist(enrich_mat1), nrow =length(enrich_mat1), byrow = T)
    enrich_count1 <- unlist(enrich_count1)
    
  }
  if(all(enrich_mat1 == 0) == T){
    enrich_mat1 <- list()
  }
  
  if(!length(enrich_mat2) == 0){
    
    enrich_mat2 <- matrix(unlist(enrich_mat2), nrow =length(enrich_mat2), byrow = T)
    enrich_count2 <- unlist(enrich_count2)
  }
  if(all(enrich_mat2 == 0) == T){
    enrich_mat2 <- list()
  }
  
  if(!is.null(row_order)){
    cluster_colors <- factor(row_order)
    names(cluster_colors) <- row_order
    c_df <- c_df[match(rowvec, c_df$color),]
  }else{
    cluster_colors <- factor(c_df$color)
    names(cluster_colors) <- c_df$color
    row_order <- unique(c_df$color)
  }
  
  
  if(length(enrich_mat1) == 0 & length(enrich_mat2) == 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = anno_simple(row_order, col = cluster_colors, 
                                                                  simple_anno_size = unit(0.5, "cm")),
                                            entities = anno_barplot(c_df$gene_no, width = unit(2.5, "cm")),
                                            
                                            
                                            which = "row", 
                                            width = unit(4.5, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"), 
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    
    lgd_list <- list(
      
    )
  }
  else if(length(enrich_mat1) > 0 & length(enrich_mat2) == 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = anno_simple(row_order, col = cluster_colors, 
                                                                  simple_anno_size = unit(0.5, "cm")),
                                            entities = anno_barplot(c_df$gene_no, width = unit(2.5, "cm")),
                                            
                                            enriched_count = anno_text(paste0(enrich_count1, "/", c_df$gene_no), width = unit(1.5, "cm")),
                                            enriched = anno_barplot(enrich_mat1,
                                                                    width = unit(3, "cm"),
                                                                    gp = gpar(fill = RColorBrewer::brewer.pal(n = 12, name = "Paired"),
                                                                              col = RColorBrewer::brewer.pal(n = 12, name = "Paired"))),
                                            which = "row", 
                                            width = unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"), 
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    
    lgd_list <- list(
      
      ComplexHeatmap::Legend(labels = unique(user_enrichment_1$cell_type), title = "enriched",
                             legend_gp = gpar(col = RColorBrewer::brewer.pal(n = 12, name = "Paired")),
                             type = "points", pch = 15)
    )
  }
  else if(length(enrich_mat1) == 0 & length(enrich_mat2) > 0){
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = anno_simple(row_order, col = cluster_colors, 
                                                                  simple_anno_size = unit(0.5, "cm")),
                                            entities = anno_barplot(c_df$gene_no, width = unit(1.5, "cm")),
                                            
                                            enriched_count = anno_text(paste0(enrich_count2, "/", c_df$gene_no), width = unit(1.5, "cm"),
                                                                       gp = gpar(fontsize = 8)),
                                            enriched = anno_barplot(enrich_mat2,
                                                                    width = unit(5, "cm"),
                                                                    gp = gpar(fill = ggsci::pal_d3(palette = "category20")(ncol(enrich_mat2)),
                                                                              col = ggsci::pal_d3(palette = "category20")(ncol(enrich_mat2)))),
                                            
                                            which = "row", 
                                            width = unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"), 
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    
    lgd_list <- list(
      
      
      ComplexHeatmap::Legend(labels = unique(user_enrichment_2$cell_type), title = "enriched",
                             legend_gp = gpar(col = ggsci::pal_d3(palette = "category20")(20)),
                             type = "points", pch = 15)
    )
  }else{
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = anno_simple(row_order, col = cluster_colors, 
                                                                  simple_anno_size = unit(0.25, "cm")),
                                            entities = anno_barplot(c_df$gene_no, width = unit(0.75, "cm")),
                                            
                                            enriched_count_1 = anno_text(c_df$gene_no,
                                                                         width = unit(0.75, "cm"),
                                                                         gp = gpar(fontsize = 8)),
                                            enriched_1 = anno_barplot(enrich_mat1,
                                                                      width = unit(2, "cm"),
                                                                      gp = gpar(fill = RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:ncol(enrich_mat1)],
                                                                                col = RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:ncol(enrich_mat1)]),
                                                                      baseline = 0),
                                            enriched_count_2 = anno_text(paste0(enrich_count2, "/", c_df$gene_no), 
                                                                         width = unit(0.75, "cm"),
                                                                         gp = gpar(fontsize = 8)),
                                            enriched_2 = anno_barplot(enrich_mat2,
                                                                      width = unit(2, "cm"),
                                                                      gp = gpar(fill = ggsci::pal_d3(palette = "category20")(20)[1:ncol(enrich_mat2)],
                                                                                col = ggsci::pal_d3(palette = "category20")(20)[1:ncol(enrich_mat2)]),
                                                                      baseline = 0),
                                            which = "row", 
                                            width = unit(12, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"), 
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    
    lgd_list <- list(
      
      
      ComplexHeatmap::Legend(labels = unique(user_enrichment_1$cell_type), title = "enriched_1",
                             legend_gp = gpar(col =  RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:ncol(enrich_mat1)]),
                             type = "points", pch = 15),
      ComplexHeatmap::Legend(labels = unique(user_enrichment_2$cell_type), title = "enriched_2",
                             legend_gp = gpar(col = ggsci::pal_d3(palette = "category20")(20)[1:ncol(enrich_mat2)]),
                             type = "points", pch = 15)
    )
  }
  
  
  anno_list <- NULL
  
  if(!length(column_anno_categorical) == 0){
    for(a in 1:length(column_anno_categorical)){
      set.seed(a)
      tmp_colour <- ggsci::pal_startrek()(7)[sample(1:7, ncol(column_anno_categorical[[a]]))]
      anno_list <- anno_list %v% ComplexHeatmap::HeatmapAnnotation(col_anno = anno_lines(column_anno_categorical[[a]]%>%as.matrix(),width = unit(2, "cm"),
                                                                                         gp = gpar(col = tmp_colour),
                                                                                         add_points = TRUE, 
                                                                                         pt_gp = gpar(col = tmp_colour), pch = 16),
                                                                   which = "column",
                                                                   height = unit(1, "cm"),
                                                                   annotation_name_side = "right",
                                                                   gap = unit(2, "mm"),
                                                                   annotation_name_rot = 0,
                                                                   annotation_name_gp = gpar(fontsize = 8),
                                                                   annotation_label = names(column_anno_categorical)[a])
      
      lgd_list <- rlist::list.append(lgd_list, ComplexHeatmap::Legend(labels = colnames(column_anno_categorical[[a]]%>%as.matrix()), title = names(column_anno_categorical)[a],
                                                                      legend_gp = gpar(col = tmp_colour),
                                                                      type = "points", pch = 15))
    }
    
  }
  
  
  
  if(!length(column_anno_numerical) == 0){
    for(a in 1:length(column_anno_numerical)){
      tmp_col_anno_2 <- column_anno_numerical[[a]]
      tmp_col_anno_2 <- tmp_col_anno_2[colnames(mat_heatmap)]
      anno_list <- anno_list %v% ComplexHeatmap::HeatmapAnnotation(cont_anno = anno_boxplot(tmp_col_anno_2, height = unit(1, "cm")), 
                                                                   which = "column",
                                                                   annotation_name_side = "right",
                                                                   gap = unit(2, "mm"),
                                                                   annotation_name_rot = 0,
                                                                   annotation_name_gp = gpar(fontsize = 8),
                                                                   annotation_label = names(column_anno_numerical)[a], show_legend = F)
      
    }
  }
  
  
  all_conditions <- NULL
  
  if(!is.null(column_anno_categorical) | !is.null(column_anno_numerical)){
    for(setnum in 1:length(layers)){
      all_conditions <- c(all_conditions, as.character(dplyr::pull(data[[paste0("set", setnum, "_anno")]], global_settings$voi)))
    }
    all_conditions <- table(all_conditions) %>% 
      as.data.frame()%>%
      dplyr::filter(., all_conditions %in% colnames(mat_heatmap))
    all_conditions <- all_conditions[match(colnames(mat_heatmap), as.character(all_conditions$all_conditions)),]
    all_conditions <- paste0(all_conditions$all_conditions, "  [", all_conditions$Freq, "]")
    anno_list <- anno_list %v% ComplexHeatmap::columnAnnotation(groups = anno_text(all_conditions))
    
  }
  
  
  
  
  Cairo(file = paste0(working_directory,global_settings$save_folder, "/module_heatmap.pdf"),
        width = 50, 
        height = 30,
        pointsize=11,
        dpi=300,
        type = "pdf",
        units = "in")

  # Define breaks across the full GFC range and generate a matching colour palette
  breaks <- seq(-global_settings$range_GFC, global_settings$range_GFC, by = .1)
  my_colour <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(breaks))
  col_fun <- colorRamp2(breaks, my_colour)
  #
  hm <- ComplexHeatmap::Heatmap(mat_heatmap %>% as.matrix(), 
                                right_annotation = ha,
                                col = col_fun,
                                clustering_distance_rows = "euclidean",
                                clustering_distance_columns = "euclidean",
                                clustering_method_rows = "complete",
                                clustering_method_columns = "complete",
                                cluster_columns = cluster_columns,
                                cluster_rows = cluster_rows,
                                column_names_rot = 90,
                                column_names_centered = F,
                                row_names_gp = gpar(fontsize = 8),
                                column_names_gp = gpar(fontsize = 8),
                                rect_gp = grid::gpar(col = "black"),
                                heatmap_legend_param = list( title = "", legend_height = unit(4, "cm")), column_km = k)
  
  hm_w_lgd <- ComplexHeatmap::draw(hm %v% anno_list, annotation_legend_list = lgd_list, merge_legends = T, 
                                   padding = unit(c(2, 2, 2, 30), "mm"))
  
  dev.off()
  
  if(return_HM){
    return(hm_w_lgd)
  }
  
}


#' Plot integrated cluster heatmap (v2)
#'
#' Draws a module heatmap for integrated clusters, optionally with module labels
#' and annotations, and saves to a file if \code{output_file} is provided.
#'
#' @param cluster_info Cluster information data frame.
#' @param GFCs Matrix/data frame of GFC values per cluster.
#' @param cluster_columns Logical; cluster columns.
#' @param cluster_rows Logical; cluster rows.
#' @param clip_quantiles Numeric vector length 2 to clip color scale.
#' @param palette Palette name (e.g., \code{"BlueRed"}).
#' @param legend_title Legend title for the heatmap.
#' @param output_file Optional path to save the heatmap (PDF).
#'
#' @return A \code{ComplexHeatmap} object when \code{return_HM = TRUE}.
#' @export
plot_cluster_heatmap_v2 <- function(cluster_info = integrated_output$cluster_calc$cluster_information,
                                    GFCs = integrated_output$GFC_all_layers,
                                    col_order = NULL,
                                    user_enrichment_1 = NULL,
                                    user_enrichment_2 = NULL,
                                    row_order = NULL,
                                    cluster_columns,
                                    cluster_rows,
                                    k = 0,
                                    column_anno_categorical = NULL,
                                    column_anno_numerical = NULL,
                                    return_HM = TRUE,
                                    range_GFC = global_settings$range_GFC,
                                    clip_quantiles = c(0.01, 0.99),
                                    palette = c("RdBu", "BlueRed", "custom"),
                                    custom_palette = NULL,
                                    legend_title = "GFC",
                                    show_row_names = TRUE,
                                    show_column_names = TRUE,
                                    row_name_size = 12,
                                    column_name_size = 12,
                                    row_labels_mode = c("name", "hex"),
                                    row_labels = NULL,
                                    modules_only = FALSE,
                                    module_label_prefix = "M",
                                    label_modules_in_bar = FALSE,
                                    module_label_color = "white",
                                    module_bar_width_cm = 1.0,
                                    display = FALSE,
                                    row_split = NULL,
                                    column_split = NULL,
                                    cat_anno_colors = NULL,
                                    num_anno_legend = TRUE,
                                    output_file = NULL,
                                    width_in = NULL,
                                    height_in = NULL,
                                    dpi = 300,
                                    na_col = "grey90"){
  gc()
  palette <- match.arg(palette)
  row_labels_mode <- match.arg(row_labels_mode)
  c_df <- dplyr::filter(cluster_info, cluster_included == "yes")
  mat_heatmap <- NULL

  if (!is.null(row_order)) {
    rownames(GFCs) <- GFCs$cluster
    GFCs$cluster <- NULL
    missing_clusters <- setdiff(row_order, rownames(GFCs))
    if (length(missing_clusters) > 0) {
      warning("Missing clusters in GFCs: ", paste(missing_clusters, collapse = ", "))
    }
    row_order <- row_order[row_order %in% rownames(GFCs)]
    c_GFC_means <- GFCs[row_order, ]
    mat_heatmap <- rbind(mat_heatmap, c_GFC_means)
    rownames(mat_heatmap) <- row_order
  } else {
    rownames(GFCs) <- GFCs$cluster
    GFCs$cluster <- NULL
    c_GFC_means <- GFCs[c_df$color, ]
    mat_heatmap <- rbind(mat_heatmap, c_GFC_means)
    rownames(mat_heatmap) <- c_df$color
  }

  colnames(mat_heatmap) <- colnames(GFCs)

  if (!is.null(col_order)) {
    mat_heatmap <- mat_heatmap %>% as.data.frame()
    mat_heatmap <- mat_heatmap[, col_order, drop = FALSE] %>% as.matrix()
  }

  enrich_mat1 <- list()
  enrich_count1 <- list()
  enrich_mat2 <- list()
  enrich_count2 <- list()

  if (!is.null(row_order)) {
    if (!is.null(user_enrichment_1)) {
      for (x in row_order) {
        enrich_mat1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
          dplyr::pull(., count)
        enrich_count1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
          dplyr::pull(., hits) %>%
          dplyr::first(.)
      }
    }
    if (!is.null(user_enrichment_2)) {
      for (x in row_order) {
        enrich_mat2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
          dplyr::pull(., count)
        enrich_count2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
          dplyr::pull(., hits) %>%
          dplyr::first(.)
      }
    }
  } else {
    if (!is.null(user_enrichment_1)) {
      for (x in unique(user_enrichment_1$cluster)) {
        enrich_mat1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
          dplyr::pull(., count)
        enrich_count1[[x]] <- dplyr::filter(user_enrichment_1, cluster == x) %>%
          dplyr::pull(., hits) %>%
          dplyr::first(.)
      }
    }
    if (!is.null(user_enrichment_2)) {
      for (x in unique(user_enrichment_2$cluster)) {
        enrich_mat2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
          dplyr::pull(., count)
        enrich_count2[[x]] <- dplyr::filter(user_enrichment_2, cluster == x) %>%
          dplyr::pull(., hits) %>%
          dplyr::first(.)
      }
    }
  }

  if (!length(enrich_mat1) == 0) {
    enrich_mat1 <- matrix(unlist(enrich_mat1), nrow = length(enrich_mat1), byrow = TRUE)
    enrich_count1 <- unlist(enrich_count1)
  }
  if (!length(enrich_mat1) == 0 && all(enrich_mat1 == 0) == TRUE) {
    enrich_mat1 <- list()
  }

  if (!length(enrich_mat2) == 0) {
    enrich_mat2 <- matrix(unlist(enrich_mat2), nrow = length(enrich_mat2), byrow = TRUE)
    enrich_count2 <- unlist(enrich_count2)
  }
  if (!length(enrich_mat2) == 0 && all(enrich_mat2 == 0) == TRUE) {
    enrich_mat2 <- list()
  }

  if (!is.null(row_order)) {
    cluster_colors <- factor(row_order)
    names(cluster_colors) <- row_order
    c_df <- c_df[match(row_order, c_df$color), ]
  } else {
    cluster_colors <- factor(c_df$color)
    names(cluster_colors) <- c_df$color
    row_order <- unique(c_df$color)
  }

  resolve_color_names <- function(hex_colors) {
    named_cols <- grDevices::colors()
    named_rgb <- grDevices::col2rgb(named_cols)
    out <- character(length(hex_colors))
    for (i in seq_along(hex_colors)) {
      col_i <- hex_colors[i]
      if (is.na(col_i) || !nzchar(col_i)) {
        out[i] <- col_i
        next
      }
      rgb <- tryCatch(grDevices::col2rgb(col_i), error = function(e) NULL)
      if (is.null(rgb)) {
        out[i] <- col_i
        next
      }
      d <- colSums(sweep(named_rgb, 1, rgb[, 1], "-") ^ 2)
      out[i] <- named_cols[which.min(d)]
    }
    out
  }

  if (is.null(row_labels)) {
    if (row_labels_mode == "name") {
      row_labels <- resolve_color_names(row_order)
    } else {
      row_labels <- row_order
    }
  }

  if (modules_only) {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      modules = anno_simple(row_order, col = cluster_colors, simple_anno_size = unit(0.5, "cm")),
      which = "row",
      width = unit(module_bar_width_cm, "cm"),
      annotation_name_side = "top",
      gap = unit(2, "mm"),
      annotation_name_rot = 0,
      annotation_name_gp = gpar(fontsize = 8)
    )
    lgd_list <- list()
    if (label_modules_in_bar) {
      show_row_names <- FALSE
      row_labels <- NULL
    }
  } else if (length(enrich_mat1) == 0 && length(enrich_mat2) == 0) {
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = anno_simple(row_order, col = cluster_colors,
                                                                  simple_anno_size = unit(0.5, "cm")),
                                            entities = anno_barplot(c_df$gene_no, width = unit(2.5, "cm")),
                                            which = "row",
                                            width = unit(4.5, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    lgd_list <- list()
  } else if (length(enrich_mat1) > 0 && length(enrich_mat2) == 0) {
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = anno_simple(row_order, col = cluster_colors,
                                                                  simple_anno_size = unit(0.5, "cm")),
                                            entities = anno_barplot(c_df$gene_no, width = unit(2.5, "cm")),
                                            enriched_count = anno_text(paste0(enrich_count1, "/", c_df$gene_no),
                                                                       width = unit(1.5, "cm")),
                                            enriched = anno_barplot(enrich_mat1,
                                                                    width = unit(3, "cm"),
                                                                    gp = gpar(fill = RColorBrewer::brewer.pal(n = 12, name = "Paired"),
                                                                              col = RColorBrewer::brewer.pal(n = 12, name = "Paired"))),
                                            which = "row",
                                            width = unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    lgd_list <- list(
      ComplexHeatmap::Legend(labels = unique(user_enrichment_1$cell_type), title = "enriched",
                             legend_gp = gpar(col = RColorBrewer::brewer.pal(n = 12, name = "Paired")),
                             type = "points", pch = 15)
    )
  } else if (length(enrich_mat1) == 0 && length(enrich_mat2) > 0) {
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = anno_simple(row_order, col = cluster_colors,
                                                                  simple_anno_size = unit(0.5, "cm")),
                                            entities = anno_barplot(c_df$gene_no, width = unit(1.5, "cm")),
                                            enriched_count = anno_text(paste0(enrich_count2, "/", c_df$gene_no),
                                                                       width = unit(1.5, "cm"),
                                                                       gp = gpar(fontsize = 8)),
                                            enriched = anno_barplot(enrich_mat2,
                                                                    width = unit(5, "cm"),
                                                                    gp = gpar(fill = ggsci::pal_d3(palette = "category20")(ncol(enrich_mat2)),
                                                                              col = ggsci::pal_d3(palette = "category20")(ncol(enrich_mat2)))),
                                            which = "row",
                                            width = unit(9, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    lgd_list <- list(
      ComplexHeatmap::Legend(labels = unique(user_enrichment_2$cell_type), title = "enriched",
                             legend_gp = gpar(col = ggsci::pal_d3(palette = "category20")(20)),
                             type = "points", pch = 15)
    )
  } else {
    ha <- ComplexHeatmap::HeatmapAnnotation(modules = anno_simple(row_order, col = cluster_colors,
                                                                  simple_anno_size = unit(0.25, "cm")),
                                            entities = anno_barplot(c_df$gene_no, width = unit(0.75, "cm")),
                                            enriched_count_1 = anno_text(c_df$gene_no,
                                                                         width = unit(0.75, "cm"),
                                                                         gp = gpar(fontsize = 8)),
                                            enriched_1 = anno_barplot(enrich_mat1,
                                                                      width = unit(2, "cm"),
                                                                      gp = gpar(fill = RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:ncol(enrich_mat1)],
                                                                                col = RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:ncol(enrich_mat1)]),
                                                                      baseline = 0),
                                            enriched_count_2 = anno_text(paste0(enrich_count2, "/", c_df$gene_no),
                                                                         width = unit(0.75, "cm"),
                                                                         gp = gpar(fontsize = 8)),
                                            enriched_2 = anno_barplot(enrich_mat2,
                                                                      width = unit(2, "cm"),
                                                                      gp = gpar(fill = ggsci::pal_d3(palette = "category20")(20)[1:ncol(enrich_mat2)],
                                                                                col = ggsci::pal_d3(palette = "category20")(20)[1:ncol(enrich_mat2)]),
                                                                      baseline = 0),
                                            which = "row",
                                            width = unit(12, "cm"),
                                            annotation_name_side = "top",
                                            gap = unit(2, "mm"),
                                            annotation_name_rot = 0,
                                            annotation_name_gp = gpar(fontsize = 8))
    lgd_list <- list(
      ComplexHeatmap::Legend(labels = unique(user_enrichment_1$cell_type), title = "enriched_1",
                             legend_gp = gpar(col = RColorBrewer::brewer.pal(n = 12, name = "Paired")[1:ncol(enrich_mat1)]),
                             type = "points", pch = 15),
      ComplexHeatmap::Legend(labels = unique(user_enrichment_2$cell_type), title = "enriched_2",
                             legend_gp = gpar(col = ggsci::pal_d3(palette = "category20")(20)[1:ncol(enrich_mat2)]),
                             type = "points", pch = 15)
    )
  }

  anno_list <- NULL

  if (!length(column_anno_categorical) == 0) {
    for (a in 1:length(column_anno_categorical)) {
      anno_cols <- column_anno_categorical[[a]] %>% as.matrix()
      if (nrow(anno_cols) == 0 || ncol(anno_cols) == 0) {
        next
      }
      if (!is.null(cat_anno_colors) && !is.null(cat_anno_colors[[a]])) {
        tmp_colour <- cat_anno_colors[[a]]
      } else {
        base_pal <- RColorBrewer::brewer.pal(3, "Set2")
        if (ncol(anno_cols) <= 3) {
          tmp_colour <- base_pal[seq_len(ncol(anno_cols))]
        } else {
          tmp_colour <- colorRampPalette(base_pal)(ncol(anno_cols))
        }
      }
      anno_list <- anno_list %v% ComplexHeatmap::HeatmapAnnotation(col_anno = anno_lines(anno_cols,
                                                                                         width = unit(2, "cm"),
                                                                                         gp = gpar(col = tmp_colour),
                                                                                         add_points = TRUE,
                                                                                         pt_gp = gpar(col = tmp_colour), pch = 16),
                                                                   which = "column",
                                                                   height = unit(1, "cm"),
                                                                   annotation_name_side = "right",
                                                                   gap = unit(2, "mm"),
                                                                   annotation_name_rot = 0,
                                                                   annotation_name_gp = gpar(fontsize = 8),
                                                                   annotation_label = names(column_anno_categorical)[a])
      lgd_list <- rlist::list.append(lgd_list, ComplexHeatmap::Legend(labels = colnames(anno_cols),
                                                                      title = names(column_anno_categorical)[a],
                                                                      legend_gp = gpar(col = tmp_colour),
                                                                      type = "points", pch = 15))
    }
  }

  if (!length(column_anno_numerical) == 0) {
    for (a in 1:length(column_anno_numerical)) {
      tmp_col_anno_2 <- column_anno_numerical[[a]]
      tmp_col_anno_2 <- tmp_col_anno_2[colnames(mat_heatmap)]
      if (length(tmp_col_anno_2) == 0) {
        next
      }
      anno_list <- anno_list %v% ComplexHeatmap::HeatmapAnnotation(cont_anno = anno_boxplot(tmp_col_anno_2,
                                                                                           height = unit(1, "cm")),
                                                                   which = "column",
                                                                   annotation_name_side = "right",
                                                                   gap = unit(2, "mm"),
                                                                   annotation_name_rot = 0,
                                                                   annotation_name_gp = gpar(fontsize = 8),
                                                                   annotation_label = names(column_anno_numerical)[a],
                                                                   show_legend = FALSE)
      if (num_anno_legend) {
        lgd_list <- rlist::list.append(lgd_list, ComplexHeatmap::Legend(labels = c("min", "max"),
                                                                        title = names(column_anno_numerical)[a],
                                                                        legend_gp = gpar(col = "black"),
                                                                        type = "points", pch = 15))
      }
    }
  }

  all_conditions <- NULL
  if (!is.null(column_anno_categorical) | !is.null(column_anno_numerical)) {
    for (setnum in 1:length(layers)) {
      all_conditions <- c(all_conditions, as.character(dplyr::pull(data[[paste0("set", setnum, "_anno")]], global_settings$voi)))
    }
    all_conditions <- table(all_conditions) %>%
      as.data.frame() %>%
      dplyr::filter(., all_conditions %in% colnames(mat_heatmap))
    all_conditions <- all_conditions[match(colnames(mat_heatmap), as.character(all_conditions$all_conditions)), ]
    all_conditions <- paste0(all_conditions$all_conditions, "  [", all_conditions$Freq, "]")
    anno_list <- anno_list %v% ComplexHeatmap::columnAnnotation(groups = anno_text(all_conditions))
  }

  if (is.null(output_file)) {
    output_file <- paste0(working_directory, global_settings$save_folder, "/module_heatmap_v2.pdf")
  }

  if (is.null(width_in) || is.null(height_in)) {
    n_cols <- ncol(mat_heatmap)
    n_rows <- nrow(mat_heatmap)
    if (is.null(width_in)) {
      width_in <- max(8, 0.5 * n_cols + 4)
    }
    if (is.null(height_in)) {
      height_in <- max(6, 0.35 * n_rows + 4)
    }
  }

  if (requireNamespace("Cairo", quietly = TRUE)) {
    Cairo::Cairo(file = output_file, width = width_in, height = height_in, pointsize = 11,
                 dpi = dpi, type = "pdf", units = "in")
  } else {
    pdf(file = output_file, width = width_in, height = height_in, pointsize = 11)
  }

  mat_heatmap_mat <- mat_heatmap %>% as.matrix()
  mat_values <- as.numeric(mat_heatmap_mat)
  mat_values <- mat_values[is.finite(mat_values)]
  if (length(mat_values) > 0 && all(is.finite(clip_quantiles))) {
    q <- quantile(mat_values, probs = clip_quantiles, na.rm = TRUE)
    min_v <- q[[1]]
    max_v <- q[[2]]
    if (!is.finite(min_v) || !is.finite(max_v) || min_v >= max_v) {
      min_v <- -range_GFC
      max_v <- range_GFC
    }
  } else {
    min_v <- -range_GFC
    max_v <- range_GFC
  }

  if (palette == "custom" && !is.null(custom_palette)) {
    pal <- custom_palette
  } else if (palette == "BlueRed") {
    pal <- colorspace::divergingx_hcl(11, palette = "RdBu")
  } else {
    pal <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  }
  my_colour <- colorRampPalette(pal)(101)
  col_fun <- colorRamp2(seq(min_v, max_v, length.out = length(my_colour)), my_colour)

  hm <- ComplexHeatmap::Heatmap(mat_heatmap_mat,
                                right_annotation = ha,
                                col = col_fun,
                                na_col = na_col,
                                clustering_distance_rows = "euclidean",
                                clustering_distance_columns = "euclidean",
                                clustering_method_rows = "complete",
                                clustering_method_columns = "complete",
                                cluster_columns = cluster_columns,
                                cluster_rows = cluster_rows,
                                column_names_rot = 90,
                                column_names_centered = FALSE,
                                show_row_names = show_row_names,
                                show_column_names = show_column_names,
                                row_labels = row_labels,
                                row_names_gp = gpar(fontsize = row_name_size),
                                column_names_gp = gpar(fontsize = column_name_size),
                                rect_gp = grid::gpar(col = "black"),
                                heatmap_legend_param = list(title = legend_title, legend_height = unit(4, "cm")),
                                column_km = k,
                                row_split = row_split,
                                column_split = column_split)

  ht <- ComplexHeatmap::draw(hm %v% anno_list, annotation_legend_list = lgd_list, merge_legends = TRUE,
                             padding = unit(c(2, 2, 2, 30), "mm"))

  if (modules_only && label_modules_in_bar) {
    module_labels <- paste0(module_label_prefix, seq_along(row_order))
    names(module_labels) <- row_order
    row_ord_list <- ComplexHeatmap::row_order(ht)
    if (is.list(row_ord_list)) {
      row_ord <- unlist(row_ord_list)
    } else {
      row_ord <- row_ord_list
    }
    row_names <- rownames(mat_heatmap_mat)
    row_names_ord <- row_names[row_ord]
    module_labels_ord <- module_labels[row_names_ord]
    n_rows <- length(module_labels_ord)
    ComplexHeatmap::decorate_annotation("modules", {
      grid::grid.text(module_labels_ord,
                      x = unit(0.5, "npc"),
                      y = unit((n_rows:1 - 0.5) / n_rows, "npc"),
                      gp = gpar(col = module_label_color, fontsize = row_name_size))
    })
  }

  dev.off()

  if (display) {
    ComplexHeatmap::draw(hm %v% anno_list, annotation_legend_list = lgd_list, merge_legends = TRUE,
                         padding = unit(c(2, 2, 2, 30), "mm"))
  }

  if (return_HM) {
    return(ht)
  }
}




entity_enrichment <- function(cluster_info){
  
  clusters <- unique(cluster_info$color)
  clusters <- clusters[!clusters == "white"]
  
  f <- list()
  for(i in 1:length(layers)){
    # tmp <- data.frame(X1 = rownames(data[[paste0("set", i, "_counts")]]))
    # colnames(tmp) <- layers_names[i]
    # f[[i]] <- tmp
    f[[layers_names[i]]] <- rownames(data[[paste0("set", i, "_counts")]])
  }
  # f <- dplyr::bind_rows(f) %>% as.data.frame()
  enrichment_keys <- names(f) #colnames
  categories_per_cluster <- NULL
  
  for(c in clusters){
    
    genes <- dplyr::filter(cluster_info, color == c)%>%
      dplyr::pull(., "gene_n")%>%
      base::strsplit(., split = ",")%>%
      BiocGenerics::unlist(.) 
    
    cell_enrich <- list(counts = list(), genes = list())
    hits <- 0
    
    for (type in enrichment_keys){
      # tmp <- genes[genes %in% f[, c(type)]] 
      # print(head(f[[c(type)]]))
      tmp <- genes[genes %in% f[[c(type)]]] 
      hits <- hits + length(tmp)
      cell_enrich$counts[[type]] <- length(tmp)
      cell_enrich$genes[[type]] <- tmp
    }
    
    if(hits > length(genes)){
      hits <- length(genes)
    }
    
    tmp <- data.frame(matrix(unlist(cell_enrich$counts), ncol= length(cell_enrich$counts), byrow=T) %>% t(),stringsAsFactors=FALSE)%>%
      cbind(., names(cell_enrich$counts))%>%
      cbind(., rep(c, length(cell_enrich$counts)))
    
    sink(file = paste0(working_directory, global_settings$save_folder, "/enrichedInKeys_",c,".txt"))
    for (n in names(cell_enrich$genes)){
      cat(NULL, sep = "\n")
      cat(n, sep = "\n\n")
      cat(cell_enrich$genes[[n]], sep = "\n")
    }
    
    sink()
    colnames(tmp) <- c("count", "cell_type", "cluster")
    if(hits == 0){
      tmp$count <- 0
    }else{
      tmp$count <- (tmp$count/hits)*100
    }
    
    tmp$hits <- rep(hits, nrow(tmp))
    
    
    categories_per_cluster <- rbind(categories_per_cluster, tmp)
    
  }
  return(categories_per_cluster)
}





compare_external_signature <- function(sample_file,
                                       anno_file,
                                       grpvar,
                                       cluster_info,
                                       range_GFC){
  
  
  c_df <- dplyr::filter(cluster_info, cluster_included == "yes")
  genes <- NULL
  for (c in unique(c_df$color)){
    n <- length(genes)
    #get genes from the original cluster
    genes <- c(genes, c_df[c_df$color == c, ] %>%
      dplyr::pull(., "gene_n") %>%
      base::strsplit(., split = ",") %>%
      unlist(.))
    print(c)
    print(length(genes)-n)
  }
  tmp_GFCS <- calculate_GFCs(sample_file, anno_file, genes =  genes, grpvar = grpvar, range_GFC = range_GFC)
  return(tmp_GFCS)
  
}
