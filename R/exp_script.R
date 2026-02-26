
run_expression_analysis_1 <- function(x){

  print(paste0("Currently processed dataset: ", layers_names[x]))
  output <- list()
  
  count_table <- data[[paste0("set", x, "_counts")]]
  
  if(ncol(count_table) == 0){
    return(output)
  }
  
  anno_table <- data[[paste0("set", x, "_anno")]]
  #print(anno_table)
  # mart <- biomaRt::useMart("ensembl")
  # mart <- biomaRt::listDatasets(mart)
  # output[["human"]] <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
  # output[["mouse"]] <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl") 
  # if(base::tolower(global_settings$organism) == "human"){
  #   output[["universe_Entrez"]] <- clusterProfiler:: bitr(row.names(count_table), 
  #                                                         fromType="SYMBOL", 
  #                                                         toType="ENTREZID", 
  #                                                         OrgDb="org.Hs.eg.db", 
  #                                                         drop = T)
  # }
  # if(base::tolower(global_settings$organism) == "mouse"){
  #   output[["universe_Entrez"]] <- clusterProfiler:: bitr(row.names(count_table), 
  #                                                         fromType="SYMBOL", 
  #                                                         toType="ENTREZID", 
  #                                                         OrgDb="org.Mm.eg.db", 
  #                                                         drop = T)
  # }
  # 
  # 
  # 
  if(!layer_settings[[paste0("set", x)]][["top_var"]] == "all"){
    print(paste0("...extracting ", layer_settings[[paste0("set", x)]][["top_var"]], " top variant genes..."))
  }
  
  
  # merge counts and controls for variance calculation:
  # count_and_ctrl <- merge(count_table, data$control, by = 0, all = F)%>%
  #   tibble::column_to_rownames(., "Row.names")
  #count_and_ctrl <- count_table
  
  #print(paste0("before removing controls: ", ncol(count_and_ctrl), " columns, of which ", ncol(data$control), " are controls."))
  # extract condition vector:
  #cond <- c(anno_table[[global_settings$voi]], data$control_anno[[global_settings$voi]])
  #cond <- c(anno_table[[global_settings$voi]])
  # calculate variance between conditions:
  # count_and_ctrl$var <- t(apply(count_and_ctrl , 1 , function(x) tapply(x , cond , mean)))%>%
  #   apply(., 1, var)
  
  
  #ds <- count_and_ctrl[order(count_and_ctrl$var, decreasing = T), colnames(count_table)]
  
  # ds <- count_and_ctrl[order(count_and_ctrl$var, decreasing = T), ]
  # ds$var <- NULL
  
  
  #print(paste0("after removing controls: ", ncol(ds), " columns"))
  #ds = count_and_ctrl[order(apply(count_and_ctrl,1,var), decreasing=T),]
  ds = count_table[order(apply(count_table,1,var), decreasing=T),]
  output[["ds"]] <- ds
  if(!layer_settings[[paste0("set", x)]][["top_var"]] == "all"){
    dd2 <- head(ds, layer_settings[[paste0("set", x)]][["top_var"]])
  }else{
    dd2 <- ds
  }
  
  output[["topvar"]] <- dd2
  dd2 = t(dd2)

  corr_calc_out <- correlation_actions(dd2 = dd2, layer_set = layer_settings[[paste0("set", x)]])
  output[["corr_calc_out"]] <- corr_calc_out
  
  print(paste0("...calculating cutoff statistics..."))
  cutoff_stats = do.call("rbind", lapply(X = corr_calc_out$range_cutoff,
                                         FUN = cutoff_prep,
                                         corrdf_r = corr_calc_out$correlation_df_filt,
                                         print.all.plots = layer_settings[[paste0("set", x)]]$print_distribution_plots,
                                         global_set = global_settings,
                                         x = x))
  
  
  output[["cutoff_stats"]] <- cutoff_stats
  
  cutoff_calc_out <- optcut_fun(cutoff_stats = cutoff_stats)
  
  output[["cutoff_calc_out"]] <- cutoff_calc_out
  knitr::kable(cutoff_calc_out$cutoff_stats_concise, caption = "Correlation cut-off stats")
  #plot(opt_cut_out$dd_plot_calculated_optimal)
  # ggsave(filename = paste0("Degree_distribution_plot_gene_expression_", opt_cut_out$optimal_cutoff,"_set",x, ".pdf"), 
  #        plot = opt_cut_out$dd_plot_calculated_optimal, 
  #        device = cairo_pdf, path = paste0(working_directory,global_settings$save_folder))
  
  
  
  return(output)
}


run_expression_analysis_2 <- function(x, grouping_v, plot_HM = T, anno_groups = global_settings$voi){
  
  print(paste0("Currently processed dataset: ", layers_names[x]))
  
  output <- list()

  heatmap_out <- heatmap_top_var(data[[paste0("set",x,"_anno")]], global_settings = global_settings, 
                                 working_directory = working_directory,
                                 input = layer_specific_outputs[[paste0("set",x)]]$part1, x = x, plot_HM = plot_HM, 
                                 anno_groups = anno_groups)
  output[["heatmap_out"]] <- heatmap_out
  
  GFC_all_genes <- GFC_calculation(info_dataset = data[[paste0("set", x, "_anno")]], 
                                   global_set = global_settings,
                                   layer_set = layer_settings[[paste0("set", x)]],
                                   input1 = layer_specific_outputs[[paste0("set",x)]]$part1,
                                   input2 = output,
                                   grouping_v = grouping_v)
  
  output[["GFC_all_genes"]] <- GFC_all_genes
  
  return(output)
}


choose_auto_cutoff <- function(min_R2 = 0.8,
                               delta_R2 = 0.01,
                               node_stability_frac = 0.99,
                               max_components_frac = 0.05,
                               min_cutoff_allowed = 0.1,
                               min_scan_points = 5,
                               use_components = TRUE){
  cutoff_vec <- c()

  for (x in seq_along(layer_specific_outputs)) {
    cco <- layer_specific_outputs[[x]]$part1$cutoff_calc_out

    if (is.null(cco) || is.null(cco$cutoff_stats_concise)) {
      cutoff_vec <- c(cutoff_vec, NA_real_)
      next
    }

    df <- cco$cutoff_stats_concise
    if (!all(c("R.squared", "no_nodes") %in% colnames(df))) {
      cutoff_vec <- c(cutoff_vec, NA_real_)
      next
    }

    # Build a data.frame with explicit cutoff column
    df <- as.data.frame(df)
    df$cutoff <- as.numeric(rownames(df))
    df <- df[stats::complete.cases(df[, c("R.squared", "no_nodes")]), , drop = FALSE]

    if (nrow(df) == 0) {
      cutoff_vec <- c(cutoff_vec, NA_real_)
      next
    }

    # If the scan is very coarse, warn by returning NA
    if (nrow(df) < min_scan_points) {
      message("Layer ", x, ": cutoff scan has fewer than ", min_scan_points, " points; no automatic cutoff selected.")
      cutoff_vec <- c(cutoff_vec, NA_real_)
      next
    }

    # 1. Target R^2: cutoffs with R.squared >= min_R2
    valid_r2_df <- df[df$R.squared >= min_R2, , drop = FALSE]

    if (nrow(valid_r2_df) == 0) {
      max_r2 <- max(df$R.squared, na.rm = TRUE)
      valid_r2_df <- df[df$R.squared >= (max_r2 - delta_R2), , drop = FALSE]
    }

    # 2. Network stability safeguard
    max_nodes <- max(df$no_nodes, na.rm = TRUE)
    node_threshold <- max_nodes * node_stability_frac
    components_col <- if ("no_of_networks" %in% colnames(df)) "no_of_networks" else NA_character_

    if (use_components && !is.na(components_col)) {
      stable_df <- valid_r2_df[
        valid_r2_df$no_nodes >= node_threshold &
          valid_r2_df[[components_col]] < valid_r2_df$no_nodes * max_components_frac,
        ,
        drop = FALSE
      ]
    } else {
      # If we don't have a components column, use only node stability
      stable_df <- valid_r2_df[valid_r2_df$no_nodes >= node_threshold, , drop = FALSE]
    }

    if (nrow(stable_df) == 0) {
      # Fallback: minimum cutoff in valid_r2_df, or cutoff with max R² overall
      if (nrow(valid_r2_df) > 0) {
        suggested_cutoff <- min(valid_r2_df$cutoff, na.rm = TRUE)
      } else {
        suggested_cutoff <- df$cutoff[which.max(df$R.squared)]
      }
      suggested_cutoff <- max(suggested_cutoff, min_cutoff_allowed)
      cutoff_vec <- c(cutoff_vec, suggested_cutoff)
      next
    }

    # 3. Choose minimum cutoff that meets R² and stability criteria
    optimal_cutoff <- min(stable_df$cutoff, na.rm = TRUE)
    optimal_cutoff <- max(optimal_cutoff, min_cutoff_allowed)
    cutoff_vec <- c(cutoff_vec, optimal_cutoff)
  }

  cutoff_vec
}


plot_deg_dist <- function(cut_off_vec){
  for(x in 1:length(layers)){
    if(!cut_off_vec[x] %in% layer_specific_outputs[[paste0("set",x)]]$part1$cutoff_stats$cutoff){
      tmp_vec <- layer_specific_outputs[[paste0("set",x)]]$part1$cutoff_stats$cutoff
      cut_off_vec[x] <- tmp_vec[which.min(abs(tmp_vec - cut_off_vec[x]))]
      cutoff_vec[x] <<- cut_off_vec[x]
    }
    stats_calculated_optimal_cutoff <- layer_specific_outputs[[paste0("set",x)]]$part1$cutoff_stats[layer_specific_outputs[[paste0("set",x)]]$part1$cutoff_stats$cutoff == cut_off_vec[x], c("degree", "Probs")]
    
    stats <- layer_specific_outputs[[paste0("set", x)]][["part1"]][["cutoff_calc_out"]][["cutoff_stats_concise"]]
    stats$cut_off <- rownames(stats)
    stats <- stats[stats$cut_off == cutoff_vec[x],]
    
    
    dd_plot_calculated_optimal = ggplot(stats_calculated_optimal_cutoff,aes(x=log(degree), y= log(Probs))) +
      geom_point() +
      geom_smooth(method="lm") +
      theme_bw() + 
      ggtitle(paste0("Cut-off: ",cut_off_vec[x], "; R: ", round(stats[1],3), "; no. edges: ",
                     stats[2], "; no. nodes: ", stats[3], "; no. networks: ", stats[4]))+
      theme(plot.title = element_text(size = 10))
    
    print(dd_plot_calculated_optimal)
    ggsave(paste0("Degree_distribution_plot_",layers_names[x], "_", cut_off_vec[x], ".pdf"),
           dd_plot_calculated_optimal, device = cairo_pdf, width = 10, height = 8, path = paste0(working_directory,global_settings$save_folder))
    
  }
  
}
