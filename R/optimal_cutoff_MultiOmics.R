optcut_fun <- function(cutoff_stats){
  
  
  print(paste0("...calculating optimal cutoff..."))
  
  output <- list()
  
  cutoff_stats_concise = cutoff_stats %>% 
    dplyr::select(R.squared, cutoff, no_edges, no_nodes, no_of_networks) %>% 
    dplyr::distinct()
  
  cutoff_stats_concise = cutoff_stats_concise %>% filter(no_of_networks!=0)
  rownames(cutoff_stats_concise) = cutoff_stats_concise$cutoff
  cutoff_stats_concise = cutoff_stats_concise[,-2]
  
  
  
  crit_minmax = c("max","max", "max", "min" )
  names(crit_minmax) = colnames(cutoff_stats_concise)
  
  
  
  normalizationTypes <- rep("percentageOfMax", ncol(cutoff_stats_concise))
  names(normalizationTypes) = colnames(cutoff_stats_concise)
  if(nrow(cutoff_stats_concise)==0){
    message("cutoff_stats_concise is empty; falling back to all cutoffs (including no_of_networks == 0).")
    cutoff_stats_concise <- cutoff_stats %>% 
      dplyr::select(R.squared, cutoff, no_edges, no_nodes, no_of_networks) %>% 
      dplyr::distinct()
    rownames(cutoff_stats_concise) <- cutoff_stats_concise$cutoff
    cutoff_stats_concise <- cutoff_stats_concise[, -2]
  }
  nPT = normalizePerformanceTable(cutoff_stats_concise[,c("R.squared", "no_edges", "no_nodes", "no_of_networks")], normalizationTypes)
  w = c(0.5,0.1,0.5, -1)
  names(w) <- colnames(nPT)
  ws<-weightedSum(nPT,w)
  ranked_ws <- rank(-ws) %>% sort()
  
  
  calculated_optimal_cutoff <- as.numeric(names(ranked_ws[1]))
  
  
  
  stats_calculated_optimal_cutoff <- cutoff_stats[cutoff_stats$cutoff == calculated_optimal_cutoff, c("degree", "Probs")]
  
  dd_plot_calculated_optimal = ggplot(stats_calculated_optimal_cutoff,aes(x=log(degree), y= log(Probs))) +
    geom_point() +
    geom_smooth(method="lm") +
    theme_bw() + 
    ggtitle(paste0("Calculated optimal correlation cut-off [",calculated_optimal_cutoff, "]"))
  
 
  output[["cutoff_stats_concise"]] <- cutoff_stats_concise
  output[["dd_plot_calculated_optimal"]] <- dd_plot_calculated_optimal
  output[["optimal_cutoff"]] <- calculated_optimal_cutoff
  return(output)

}

