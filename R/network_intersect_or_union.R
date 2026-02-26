get_intersection <- function(with){
  
  combined_edgelist <- layer_specific_outputs[[with]]$part2$heatmap_out$filt_cutoff_data
  combined_edgelist$merged <- paste0(combined_edgelist$V1 %>% as.character(), 
                                     combined_edgelist$V2 %>% as.character())
  for (x in 1:length(layer_specific_outputs)){
    if(!paste0("set", x) == with){
      tmp <- layer_specific_outputs[[x]]$part2$heatmap_out$filt_cutoff_data
      tmp$merged <- paste0(tmp$V1 %>% as.character(), 
                           tmp$V2 %>% as.character())
      tmp <- dplyr::filter(tmp, merged %in% combined_edgelist$merged)
      combined_edgelist <- rbind(combined_edgelist, tmp)
    }
  }
  combined_edgelist <- unique(combined_edgelist)%>%
    dplyr::arrange(., desc(rval))
  combined_edgelist <- combined_edgelist[!duplicated(combined_edgelist$merged),]
  combined_edgelist$merged <- NULL
  return(combined_edgelist)
}



get_union <- function(){
  
  combined_edgelist <- NULL
  
  for (x in 1:length(layer_specific_outputs)){
    combined_edgelist <- rbind(combined_edgelist, layer_specific_outputs[[x]]$part2$heatmap_out$filt_cutoff_data)
  }
  combined_edgelist$weight <- combined_edgelist$rval
  combined_edgelist$rval <- NULL
  combined_edgelist$pval <- NULL
  
  return(combined_edgelist)
}

build_integrated_network <- function(multi_edges = "min"){
  merged_net <- igraph::graph_from_data_frame(integrated_output$combined_edgelist, directed=FALSE)
  merged_net <- igraph::simplify(merged_net, edge.attr.comb=list(weight=multi_edges, "ignore"))
  return(merged_net)
}