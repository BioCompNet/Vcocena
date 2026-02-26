merge_GFCs <- function(GFC_when_not_expressed){
  col_names_new_GFC <- NULL
  for (x in 1:length(layer_specific_outputs)){
    new_col <- colnames(layer_specific_outputs[[x]]$part2$GFC_all_genes)[-ncol(layer_specific_outputs[[x]]$part2$GFC_all_genes)]

    col_names_new_GFC <- c(col_names_new_GFC, new_col)
  }
  new_GFC <- NULL
  for(y in igraph::get.vertex.attribute(integrated_output$merged_net)$name){
    
    line <- NULL
    
    
    for(z in 1:length(layer_specific_outputs)){
      if(y %in% layer_specific_outputs[[z]]$part2$GFC_all_genes$Gene){

        GFC_tmp <- layer_specific_outputs[[z]]$part2$GFC_all_genes
        
        line <- c(line, GFC_tmp[GFC_tmp$Gene == y, colnames(GFC_tmp)[1:(ncol(GFC_tmp)-1)]]) %>%
          unlist(.)
        
      }else{
        line <- c(line, rep(GFC_when_not_expressed, (length(colnames(layer_specific_outputs[[z]]$part2$GFC_all_genes))-1)))
      }
    }
    
    line <- as.data.frame(line) %>% t()
    new_GFC <- rbind(new_GFC, line)
  }
  new_GFC <- as.data.frame(new_GFC)
  
  colnames(new_GFC) <- c(col_names_new_GFC)
  new_GFC$Gene <- igraph::get.vertex.attribute(integrated_output$merged_net)$name
  rownames(new_GFC) <- new_GFC$Gene
  return(new_GFC)
}
