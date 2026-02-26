plot_integrated_network <- function(network, layout = NULL, gene_labels = NULL, save = T, label_offset = 50){
  gene_to_cluster <- do.call(rbind, apply(integrated_output$cluster_calc$cluster_information,1, function(x){
    tmp <- x["gene_n"] %>%
      base::strsplit(., split = ",")%>%
      unlist(.)
    data.frame(name = tmp, color = rep(x["color"], length(tmp)))
  }))
  gene_to_cluster$name <- as.character(gene_to_cluster$name)
  gene_to_cluster$color <- as.character(gene_to_cluster$color)
  gene_to_cluster <- gene_to_cluster[match(V(network)$name, gene_to_cluster$name),]
  V(network)$color <- gene_to_cluster$color
  network <- delete.vertices(network, gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name"))
  #merged_net <<- network
  if(is.null(layout)){
    set.seed(123)
    # Default to igraph layout if layout_algorithm is not set or not "cytoscape"
    if(!is.null(global_settings$layout_algorithm) &&
       global_settings$layout_algorithm == "cytoscape"){
      RCy3::createNetworkFromIgraph(network, "my igraph 2", )
      RCy3::setNodeColorMapping(table.column = "color", mapping.type = "passthrough")
      RCy3::layoutNetwork("force-directed")
      l <- RCy3::getNodePosition()
      l <-l[match(get.vertex.attribute(network)$name, rownames(l)),]
      l[] <- lapply(l, as.numeric)
      l <- as.matrix(l)
    }else{
      l <- igraph::layout.fruchterman.reingold(network) 
      rownames(l) <- V(network)$name
    }
    
  }else{
    l <- layout
  }
  
  if(!is.null(gene_labels)){
    
    
    new_genes_df <- data.frame(name = paste0("label_", c(1:length(gene_labels))), label = gene_labels %>% as.character())
    mean_coord_x <- median(l[,1])
    mean_coord_y <- median(l[,2])
    new_indeces <- match(new_genes_df$label, get.vertex.attribute(network)$name)
    new_genes_df$coords_x <- l[new_indeces,1]
    new_genes_df$coords_y <- l[new_indeces,2]
    left_up <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y >= mean_coord_y)%>%
      dplyr::pull(., label)
    left_down <- dplyr::filter(new_genes_df, coords_x <= mean_coord_x & coords_y < mean_coord_y)%>%
      dplyr::pull(., label)
    right_up <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y >= mean_coord_y)%>%
      dplyr::pull(., label)
    right_down <- dplyr::filter(new_genes_df, coords_x > mean_coord_x & coords_y < mean_coord_y)%>%
      dplyr::pull(., label)
    

    
    new_position_l <- matrix(cbind(rep(ceiling(min(l[,1]))-label_offset, length(left_up)+ length(left_down)), 
                                 seq(from=ceiling(max(l[,2])), to = ceiling(min(l[,2])), 
                                     length.out = length(left_up)+ length(left_down))), ncol = 2)
    new_position_r <- matrix(cbind(rep(ceiling(max(l[,1]))+label_offset, length(right_up)+ length(right_down)), 
                                   seq(from=ceiling(max(l[,2])), to = ceiling(min(l[,2])), 
                                       length.out = length(right_up)+ length(right_down))), ncol = 2)
    new_position <- rbind(new_position_l, new_position_r)
    
    colnames(new_position) <- colnames(l)
    l2 <- rbind(l, new_position)%>%as.matrix()
    new_genes_df_l <- new_genes_df[new_genes_df$label %in% c(left_up, left_down),]
    new_genes_df_l <- new_genes_df_l[order(new_genes_df_l$coords_y, decreasing = T),]
    new_genes_df_r <- new_genes_df[new_genes_df$label %in% c(right_up, right_down),]
    new_genes_df_r <- new_genes_df_r[order(new_genes_df_r$coords_y, decreasing = T),]
    #new_genes_df <- new_genes_df[match(c(left_up, left_down, right_up, right_down), new_genes_df$label),]
    new_genes_df <- rbind(new_genes_df_l, new_genes_df_r)
    
    
    
    network2 <- igraph::add.vertices(network, nv = length(new_genes_df$name), 
                            attr = list(name = new_genes_df$name))
    new_labels <- lapply(get.vertex.attribute(network2)$name, function(x){
      if(x %in% get.vertex.attribute(network)$name){
        NA
      }else{
        dplyr::filter(new_genes_df, name == x)%>%
          dplyr::pull(., "label")
      }
    })%>% unlist()
    
    V(network2)$label <- new_labels
    
    new_edges <- matrix(c(match(new_genes_df$label, get.vertex.attribute(network2)$name), 
                        match(new_genes_df$name, get.vertex.attribute(network2)$name)), 
                        ncol = 2, byrow = F)
    new_edges <- as.vector(t(new_edges))
      
    network2 <- igraph::add.edges(network2, new_edges)
    
    new_edge_color <- apply(igraph::get.edgelist(network2),1, function(x){
      if(x[2] %in% new_genes_df$name){
        igraph::get.vertex.attribute(network2, name = "color", index = x[1])
      }else{
        "lightgrey"
      }
    })
    
    new_label_color <- lapply(get.vertex.attribute(network2)$name, function(x){
      if(x %in% new_genes_df$name){
        tmp_gene_n <- dplyr::filter(new_genes_df, name == x)%>%
          dplyr::pull(., label)
        igraph::get.vertex.attribute(network2, name = "color", index = tmp_gene_n)
      }else{
        NA
      }
    })%>% unlist()
    
    
    new_frame_color <- lapply(get.vertex.attribute(network2)$name, function(x){
      if(x %in% new_genes_df$name){
        "white"
      }else{
        "black"
      }
    })%>% unlist()
    
    if(save == T){
      Cairo::CairoPDF(file = paste0(working_directory, global_settings$save_folder, "/network_col_by_module.pdf"), 
                      width = 15, height = 15)
      
      igraph::plot.igraph(network2, vertex.size = 3, vertex.label = new_labels, vertex.label.cex = 0.75,
                          layout = l2, main = "co-expression network coloured by module",
                          edge.color = new_edge_color, vertex.label.color = new_label_color,
                          vertex.frame.color = new_frame_color)
      
      dev.off()
    }
    
    
    igraph::plot.igraph(network2, vertex.size = 3, vertex.label = new_labels, vertex.label.cex = 0.75,
                        layout = l2, main = "co-expression network coloured by module",
                        edge.color = new_edge_color, vertex.label.color = new_label_color,
                        vertex.frame.color = new_frame_color)
    
    integrated_output$cluster_calc[["labelled_network"]] <<- network2
    integrated_output$cluster_calc[["network_col_by_module"]] <<- network
  }else{
    
    if(save == T){
      Cairo::CairoPDF(file = paste0(working_directory, global_settings$save_folder, "/network_col_by_module.pdf"), 
                      width = 15, height = 15)
      
      igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3, 
                          layout = l, main = "co-expression network coloured by cluster")
      
      dev.off()
    }
    
    
    igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3, 
                        layout = l, main = "co-expression network coloured by cluster")
    
    integrated_output$cluster_calc[["labelled_network"]] <<- list()
    integrated_output$cluster_calc[["network_col_by_module"]] <<- network
  }
  
  
  if(is.null(layout)){
    return(l)
  }
  
}


plot_GFC_network <- function(GFCs = integrated_output$GFC_all_layers, network = integrated_output$cluster_calc$network_col_by_module){
  GFCs_wog <- GFCs[,colnames(GFCs)[!colnames(GFCs) == "Gene"]]
  colors <- apply(GFCs_wog,2, function(x){
    y <- ((round(x, digits = 1) +2) *10)+1
    colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, by = .1)))[y] 
  }) %>% as.data.frame() 
  colors$name <- GFCs$Gene
  colors[] <- lapply(colors, as.character)
  
  gene_to_cluster <- do.call(rbind, apply(integrated_output$cluster_calc$cluster_information,1, function(x){
    tmp <- x["gene_n"] %>%
      base::strsplit(., split = ",")%>%
      unlist(.)
    data.frame(name = tmp, color = rep(x["color"], length(tmp)))
  }))
  
  gene_to_cluster$name <- as.character(gene_to_cluster$name)
  gene_to_cluster$color <- as.character(gene_to_cluster$color)
  if(length(V(network)$name[!V(network)$name %in% gene_to_cluster$name]) > 0){
    gene_to_cluster <- rbind(gene_to_cluster, data.frame(name = V(network)$name[!V(network)$name %in% gene_to_cluster$name],
                                                         color = "white"))
    gene_to_cluster$name <- as.character(gene_to_cluster$name)
    gene_to_cluster$color <- as.character(gene_to_cluster$color)
  }
  colors <- colors[match(V(network)$name, colors$name),]
  
  if(any((gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name"))%in% 
         igraph::get.vertex.attribute(network)$name)){
    remove_v <- gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name")
    remove_v <- remove_v[remove_v %in% V(network)$name]
    network <- delete.vertices(network, remove_v)
  }
  
  colors <- colors[!colors$name %in% (gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name")), ]
  #l <- integrated_output$cluster_calc$layout[match(V(network)$name, rownames(integrated_output$cluster_calc$layout)), ]
  for(x in colnames(colors)){
    if(! x == "name"){
      V(network)$color <- colors %>% dplyr::pull(., var = x)
      Cairo::CairoPDF(file = paste0(working_directory, global_settings$save_folder, "/GFC_network_", x, ".pdf"), 
                      width = 15, height = 15)
      igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3,
                          layout = integrated_output$layout,
                          main = x)
      dev.off()
      igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3,
                          layout = integrated_output$layout, 
                          main = x)
    }
    
  }
}




plot_layer_network <- function(network){
  
  gene_to_layer <- NULL
  for(i in 1:length(layers)){
    gene_to_layer <- rbind(gene_to_layer, data.frame(entity = V(layer_specific_outputs[[i]]$network)$name, layer = layers_names[i]))
  }
  
  
  gene_to_layer <- gene_to_layer[gene_to_layer$entity %in% V(network)$name, ]
  
  gene_to_layer <- gene_to_layer[match(V(network)$name, gene_to_layer$entity),]
  
  
  
  
  colors <- dplyr::pull(gene_to_layer, "layer") %>% as.factor() 
  V(network)$color <- colors 
  Cairo::CairoPDF(file = paste0(working_directory, global_settings$save_folder, "/layer_network.pdf"), 
                  width = 15, height = 15)
  igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3,
                      layout = integrated_output$layout,
                      main = "network coloured by layer")
  dev.off()
  igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3,
                      layout = integrated_output$layout, 
                      main = "network coloured by layer")
    
}

# plot_GFC_network <- function(GFCs, network){
#   GFCs_wog <- GFCs[,colnames(GFCs)[!colnames(GFCs) == "Gene"]]
#   colors <- apply(GFCs_wog,2, function(x){
#     y <- ((round(x, digits = 1) +2) *10)+1
#     colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")))(length(seq(-2, 2, by = .1)))[y] 
#   }) %>% as.data.frame() 
#   colors$name <- GFCs$Gene
#   colors[] <- lapply(colors, as.character)
#   
#   gene_to_cluster <- do.call(rbind, apply(integrated_output$cluster_calc$cluster_information,1, function(x){
#     tmp <- x["gene_n"] %>%
#       base::strsplit(., split = ",")%>%
#       unlist(.)
#     data.frame(name = tmp, color = rep(x["color"], length(tmp)))
#   }))
#   
#   gene_to_cluster$name <- as.character(gene_to_cluster$name)
#   gene_to_cluster$color <- as.character(gene_to_cluster$color)
#   colors <- colors[match(V(network)$name, colors$name),]
#   if(any((gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name"))%in% 
#          igraph::get.vertex.attribute(integrated_output$merged_net)$name)){
#     remove_v <- gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name")
#     remove_v <- remove_v[remove_v %in% igraph::get.vertex.attribute(integrated_output$merged_net)$name]
#     network <- delete.vertices(network, remove_v)
#   }
#   
#   colors <- colors[!colors$name %in% (gene_to_cluster[gene_to_cluster$color == "white",] %>% dplyr::pull(., "name")), ]
#   for(x in colnames(colors)){
#     if(! x == "name"){
#       V(network)$color <- colors %>% dplyr::pull(., var = x)
#       igraph::plot.igraph(network, vertex.label = NA, vertex.size = 3,
#                           layout = integrated_output$cluster_calc$layout, main = x)
#     }
#     
#   }
# }
