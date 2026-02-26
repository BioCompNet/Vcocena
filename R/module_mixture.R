
cumulative_dist <- function(P, x){
  return(P(x))
}

mixture_score <- function(P, freq_table){
  score <- 0
  for(i in 1:length(freq_table)){
    # total entities in cluster:
    total <- sum(freq_table[[i]]$Freq)
    
    for(l in layers_names){
      if(l %in% freq_table[[i]]$Var1){
        ratio = (dplyr::filter(freq_table[[i]], Var1 == l)%>%
                   dplyr::pull(., Freq)%>%
                   as.numeric())/total
        if(ratio > 1/length(layers)){
          score <- score + cumulative_dist(P = P, x = 1 - ratio)
        }else{
          score <- score + cumulative_dist(P = P, x = ratio)
        }
      }else{
        score <- score + 0
      }
    }
  }
  return(score)
}

# create_integrated_edgelist2 <- function(el){
#   
#   # create nodes for condition and their edges to entities above set GFC cutoff:
#   #pheno_edgelist <- get_pheno_edgelist(GFC_boundary = GFC_boundary)
#   integrated_edgelist <- get_layer_edgelists()
#   #integrated_edgelist <- rbind(pheno_edgelist, layer_edgelist)
#   cross_edgelist <- el
#   cross_edgelist$weight <- rescale(cross_edgelist$weight, to = c(0.5,1), from = c(min(cross_edgelist$weight), max(cross_edgelist$weight)))
#   integrated_edgelist <- rbind(integrated_edgelist, cross_edgelist)
#   return(integrated_edgelist)
# }

# 
# cluster_integrated_network2 <- function(network){
#   
#   cluster_calc<- list()
#   intergrated_GFCs <- merge_GFCs()
#   info <- cluster_calculation(igraph = network,
#                               cluster_algo = "cluster_louvain",
#                               no_of_iterations = 1,
#                               max_cluster_count_per_gene = 1,
#                               min_cluster_size = global_settings$min_nodes_number_for_cluster,
#                               GFC = intergrated_GFCs, 
#                               mute = T)
#   return(info)
# }

find_cross_corr_cutoff <- function(n, threshold){
  edgelists <- calculate_cross_corrs_identity(threshold = threshold, n = n)
  
  entity_to_layer <- NULL
  for(x in 1:length(layers)){
    tmp <- data[[paste0("set", x, "_counts")]]%>%
      row.names()
    entity_to_layer <- rbind(entity_to_layer, data.frame(entity = tmp, layer = layers_names[x]))
  }
  
  num_edges <- c()
  mix_score <- list()
  num_clusters <- c()
  
  P = ecdf(rnorm(n = 1000, mean = 1/length(layers), sd = 1/(2*length(layers)))) 
  
  for(i in names(edgelists)){
    num_edges <- c(num_edges, nrow(edgelists[[i]]))
    el <- create_integrated_edgelist(el = edgelists[[i]])
    mergednet <- igraph::graph_from_data_frame(el, directed = F)%>%
      igraph::simplify(., remove.multiple = TRUE, remove.loops = TRUE)
    info <- cluster_integrated_network(network = mergednet)
    cluster_df <- dplyr::filter(info, !color == "white")
    stats_list <- apply(cluster_df, 1, function(y){
      ents <- base::strsplit(y[3], ",")%>%
        unlist()
      cluster_entity_to_layer <- entity_to_layer[entity_to_layer$entity %in% ents,]
      df <- table(cluster_entity_to_layer$layer)%>%as.data.frame()
    })
    mix_score[[i]] <- mixture_score(P = P, freq_table = stats_list)
    num_clusters <- c(num_clusters, nrow(cluster_df))
  }
  num_edges <- as.numeric(num_edges)
  mix_score <- list.rbind(mix_score) %>% as.data.frame()
  colnames(mix_score) <- c("mix_score")
  mix_score$cutoff <- rownames(mix_score)
  mix_score$num_clusters <- num_clusters
  mix_score$num_cross_edges <- num_edges
  max_score <- max(mix_score$mix_score)
  print(paste0("The prodposed cross-correlation cutoff based on the highest mixture-score is: ", pull(mix_score[mix_score$mix_score == max_score,], "cutoff"), "."))
  return(mix_score)
}
# test <- find_cross_corr_cutoff(n = 100)

# p1 <- ggplot(test)+
#   geom_point(aes(x = cutoff, y = mix_score), col = "lightgreen", pch = 15, size = 2)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90))
# # p1
# p2 <- ggplot(test)+
#   geom_point(aes(x = cutoff, y = num_clusters), col = "magenta", pch = 15, size = 2)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90))
# # p2
# p3 <- ggplot(test)+
#   geom_point(aes(x = cutoff, y = num_cross_edges), col = "lightblue", pch = 15, size = 2)+
#   theme_bw()+
#   theme(axis.text.x = element_text(angle = 90))
# # p3
# cp <- cowplot::plot_grid(p1, p2, p3, align = "v", ncol = 1)
# cp
