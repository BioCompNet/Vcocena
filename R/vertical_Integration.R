create_network_per_layer <- function() {
  for (x in seq_along(layer_specific_outputs)) {
    edgelist <- layer_specific_outputs[[x]]$part2$heatmap_out$filt_cutoff_data
    edgelist$weight <- edgelist$rval
    edgelist$rval <- NULL
    edgelist$pval <- NULL
    layer_specific_outputs[[x]][["network"]] <<-
      igraph::graph_from_data_frame(edgelist, directed = FALSE)
  }
}


cluster_networks <- function(cluster_algo = "auto",
                             no_of_iterations = 10,
                             max_cluster_count_per_gene = 1,
                             min_cluster_size = global_settings$min_nodes_number_for_cluster,
                             mute = FALSE,
                             min_nodes_for_clustering = 3,
                             min_edges_for_clustering = 1) {
  summaries <- vector("list", length(layer_specific_outputs))
  names(summaries) <- paste0("set", seq_along(layer_specific_outputs))
 
  for (x in seq_along(layer_specific_outputs)) {
    network <- layer_specific_outputs[[x]]$network
    layer_specific_outputs[[x]][["cluster_calc"]] <<- list()

    # Skip clustering if network is too small
    if (igraph::gorder(network) < min_nodes_for_clustering ||
        igraph::gsize(network) < min_edges_for_clustering) {
      if (!mute) {
        message("Layer ", x, " (", layers_names[x],
                "): network too small for clustering (nodes = ",
                igraph::gorder(network), ", edges = ", igraph::gsize(network), ").")
      }
      empty_df <- data.frame(
        clusters         = character(0),
        gene_no          = integer(0),
        gene_n           = character(0),
        cluster_included = character(0),
        color            = character(0),
        conditions       = character(0),
        grp_means        = character(0),
        vertexsize       = integer(0),
        stringsAsFactors = FALSE
      )
      layer_specific_outputs[[x]]$cluster_calc[["cluster_information"]] <<- empty_df
      summaries[[x]] <- list(
        layer_name         = layers_names[x],
        n_clusters         = 0L,
        n_included_clusters = 0L,
        n_white_genes      = igraph::gorder(network)
      )
      next
    }

    info <- cluster_calculation(
      igraph                   = network,
      cluster_algo             = cluster_algo,
      no_of_iterations         = no_of_iterations,
      max_cluster_count_per_gene = max_cluster_count_per_gene,
      min_cluster_size         = min_cluster_size,
      GFC                      = layer_specific_outputs[[x]]$part2$GFC_all_genes,
      mute                     = mute
    )

    layer_specific_outputs[[x]]$cluster_calc[["cluster_information"]] <<- info

    # Basic diagnostics
    total_clusters   <- nrow(info)
    included_clusters <- sum(info$cluster_included == "yes")
    white_genes <- if ("color" %in% colnames(info)) {
      # Count genes in white clusters
      white_rows <- info$color == "white"
      if (any(white_rows)) {
        sum(info$gene_no[white_rows])
      } else {
        0L
      }
    } else {
      NA_integer_
    }

    summaries[[x]] <- list(
      layer_name          = layers_names[x],
      n_clusters          = total_clusters,
      n_included_clusters = included_clusters,
      n_white_genes       = white_genes
    )

    if (!mute) {
      message("Layer ", x, " (", layers_names[x], "): ",
              included_clusters, " included clusters out of ",
              total_clusters, "; white genes (small/ambiguous clusters): ",
              white_genes, ".")
    }
  }

  invisible(summaries)
}







cluster_GFC_df_per_layer <- function(){
  
  for (x in 1:length(layer_specific_outputs)){
    
    cluster_df_tmp <- layer_specific_outputs[[x]]$cluster_calc$cluster_information %>% 
      dplyr::filter(., !color == "white")
    
    conditions <- cluster_df_tmp$conditions[1] %>%
      base::strsplit(., split = "#")%>%
      unlist()
    
    vals <- apply(cluster_df_tmp, 1, function(y){
              y["grp_means"] %>% base::strsplit(., split = ",")%>%
                unlist()%>%
                as.numeric()
            })%>%
            t()%>%
            as.data.frame()
    
    df <- cbind(cluster_df_tmp$color, vals)
    
    colnames(df) <- c("cluster", conditions)
    
    layer_specific_outputs[[x]]$cluster_calc[["GFC_per_cluster"]] <<- df
    
  }
}


cluster_GFC_df_integrated <- function(){
  cluster_df_tmp <- integrated_output$cluster_calc$cluster_information %>% 
    dplyr::filter(., !color == "white")
  
  conditions <- cluster_df_tmp$conditions[1] %>%
    base::strsplit(., split = "#")%>%
    unlist()
  
  vals <- apply(cluster_df_tmp, 1, function(y){
    y["grp_means"] %>% base::strsplit(., split = ",")%>%
      unlist()%>%
      as.numeric()
  })%>%
    t()%>%
    as.data.frame()
  
  df <- cbind(cluster_df_tmp$color, vals)
  
  colnames(df) <- c("cluster", conditions)
  
  integrated_output$cluster_calc[["GFC_per_cluster"]] <<- df
}




create_pheno_vectors <- function(GFC_boundary = "auto"){
  
  if(GFC_boundary == "auto"){
    GFC_bvec <- NULL
    # determine GFC-cut-offs by using the mean positive and mean negative GFC for each layer
    for( x in 1:length(layers)){
      # e1 <- c(layer_specific_outputs[[paste0("set",x)]]$part2$GFC_all_genes$covid_No, layer_specific_outputs[[paste0("set",x)]]$part2$GFC_all_genes$covid_Yes)
      e1 <- unlist(layer_specific_outputs[[paste0("set",x)]]$part2$GFC_all_genes[, 1:(ncol(layer_specific_outputs[[paste0("set",x)]]$part2$GFC_all_genes)-1)]) %>%
        as.numeric()
      e1p <- e1[e1>=0]
      e1m <- e1[e1<0]
      GFC_bvec <- c(GFC_bvec, paste0(mean(as.numeric(e1p)), "_", mean(as.numeric(e1m))))
    }
  }else if(length(GFC_boundary) == 1){
    GFC_bvec <- rep(paste0(GFC_boundary, "_-", GFC_boundary), length(layers))
  }else{
    GFC_bvec <- rep(("1.5_-1.5"), length(layers))
  }
  
  
  for (x in 1:length(layer_specific_outputs)){
    vec <- list()
    gfc <- GFC_bvec[x]%>%
      base::strsplit(., split = "_")%>%
      unlist()
    print(gfc)
    apply(layer_specific_outputs[[x]]$cluster_calc$GFC_per_cluster, 1, function(y){
      tmp <- NULL
      for(i in 1:(length(y)-1)){
        
        if(as.numeric(y[i+1]) <= as.numeric(gfc[2])){
          tmp <- c(tmp, -1)
        }else if(as.numeric(y[i+1]) >= as.numeric(gfc[1])){
          tmp <- c(tmp, 1)
        }else{
          tmp <- c(tmp, 0) 
        }
      }
      vec[[y[1]]] <<- tmp
    })
    
    layer_specific_outputs[[x]][["cluster_vec"]] <<- vec
  }
}









find_identicals <- function(){
  output <- list()
  for(i in 1:(length(layers)-1)){
    GFC_list_1 <- layer_specific_outputs[[i]]$cluster_vec
    for(j in (i+1):length(layers)){
      tmp <- NULL
      GFC_list_2 <- layer_specific_outputs[[j]]$cluster_vec
      for(k in 1:length(GFC_list_1)){
        for(l in 1:length(GFC_list_2)){
          if(all(GFC_list_1[[k]] == GFC_list_2[[l]]) & all(GFC_list_1[[k]] == 0) == F){
            tmp <- rbind(tmp, data.frame(c1 = names(GFC_list_1)[k], c2 = names(GFC_list_2)[l]))
          }
        }
      }
    }
    output[[paste0(i, "_", j)]] <- tmp
  }
  return(output)
}

calc_vector_sum <- function(vec1, vec2){
  m <- cbind(vec1, vec2)
  sums <- apply(m, 1, sum)
  return(all(sums == 0))
}










find_zero_sums <- function(){
  output <- list()
  for(i in 1:(length(layers)-1)){
    GFC_list_1 <- layer_specific_outputs[[i]]$cluster_vec
    for(j in (i+1):length(layers)){
      tmp <- NULL
      GFC_list_2 <- layer_specific_outputs[[j]]$cluster_vec
      for(k in 1:length(GFC_list_1)){
        for(l in 1:length(GFC_list_2)){
          if(calc_vector_sum(GFC_list_1[[k]], GFC_list_2[[l]]) == T & all(GFC_list_1[[k]] == 0) == F){
            tmp <- rbind(tmp, data.frame(c1 = names(GFC_list_1)[k], c2 = names(GFC_list_2)[l]))
          }
        }
      }
    }
    output[[paste0(i, "_", j)]] <- tmp
  }
  return(output)
}










cross_corr_pair<- function(set1, set2, c1, c2){
  c_df_1 <- layer_specific_outputs[[paste0("set", set1)]]$cluster_calc$cluster_information %>%
    dplyr::filter(., color == c1) %>%
    dplyr::pull(., "gene_n") %>%
    base::strsplit(., split = ",") %>%
    unlist()
  c_df_2 <- layer_specific_outputs[[paste0("set", set2)]]$cluster_calc$cluster_information %>%
    dplyr::filter(., color == c2) %>%
    dplyr::pull(., "gene_n") %>%
    base::strsplit(., split = ",") %>%
    unlist()
  
  count_df_1 <- data[[paste0("set", set1, "_counts")]]
  count_df_1 <- count_df_1[rownames(count_df_1) %in% c_df_1, ]%>%
    t()
  count_df_2 <- data[[paste0("set", set2, "_counts")]]
  count_df_2 <- count_df_2[rownames(count_df_2) %in% c_df_2, ]%>%
    t()
  
  
  correlations <- lineup::corbetw2mat(count_df_1, count_df_2, what = "all")%>%
    reshape2::melt(.)
  return(correlations)
  
}







calculate_cross_corrs_identity <- function(threshold = 0.7, n = 100){
  
  if (length(integrated_output$identicals) == 0) {
    return(list())
  }
  
  cross_corrs <- NULL
  
  for(i in seq_len(length(integrated_output$identicals))){
    nm <- names(integrated_output$identicals)[i]
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) {
      next
    }
    sets <- base::strsplit(nm, split = "_")%>%
      unlist()
    if (length(sets) < 2) {
      next
    }
    if (is.null(integrated_output$identicals[[i]]) || nrow(integrated_output$identicals[[i]]) == 0) {
      next
    }
    cross_corrs <- c(cross_corrs, apply(integrated_output$identicals[[i]], 1, function(x){
                                            cross_corr_pair(set1 = sets[1], set2 = sets[2],
                                                            c1 = x[1], c2 = x[2])
                                      }))
  }
  
  
  #new:
  cross_corrs <- dplyr::bind_rows(cross_corrs)
  
  if (is.null(cross_corrs) || nrow(cross_corrs) == 0 || !("value" %in% colnames(cross_corrs))) {
    return(list())
  }
  max_val <- max(cross_corrs$value, na.rm = TRUE)
  if (!is.finite(max_val)) {
    return(list())
  }
  
  cutoff_sequence <- seq(from = threshold, to = max_val, length.out = n) %>%
    round(., digits = 3)%>%
    unique()
  
  output <- lapply(cutoff_sequence, function(x){
    tmp <- dplyr::filter(cross_corrs, value >= x)
    
    colnames(tmp) <- c("V1", "V2", "weight")
    return(tmp)
  })
  
  names(output) <- as.character(cutoff_sequence)
  # cross_corrs <- dplyr::bind_rows(cross_corrs)%>%
  #   dplyr::filter(., value >= threshold)
  # 
  # colnames(cross_corrs) <- c("V1", "V2", "weight")
  
  # return(cross_corrs)
  return(output)
}



cross_corrs <- function(threshold){
  if (length(integrated_output$identicals) == 0) {
    return(data.frame(V1 = character(), V2 = character(), weight = numeric()))
  }
  
  cross_corrs <- NULL
  
  for(i in seq_len(length(integrated_output$identicals))){
    nm <- names(integrated_output$identicals)[i]
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) {
      next
    }
    sets <- base::strsplit(nm, split = "_")%>%
      unlist()
    if (length(sets) < 2) {
      next
    }
    if (is.null(integrated_output$identicals[[i]]) || nrow(integrated_output$identicals[[i]]) == 0) {
      next
    }
    cross_corrs <- c(cross_corrs, apply(integrated_output$identicals[[i]], 1, function(x){
      cross_corr_pair(set1 = sets[1], set2 = sets[2],
                      c1 = x[1], c2 = x[2])
    }))
  }
  
  cross_corrs <- dplyr::bind_rows(cross_corrs)
  if (is.null(cross_corrs) || nrow(cross_corrs) == 0 || !("value" %in% colnames(cross_corrs))) {
    return(data.frame(V1 = character(), V2 = character(), weight = numeric()))
  }
  cross_corrs <- dplyr::filter(cross_corrs, value >= threshold)

  colnames(cross_corrs) <- c("V1", "V2", "weight")
  
  return(cross_corrs)
}



calculate_cross_corrs_zero_sums <- function(threshold =  -0.7){
  cross_corrs <- NULL 
  if (length(integrated_output$zero_sums) == 0) {
    return(data.frame(V1 = character(), V2 = character(), weight = numeric()))
  }
  for( i in seq_len(length(integrated_output$zero_sums))){
    nm <- names(integrated_output$zero_sums)[i]
    if (is.null(nm) || is.na(nm) || !nzchar(nm)) {
      next
    }
    sets <- base::strsplit(nm, split = "_")%>%
      unlist()
    if (length(sets) < 2) {
      next
    }
    if (is.null(integrated_output$zero_sums[[i]]) || nrow(integrated_output$zero_sums[[i]]) == 0) {
      next
    }
    cross_corrs <- c(cross_corrs, apply(integrated_output$zero_sums[[i]], 1, function(x){
      cross_corr_pair(set1 = sets[1], set2 = sets[2],
                      c1 = x[1], c2 = x[2])
    }))
  }
  cross_corrs <- dplyr::bind_rows(cross_corrs)
  if (is.null(cross_corrs) || nrow(cross_corrs) == 0 || !("value" %in% colnames(cross_corrs))) {
    return(data.frame(V1 = character(), V2 = character(), weight = numeric()))
  }
  cross_corrs <- dplyr::filter(cross_corrs, value <= threshold)
  
  colnames(cross_corrs) <- c("V1", "V2", "weight")
  cross_corrs$weight <- abs(cross_corrs$weight)
  
  return(cross_corrs)
  
  
}








get_pheno_edgelist <- function(GFC_boundary){
  edgelist <- NULL
  for(x in 1:length(layers)){
    apply(layer_specific_outputs[[x]]$cluster_calc$GFC_per_cluster, 1, function(l){
      entities <- dplyr::filter(layer_specific_outputs[[x]]$cluster_calc$cluster_information, color == l[1])%>%
        dplyr::pull(., "gene_n")%>%
        base::strsplit(., split = ",")%>%
        unlist(.)
      
      
      GFCs <- layer_specific_outputs[[x]]$part2$GFC_all_genes[layer_specific_outputs[[x]]$part2$GFC_all_genes$Gene %in% entities, ]
      
      
      for (i in 2:length(l)){
        if(as.numeric(l[i]) >= GFC_boundary){
          
          filt_entities <- layer_specific_outputs[[x]]$part2$GFC_all_genes[layer_specific_outputs[[x]]$part2$GFC_all_genes[colnames(layer_specific_outputs[[x]]$cluster_calc$GFC_per_cluster)[i]] >= GFC_boundary,]%>%
            dplyr::pull(., "Gene")
          
          el_tmp <- data.frame(V1 = filt_entities)
          el_tmp$V2 <- rep(l[1], nrow(el_tmp))
          #el_tmp$V2 <- rep(colnames(layer_specific_outputs[[x]]$cluster_calc$GFC_per_cluster)[i], nrow(el_tmp))
          el_tmp <- rbind(el_tmp, data.frame(V1 = l[1], V2 = colnames(layer_specific_outputs[[x]]$cluster_calc$GFC_per_cluster)[i]))
          
          el_tmp$weight <- 1
          el_tmp <- rbind(el_tmp, data.frame(V1 = l[1], V2 = colnames(layer_specific_outputs[[x]]$cluster_calc$GFC_per_cluster)[i], weight = 50))
          edgelist <<- rbind(edgelist, el_tmp)
          
        }else if(as.numeric(l[i])<= -GFC_boundary){
          
          filt_entities <- layer_specific_outputs[[x]]$part2$GFC_all_genes[layer_specific_outputs[[x]]$part2$GFC_all_genes[colnames(layer_specific_outputs[[x]]$cluster_calc$GFC_per_cluster)[i]] <= -GFC_boundary,]%>%
            dplyr::pull(., "Gene")
          
          el_tmp <- data.frame(V1 = filt_entities)
          el_tmp$V2 <- rep(l[1], nrow(el_tmp))
          
          #el_tmp$V2 <- rep(colnames(layer_specific_outputs[[x]]$cluster_calc$GFC_per_cluster)[i], nrow(el_tmp))
          el_tmp$weight <- 1
          el_tmp <- rbind(el_tmp, data.frame(V1 = l[1], V2 = colnames(layer_specific_outputs[[x]]$cluster_calc$GFC_per_cluster)[i], weight = 50))
          edgelist <<- rbind(edgelist, el_tmp)
        }
      }
      
    })
  }
  return(edgelist)
}










get_layer_edgelists <- function(){
  edgelist <- NULL
  rescale_min <- min(edge_attr(layer_specific_outputs[[1]]$network)[[1]])
  rescale_max <- max(edge_attr(layer_specific_outputs[[1]]$network)[[1]])
  for(i in 1:length(layers)){
    tmp <- igraph::as_edgelist(layer_specific_outputs[[i]]$network) %>%
      as.data.frame()
    colnames(tmp) <- c("V1", "V2")
    tmp$weight <- rescale(edge_attr(layer_specific_outputs[[i]]$network)[[1]], to = c(rescale_min,rescale_max))
    edgelist <- rbind(edgelist, tmp)
  }
  return(edgelist)
}







create_integrated_edgelist <- function(el){
  
  # create nodes for condition and their edges to entities above set GFC cutoff:
  #pheno_edgelist <- get_pheno_edgelist(GFC_boundary = GFC_boundary)
  integrated_edgelist <- get_layer_edgelists()
  #integrated_edgelist <- rbind(pheno_edgelist, layer_edgelist)
  cross_edgelist <- el
  
  if (is.null(cross_edgelist) || nrow(cross_edgelist) == 0) {
    return(integrated_edgelist)
  }
  
  rescale_min <- min(edge_attr(layer_specific_outputs[[1]]$network)[[1]])
  rescale_max <- max(edge_attr(layer_specific_outputs[[1]]$network)[[1]])
  
  cross_edgelist$weight <- rescale(cross_edgelist$weight, to = c(rescale_min,rescale_max))
  integrated_edgelist <- rbind(integrated_edgelist, cross_edgelist)
  return(integrated_edgelist)
}




merge_GFCs <- function(){
  merged_GFCs <- NULL
  for(x in 1:length(layers)){
    merged_GFCs <- rbind(merged_GFCs, layer_specific_outputs[[x]]$part2$GFC_all_genes)
  }
  return(merged_GFCs)
}




cluster_integrated_network <- function(network, mute = F){
  
    intergrated_GFCs <- merge_GFCs()
    info <- invisible(cluster_calculation(igraph = network,
                                cluster_algo = "cluster_louvain", #cluster_louvain
                                no_of_iterations = 1,
                                max_cluster_count_per_gene = 1,
                                min_cluster_size = global_settings$min_nodes_number_for_cluster,
                                GFC = intergrated_GFCs,
                                mute = mute))
    return(info)
}


# cluster_integrated_network2 <- function(network){
#   
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


cluster_anti_network <- function(){
  
  network <- integrated_output$anti_net
  integrated_output[["cluster_calc_anti"]] <<- list()
  intergrated_GFCs <- merge_GFCs()%>%dplyr::filter(., Gene %in% integrated_output$cross_anti_correlations$V1 | Gene %in% integrated_output$cross_anti_correlations$V2)
  
  integrated_output$cluster_calc[["cluster_information_anti"]] <<- invisible(cluster_calculation(igraph = network,
                                                                                  cluster_algo = "auto",
                                                                                  no_of_iterations = 50,
                                                                                  max_cluster_count_per_gene = 1,
                                                                                  min_cluster_size = global_settings$min_nodes_number_for_cluster,
                                                                                  GFC = intergrated_GFCs))
}



layer_distribution_in_clusters <- function(network_type){
  entity_to_layer <- NULL
  for(x in 1:length(layers)){
    tmp <- data[[paste0("set", x, "_counts")]]%>%
      row.names()
    entity_to_layer <- rbind(entity_to_layer, data.frame(entity = tmp, layer = layers_names[x]))
  }
  stats_df <- list()
  if(network_type == "co"){
    cluster_df <- dplyr::filter(integrated_output$cluster_calc$cluster_information, !color == "white")
  }else if(network_type == "anti"){
    cluster_df <- dplyr::filter(integrated_output$cluster_calc$cluster_information_anti, !color == "white")
  }else{
    print("invalid network type")
    stop()
  }
  
  
  stats_df <- apply(cluster_df, 1, function(y){
    ents <- base::strsplit(y[3], ",")%>%
      unlist()
    cluster_entity_to_layer <- entity_to_layer[entity_to_layer$entity %in% ents,]
    df <- table(cluster_entity_to_layer$layer)%>%as.data.frame()
    g <- ggplot(data=df, aes(x=Var1, y=Freq)) +
      geom_bar(stat="identity")+ggtitle(y[5])
    return(g)
  })
  return(stats_df)
}





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
    info <- cluster_integrated_network(network = mergednet, mute = T)
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


cumulative_dist <- function(P, x){
  return(P(x))
}

mixture_score <- function(P, freq_table){
  score <- 0
  for(i in 1:length(freq_table)){
    # total entities in cluster:
    total <- sum(freq_table[[i]]$Freq)
    # print(total)
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


find_all_cluster_corrs <- function(){  
  output <- list()
  for(i in 1:(length(layers)-1)){
    for(line_i in 1:nrow(layer_specific_outputs[[paste0("set", i)]]$cluster_calc$cluster_information)){
      if(!layer_specific_outputs[[paste0("set", i)]]$cluster_calc$cluster_information$color[line_i] == "white"){
        c1 <- layer_specific_outputs[[paste0("set", i)]]$cluster_calc$cluster_information$color[line_i]
        v1 <- layer_specific_outputs[[paste0("set", i)]]$cluster_calc$cluster_information$grp_means[line_i] %>%
          base::strsplit(., split = ",")%>%
          unlist()%>%
          as.numeric()
      }
      for(j in (i+1):length(layers)){
        for(line_j in 1:nrow(layer_specific_outputs[[paste0("set", j)]]$cluster_calc$cluster_information)){
          if(!layer_specific_outputs[[paste0("set", j)]]$cluster_calc$cluster_information$color[line_j] == "white"){
            c2 <- layer_specific_outputs[[paste0("set", j)]]$cluster_calc$cluster_information$color[line_j]
            v2 <- layer_specific_outputs[[paste0("set", j)]]$cluster_calc$cluster_information$grp_means[line_j] %>%
              base::strsplit(., split = ",")%>%
              unlist()%>%
              as.numeric()
            output[[paste0(i, "_", j)]][[paste0(c1, "_", c2)]] <- cor(v1, v2, method = "pearson")
          }
        }
      }
    }
    tmp <- bind_rows(output[[paste0(i,"_", j)]]) %>%
      t()
    tmp <- data.frame(clusters = rownames(tmp), correlation = tmp[,1])
    output[[paste0(i,"_", j)]] <- tmp
  }
  return(output)
}



ids <- function(crosscorrlist = integrated_output[["all_cluster_corrs"]], min_corr = 0.97){
  output <- list()
  for(i in names(crosscorrlist)){
    tmp <- crosscorrlist[[i]]%>%
      dplyr::filter(., correlation > min_corr)
    tmp <- lapply(tmp$clusters, function(x){
      base::strsplit(x, split = "_")%>%
        unlist()
    })%>%rlist::list.rbind()%>%
      as.data.frame()
    colnames(tmp) <- c("c1", "c2")
    output[[i]] <- tmp
  }
  return(output)
}

ops <- function(crosscorrlist = integrated_output[["all_cluster_corrs"]], min_corr = 0.97){
  output <- list()
  for(i in names(crosscorrlist)){
    tmp <- crosscorrlist[[i]]%>%
      dplyr::filter(., correlation < -min_corr)
    tmp <- lapply(tmp$clusters, function(x){
      base::strsplit(x, split = "_")%>%
        unlist()
    })%>%rlist::list.rbind()%>%
      as.data.frame()
    colnames(tmp) <- c("c1", "c2")
    output[[i]] <- tmp
  }
  return(output)
}


find_entity <- function(entity =  NULL){
  ctg <- do.call(rbind, apply(integrated_output$cluster_calc$cluster_information, 1, function(x){
    if(!x[5] == "white"){
      g <- x[3] %>%
        base::strsplit(., split = ",")%>%
        unlist()
      tmp <- data.frame(entity = g, cluster = x[5])
      return(tmp)
    }
  }))
  return(ctg[ctg$entity == entity,] %>% dplyr::pull(., "cluster"))
}
