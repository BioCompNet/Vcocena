merge_clusters <- function(to_merge, cluster_info){
  for (entry in to_merge){
    entry <- base::strsplit(entry, split = ",") %>%
      BiocGenerics::unlist(.)
    
    tmp <- cluster_info[cluster_info$color %in% entry,]
    cluster_info <- cluster_info[!cluster_info$color %in% entry,]
    grp_means <- tmp$grp_means %>%
      base::strsplit(., split = ",")%>%
      BiocGenerics::unlist(.) %>%
      as.numeric(.) %>%
      matrix(., ncol = (length(.)/length(tmp$grp_means)), byrow = T)%>%
      apply(., 2, mean) %>%
      as.character()%>%
      BiocGenerics::paste(., collapse = ",")
    

    merged_clusters <- data.frame(clusters = tmp$clusters[1], 
                                  gene_no = sum(tmp$gene_no),
                                  gene_n = BiocGenerics::paste(tmp$gene_n, collapse = ","),
                                  cluster_included = "yes",
                                  color = tmp$color[1],
                                  conditions = tmp$conditions[1],
                                  grp_means = grp_means,
                                  vertexsize = 3)
    cluster_info <- rbind(cluster_info, merged_clusters)
    
  }
  return(cluster_info)
}



combine_graphs <- function(layers){
  output <- list()
  output[["combined_graphs"]]<- igraph::graph_from_literal()
  
  if(!is.null(layers$gene_exp)){
    output$combined_graphs <- igraph::disjoint_union(output$combined_graphs,
                                                           layer_specific_outputs$gene_exp$heatmap_out$filt_cutoff_graph)
  }
  if(!is.null(layers$lnc)){
    output$combined_graphs <- igraph::union(output$combined_graphs, 
                                                  layer_specific_outputs$lnc$heatmap_out$filt_cutoff_graph)
  }
  
  return(output)
  
}


algo_alluvial <- function(network){
  cluster_algo_list =c("cluster_fast_greedy",
                       "cluster_infomap",
                       "cluster_walktrap",
                       "cluster_leading_eigen",
                       "cluster_label_prop")
    
  
  color.cluster <- c("orchid", "maroon", "darkgreen",  "darkorange", "darkgrey", "gold", "steelblue", "indianred",
                     "pink", "lightgreen", "lightblue","sandybrown",   "khaki",  "turquoise","darkblue",
                     "cadetblue","greenyellow","cyan", "thistle", "darkmagenta", "coral", "red", "blue",
                     "green", "yellow", "brown", "black", "darkgoldenrod", "cornsilk", "firebrick", "deeppink",
                     "dodgerblue", "lightpink", "midnightblue", "slategray", "aquamarine",
                     "chocolate", "darkred", "navy", "olivedrab", "peachpuff",
                     "seagreen", "plum", "tomato", "snow", "wheat")
  
  # louvain as base:
  # set.seed(168575)
  # cfg = cluster_louvain(network)
  # louvain <- data.frame(gene = get.vertex.attribute(network)$name, cluster = cfg$membership)
  # louvain$cluster <- lapply(louvain$cluster, function(x){
  #   color.cluster[x]
  # })%>% unlist()
  # tmp2 <- table(louvain$cluster) %>% as.data.frame()
  # louvain <- louvain[louvain$cluster %in% (dplyr::filter(tmp2, Freq >= global_settings$min_nodes_number_for_cluster)%>% 
  #                                dplyr::pull(., "Var1")),]
  #not_in <- get.vertex.attribute(network)$name[!get.vertex.attribute(network)$name %in% louvain$gene]
  #louvain <- rbind(louvain, data.frame(gene = not_in, cluster = "white" ))
  louvain <- NULL
  for(c in unique(integrated_output$cluster_calc$cluster_information$color)){
    tmp <- dplyr::filter(integrated_output$cluster_calc$cluster_information, color == c)%>%
      dplyr::pull(., "gene_n")%>%
      base::strsplit(., split = ",")%>%
      unlist(.)
    louvain <- rbind(louvain, data.frame(gene = tmp, cluster = rep(c, length(tmp))))
  }
  
  louvain$cluster <- as.factor(louvain$cluster)
  
  output <- list()
  output[["cluster_louvain"]] <- louvain
  for(a in cluster_algo_list){
    # set.seed(168575)
    # new:
    set.seed(1)
    set.seed(.Random.seed[1])



    cfg = get(a)(network)
    tmp <- data.frame(gene = get.vertex.attribute(network)$name, cluster = cfg$membership)
    tmp$cluster <- lapply(tmp$cluster, function(x){
      color.cluster[x]
    })%>% unlist()
    tmp2 <- table(tmp$cluster) %>% as.data.frame()
    tmp <- tmp[tmp$cluster %in% (dplyr::filter(tmp2, Freq >= global_settings$min_nodes_number_for_cluster)%>% 
                 dplyr::pull(., "Var1")),]
    tmp$cluster <- as.factor(tmp$cluster)
    
    
    test3 <- merge(louvain, tmp , by = "gene", all = T)
    colnames(test3) <- c("gene", "cluster1", "cluster2")
    test3$cluster1 <- as.character(test3$cluster1)
    test3$cluster2 <- as.character(test3$cluster2)
    test3[is.na(test3)] <- "white"
    test3$Freq <- 1
    test3$merged <- paste0(test3$cluster1, test3$cluster2)
    
    
    
    links <- test3[,2:ncol(test3)] %>% unique()
    links$Freq <- lapply(links$merged, function(x){
      length(test3$merged[test3$merged == x])
    })%>% unlist()
    links$merged <- NULL
    links$cluster1_name <- paste0(links$cluster1, "_louvain")
    links$cluster2_name <- paste0(links$cluster2, base::strsplit(a, split = "cluster")[[1]][2])
    
    nodes <- data.frame(name=c(as.character(links$cluster1_name)%>% unique(), as.character(links$cluster2_name)%>% unique()))
    
    links$IDsource <- match(links$cluster1_name, nodes$name)-1 
    links$IDtarget <- match(links$cluster2_name, nodes$name)-1
    colnames(links) <- c("source_col", "target_col", "value", "source", "target","IDsource", "IDtarget")
    node_col <- data.frame(name=c(as.character(links$source_col)%>% unique(), as.character(links$target_col)%>% unique()))
    #my_color <- 'd3.scaleOrdinal() .domain(["group_A", "group_B","group_C", "group_D", "group_E", "group_F", "group_G", "group_H"]) .range(["blue", "blue" , "blue", "red", "red", "yellow", "purple", "purple"])'
    my_color <- paste0('d3.scaleOrdinal() .domain([',paste0(paste0('"',as.character(nodes$name), collapse = '", '),'"', collapse = '"'),
                       ']) .range([', paste0(paste0('"', as.character(node_col$name), collapse = '", '), '"',collapse = '"'), '])')
    p <- sankeyNetwork(Links = links, Nodes = nodes,
                       Source = "IDsource", Target = "IDtarget",
                       Value = "value", NodeID = "name", 
                       sinksRight=FALSE, colourScale=my_color,
                       nodeWidth=40, fontSize=13, nodePadding=10)
    #%>% 
     # htmlwidgets::prependContent(htmltools::tags$h1(paste0("cluster_louvain vs. ", a)))
    print(p)
    output[[a]] <- tmp
  }
  return(output)
}



cluster_calculation <- function(igraph,
                                cluster_algo,
                                no_of_iterations,
                                max_cluster_count_per_gene,
                                min_cluster_size,
                                GFC, 
                                mute = F,
                                leiden_resolution = 1.0) {
  
  
  g = igraph
  # If the graph has no vertices, return an empty cluster-information
  # data.frame with the expected columns.
  if (igraph::gorder(g) == 0) {
    if (!mute) {
      message("Graph has no vertices; returning empty cluster information.")
    }
    return(data.frame(
      clusters = character(0),
      gene_no = integer(0),
      gene_n = character(0),
      cluster_included = character(0),
      color = character(0),
      conditions = character(0),
      grp_means = character(0),
      vertexsize = integer(0),
      stringsAsFactors = FALSE
    ))
  }

  comps = count_components(g)
  
  # Base set of igraph cluster algorithms considered for "auto" mode
  cluster_algo_list <- c(
    "cluster_label_prop",
    "cluster_fast_greedy",
    "cluster_louvain",
    "cluster_infomap",
    "cluster_walktrap",
    "cluster_leading_eigen"
  )
  # Optionally include Leiden if available in igraph
  if (exists("cluster_leiden", where = asNamespace("igraph"), inherits = FALSE)) {
    cluster_algo_list <- c(cluster_algo_list, "cluster_leiden")
  }
  
  algos_to_use = switch(cluster_algo=="auto", cluster_algo_list, cluster_algo)
  
  

  
  cluster_calculations =function(graph_obj, algo, case, iter) {

    set.seed(168575)
    # Special handling for Leiden if available
    if (algo == "cluster_leiden" &&
        exists("cluster_leiden", where = asNamespace("igraph"), inherits = FALSE)) {
      cfg <- igraph::cluster_leiden(
        graph_obj,
        resolution_parameter = leiden_resolution
      )
    } else {
      cfg <- get(algo)(graph_obj)
    }
    
    mod_score =modularity(graph_obj, cfg$membership)
    
    mod_df= data.frame(modularity_score=mod_score, cluster_algorithm=algo, stringsAsFactors = F)
    
    ##making switch so that in the end when only the best algorithm is to be used then the same function can be used
    output = switch(case, best= cfg$membership, test= mod_df, final=cfg)
    if(!mute){
      print(paste0(algo," algorithm tested"))
    }
    
    return(output)
  }
  
  if(cluster_algo=="auto"){
    if(!mute){
      print(algos_to_use)
    }
    
    df_modularity_score = do.call("rbind", lapply(algos_to_use,
                                                  cluster_calculations,
                                                  graph_obj=g,
                                                  case="test", iter=1))
    
    
    cluster_algo_used <- df_modularity_score %>%
      dplyr::filter(modularity_score == max(modularity_score)) %>%
      dplyr::pull(cluster_algorithm)
    if (length(cluster_algo_used) > 1) {
      # Resolve ties deterministically using the original algorithm order
      cluster_algo_used <- cluster_algo_list[cluster_algo_list %in% cluster_algo_used][1]
    }
    if(!mute){
      print(paste(cluster_algo_used, " will be used based on your input (if not auto option was specified) or the highest modularity score "))
    }
    
    
  }else{
    cluster_algo_used <- cluster_algo
    if(!mute){
      print(paste(cluster_algo_used, " will be used."))
    }
    
    
  }
  
  
  igraph_list = list()
  igraph_list[[1]] = g
  
  ###apply the best clustering algorithm
  gene_which_cluster=do.call("cbind", lapply(1:no_of_iterations,
                                             cluster_calculations,
                                             algo=cluster_algo_used,
                                             case="best",
                                             graph_obj=g))
  
  
  
  ##frequency and identity of cluster assingment of genes
  if(base::ncol(gene_which_cluster) > 1) {
    gene_cluster_ident = apply(gene_which_cluster,1, function(x){
      if(length(unique(x)) > max_cluster_count_per_gene) {   #LISA: was >=
        0
      }else{
        names(which(table(x) == max(table(x))))[1]
      }
    })
  } else{ gene_cluster_ident = gene_which_cluster[,1]}
  
  
  
  white_genes_clustercounts <- as.integer(grep(gene_cluster_ident, pattern = "\\b0\\b") %>%  
                                            length() %>% as.character())
  
  if(!mute){
    print(paste(white_genes_clustercounts, "genes were assigned to more than", max_cluster_count_per_gene,
                "clusters. These genes are assigned to Cluster 0 and will be painted white in the network."))
  }
  
  
  
  
  cluster_Data = data.frame(genes=vertex_attr(g, "name"),
                            clusters= paste0("Cluster ",gene_cluster_ident),
                            stringsAsFactors = FALSE)
  
  #summarize the data
  #produces a table where col are cluster name, number of components,
  #names of genes in cluster
  
  dfk=cluster_Data %>%
    dplyr::count(clusters,genes) %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(gene_no= sum(n), gene_n = paste0(genes,collapse = ","), dplyr.summarise.inform = F) %>%
    dplyr::mutate(cluster_included=ifelse(gene_no>=min_cluster_size, "yes", "no"), color="white")
  
  
  
  dfk <- dfk[order(dfk$gene_no, decreasing = T),]
  
  ##ggplot

  # Dynamic colour palette for clusters: generate enough distinct colours
  base_cols <- c(
    "orchid", "maroon", "darkgreen", "darkorange", "darkgrey", "gold", "steelblue", "indianred",
    "pink", "lightgreen", "lightblue", "sandybrown", "khaki", "turquoise", "darkblue",
    "cadetblue", "greenyellow", "cyan", "thistle", "darkmagenta", "coral", "red", "blue",
    "green", "yellow", "brown", "black", "darkgoldenrod", "cornsilk", "firebrick", "deeppink",
    "dodgerblue", "lightpink", "midnightblue", "slategray", "aquamarine",
    "chocolate", "darkred", "navy", "olivedrab", "peachpuff",
    "seagreen", "plum", "tomato", "snow", "wheat"
  )
  included_clusters <- dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", "clusters", drop = TRUE]
  n_clusters <- length(unique(included_clusters))
  if (n_clusters == 0) {
    color.cluster <- base_cols[1]
  } else {
    color.cluster <- grDevices::colorRampPalette(base_cols)(n_clusters)
  }

  plot_clusters <- ggplot(data = dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", ],
                          aes(x = clusters)) +
    geom_bar(aes(fill = clusters)) +
    scale_fill_manual(values = color.cluster)
  
  if (nrow(dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", ]) > 0) {
    plot_clust <- ggplot_build(plot_clusters)
    dfk[dfk$cluster_included == "yes" & dfk$clusters != "Cluster 0", "color" ] <- plot_clust$data[[1]]["fill"]
  }
  
  
  
  
  white_genes_clustersize <- as.integer(dfk %>% dplyr::filter(cluster_included=="no")%>%
                                          dplyr::summarise(n=sum(gene_no), dplyr.summarise.inform = F) %>% purrr::map(1))
  
  if(!mute){
    print(paste0(white_genes_clustersize, " genes were assigned to clusters with a smaller size than the defined minimal cluster size of ",
                 min_cluster_size, " genes per cluster. These genes will also be painted white in the network."))
  }
  
  
  
  
  
  # for each cluster produces a row of means per condition (info data) for all genes within the cluster
  #included_clusters = subset(dfk, included=="yes")
  cluster_df=dfk
  gfc_dat = GFC
  
  gfc_mean_clustergene=function(rownum, cluster_df, gfc_dat){
    d1 = cluster_df[rownum,]
    gene_names= d1["gene_n"] %>%
      stri_split_regex(pattern = ",") %>%
      unlist()
    gfc_means = gfc_dat[gfc_dat$Gene%in%gene_names,] %>%
      dplyr::select(-Gene) %>%
      colMeans()
    d1$conditions = paste0(names(gfc_means), collapse = "#")
    d1$grp_means = paste0(round(gfc_means,3) , collapse = ",")
    return(d1)
  }
  
  
  dfk_allinfo=do.call("rbind", lapply(1:nrow(dfk), gfc_mean_clustergene, cluster_df= dfk,
                                      gfc_dat=GFC))
  dfk_allinfo$vertexsize = ifelse(dfk_allinfo$cluster_included=="yes",3,1)
  
  
  return(dfk_allinfo)
  
}

update_clustering_algorithm <- function(alluvials = NULL, new_algo = NULL, gtc = NULL){
  if(is.null(gtc)){
    gtc <- alluvials[[new_algo]]
  }else{
    gtc <- gtc
    colnames(gtc) <- c("gene", "cluster")
  }
  
  gtc[] <- lapply(gtc, as.character)
  new_cluster_info <- NULL
  for(c in unique(gtc$cluster)){
    tmp <- gtc[gtc$cluster == c, ]
    gene_no <- nrow(tmp)
    gene_n <- paste0(tmp$gene, collapse = ",")
    if(c == "white"){
      cluster_included <- "no"
      vertexsize <- 1
    }else{
      cluster_included <- "yes"
      vertexsize <- 3
    }
    color <- c
    conditions <- paste0(colnames(integrated_output$GFC_all_layers)[1:(ncol(integrated_output$GFC_all_layers)-1)], collapse = "#")
    gfc_means = integrated_output$GFC_all_layers[integrated_output$GFC_all_layers$Gene %in% tmp$gene,] %>%
      dplyr::select(-Gene) %>%
      colMeans()
    grp_means = paste0(round(gfc_means,3) , collapse = ",")
    
    new_cluster_info <- rbind(new_cluster_info,
                              data.frame(clusters = c,
                                         gene_no = gene_no,
                                         gene_n = gene_n,
                                         cluster_included = cluster_included,
                                         color = c,
                                         conditions = conditions,
                                         grp_means = grp_means,
                                         vertexsize = vertexsize, 
                                         stringsAsFactors = F))
    
  }
  return(new_cluster_info)
}
