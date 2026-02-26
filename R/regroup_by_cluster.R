GeneToCluster <- function(){
  gtc <- do.call(rbind, apply(integrated_output$cluster_calc$cluster_information,1, function(x){
    tmp <- x["gene_n"] %>%
      base::strsplit(., split = ",")%>%
      unlist(.)
    data.frame(gene = tmp, color = rep(x["color"], length(tmp)))
  }))
  return(gtc)
}


intra_sample_FC <- function(l, gtc){
  
  if(is.null(gtc)){
    gene_to_cluster <- GeneToCluster()
  }else{
    colnames(gtc) <- c("gene", "cluster")
    gene_to_cluster <- gtc
    gene_to_cluster$gene <- as.character(gene_to_cluster$gene)
    gene_to_cluster$cluster <- as.character(gene_to_cluster$cluster)
    colnames(gene_to_cluster) <- c("gene", "color")
  }
  
  
  counts <- data[[paste0("set", l , "_counts")]]
  
  mean_expression_per_sample <- apply(counts, 2, mean)
  
  mean_expression_per_cluster <- NULL
  
  for(c in unique(gene_to_cluster$color)){
    if(!c == "white"){
      genes <- gene_to_cluster[gene_to_cluster$color == c, ] %>%
        dplyr::pull(., "gene")
      
      filt_counts <- counts[rownames(counts) %in% genes, ]
      tmp <- apply(filt_counts, 2, mean)%>%
        as.data.frame()%>%
        t()%>%
        as.data.frame()
      rownames(tmp) <- c
      colnames(tmp) <- colnames(counts)
      mean_expression_per_cluster <- rbind(mean_expression_per_cluster, tmp)
    }
  }
  mean_expression_per_cluster <- rbind(mean_expression_per_cluster, mean_expression_per_sample)
  
  FC_from_mean_per_cluster <- apply(mean_expression_per_cluster,2, function(x){
    gtools::foldchange(x, x[length(x)])
  }) %>% as.data.frame()
  colnames(FC_from_mean_per_cluster) <- colnames(mean_expression_per_cluster)
  rownames(FC_from_mean_per_cluster) <- rownames(mean_expression_per_cluster)
  FC_from_mean_per_cluster <- FC_from_mean_per_cluster[1:(nrow(FC_from_mean_per_cluster)-1),]
  # rownames(mean_expression_per_cluster) <- c(rownames(mean_expression_per_cluster)[1:(nrow(mean_expression_per_cluster)-1)],
  #                                            "mean")
  return(FC_from_mean_per_cluster)
}   


plot_PCA_cluster <- function(gtc = NULL, algo = NULL){
  plotlist <- list()
  pca_list <- list()
  for(l in 1:length(layers)){
    tmp_out <- list()
    if(is.null(gtc)){
      test <- intra_sample_FC(l)
    }else{
      test <- intra_sample_FC(l, gtc = gtc)
      
    }
    
    pca <- prcomp(t(test), scale = T)
    tmp_out[[1]] <- pca
    pca.var <- pca$sdev^2
    pca.var.per <- data.frame(pc = 1:length(pca.var), val = round(pca.var/sum(pca.var)*100,1))
    tmp_out[[2]] <- find_top85(pca.var.per)
    anno <- data[[paste0("set", l, "_anno")]]
    pca.data <- data.frame(Sample = rownames(pca$x),
                           X = pca$x[,1],
                           Y = pca$x[,2],
                           Group = dplyr::pull(anno , global_settings$voi))
    p <- ggplot(pca.data, aes(x = X, y = Y,  label= Sample, color = Group))+
      geom_point(size = 4)+
      #geom_text(col = "darkblue")+
      scale_color_jama()+
      ggtitle(paste0(layers_names[l], " by module - ", algo))+
      ylab(paste0("PC 2", " (", pca.var.per[2, 2], "%)"))+
      xlab(paste0("PC 1", " (", pca.var.per[1, 2], "%)"))
    plotlist[[l]] <- p
    pca_list[[l]] <- pca
  }
 cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 2, align = "hv")
 plot(cp)
 return(pca_list)
}


k_means_prep_regrouping_from_cluster <- function(){
  output <- list()
  for(l in 1:length(layers)){
    tmp_out <- list()
    test <- intra_sample_FC(l)
    pca <- prcomp(t(test), scale = T)
    tmp_out[[1]] <- pca
    pca.var <- pca$sdev^2
    pca.var.per <- data.frame(pc = 1:length(pca.var), val = round(pca.var/sum(pca.var)*100,1))
    tmp_out[[2]] <- find_top85(pca.var.per)
    pca.data <- data.frame(Sample = rownames(pca$x),
                           X = pca$x[,1],
                           Y = pca$x[,2])
    p <- ggplot(pca.data, aes(x = X, y = Y,  label= Sample))+
      geom_point(size = 4)+
      #geom_text(col = "darkblue")+
      scale_color_jama()+
      ggtitle(layers_names[l])+
      ylab("PC 2")+
      xlab("PC 1")
    plot(p)
    
    ellbow <- ellbowplot(pca)
    silhouette <- silhouette_plot(pca)
    output[[paste0("set_", l)]] <- tmp_out
  }
  return(output)
}


k_means_regrouping_from_cluster <- function(kvec, prep_out){
  for(l in 1:length(layers)){
    set.seed(123)
    i <- prep_out[[l]][[2]]
    pca <- prep_out[[l]][[1]]
    pca.data <- data.frame(Sample = rownames(pca$x),
                           X = pca$x[,1],
                           Y = pca$x[,2])
    fit <- kmeans(x = pca$x[,1:i], centers = kvec[l], nstart = 25, iter.max = 1000) 
    # get cluster means
    aggregate(pca$x[,1:i],by=list(fit$cluster),FUN=mean)
    
    dat <- data.frame(Sample = rownames(pca$x[,1:i]), cluster =  paste0(layers_names[l], "_", fit$cluster))
    dat$cluster <- as.factor(dat$cluster)
    colnames(dat) <- c("Sample", global_settings$voi)
    
    anno <- data[[paste0("set", l, "_anno")]]
    anno[paste0(global_settings$voi, "_old")] <- anno[global_settings$voi]
    anno[global_settings$voi] <- NULL
    anno$Sample <- rownames(anno)
    
    anno <- merge(anno, dat, by = "Sample")
    rownames(anno) <- anno$Row.names
    anno$Row.names <- NULL
    #colnames(anno) <- colnames(data[[paste0("set", l, "_anno")]])
    
    
    
    pca.data <- data.frame(Sample = rownames(pca$x),
                           X = pca$x[,1],
                           Y = pca$x[,2],
                           Group = dplyr::pull(anno, global_settings$voi))
    #print(head(pca.data))
    p <- ggplot(pca.data, aes(x = X, y = Y,  label= Sample, color = Group))+
      geom_point(size = 4)+
      #geom_text(col = "darkblue")+
      scale_color_jama()+
      ggtitle(layers_names[l])+
      ylab("PC 2")+
      xlab("PC 1")
    plot(p)
    
    if(!global_settings$control == "none"){
      rownames(anno) <- anno$Sample
      data[[paste0("set", l, "_anno")]] <<- determine_new_controls(anno, i = l)
      
    }else{
      rownames(anno) <- anno$Sample
      data[[paste0("set", l, "_anno")]] <<- anno
      
    }
  }
}



regrouping_preparations <- function(by = NULL, algo = NULL){
  
  if(is.null(by)){
    print("The data will not be regrouped.")
  }else if(by == "topvar"){
    for(i in 1:length(integrated_output$PCAs$topvar)){
      ellbowplot(integrated_output$PCAs$topvar[[i]], title = layers_names[i])
      silhouette_plot(integrated_output$PCAs$topvar[[i]], title = layers_names[i])
    }
    
  }else{
    for(i in 1:length(integrated_output$PCAs[[paste0("module ", algo)]])){
      ellbowplot(integrated_output$PCAs[[paste0("module ", algo)]][[i]], title = layers_names[i])
      silhouette_plot(integrated_output$PCAs[[paste0("module ", algo)]][[i]], title = layers_names[i])
    }
  }
}


regroup_data_module <- function(dat, k, l){
  set.seed(123)
  fit <- kmeans(x = dat, centers = k, nstart = 25, iter.max = 1000) # 5 cluster solution
  # get cluster means
  aggregate(dat,by=list(fit$cluster),FUN=mean)
  
  dat2 <- data.frame(Sample = rownames(dat), cluster =  paste0(layers_names[l], "_", fit$cluster))
  dat2$cluster <- as.factor(dat2$cluster)
  colnames(dat2) <- c("Sample", global_settings$voi)
  
  anno <- data[[paste0("set", l, "_anno")]]
  anno[paste0(global_settings$voi, "_old")] <- anno[global_settings$voi]
  anno[global_settings$voi] <- NULL
  anno$Sample <- rownames(anno)
  
  anno <- merge(anno, dat2, by = "Sample")
  rownames(anno) <- anno$Row.names
  anno$Row.names <- NULL
  
  
  if(!global_settings$control == "none"){
    rownames(anno) <- anno$Sample
    data[[paste0("set", l, "_anno")]] <<- determine_new_controls(anno, i = l)
    
  }else{
    rownames(anno) <- anno$Sample
    data[[paste0("set", l, "_anno")]] <<- anno
    
  }
}


regroup_data <- function(by = NULL, algo = NULL, kvec){
  if(is.null(by)){
    print("The data will not be regrouped.")
  }else if(by == "topvar"){
    for(l in 1:length(layers)){
      if(!kvec[l] == 0){
        pca <- integrated_output$PCAs$topvar[[l]]
        pca.var <- pca$sdev^2
        pca.var.per <- data.frame(pc = 1:length(pca.var), val = round(pca.var/sum(pca.var)*100,1))
        i <- find_top85(pca.var.per)
        print(paste0("The first ", i, " PCs accounts for at lest 85% of the variance."))
        dat <- pca$x[,1:i]
        regroup_data_topvar(dat = dat,
                            k = kvec[l], x = l)
      }
    }
    
  }else{
    for(l in 1:length(layers)){
      if(!kvec[l] == 0){
        pca <- integrated_output$PCAs[[paste0("module ", algo)]][[l]]
        pca.var <- pca$sdev^2
        pca.var.per <- data.frame(pc = 1:length(pca.var), val = round(pca.var/sum(pca.var)*100,1))
        i <- find_top85(pca.var.per)
        print(paste0("The first ", i, " PCs accounts for at lest 85% of the variance."))
        dat <- pca$x[,1:i]
        regroup_data_module(dat = dat, k = kvec[l], l = l) 
      }
    }
  }
}
