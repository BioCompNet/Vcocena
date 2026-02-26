plot_PCAs <- function(PC_a, PC_b, i, es){
  output <- list()
  pca <- prcomp(t(layer_specific_outputs[[paste0("set",i)]][["part1"]][["topvar"]]), scale = T)
  pca.var <- pca$sdev^2
  pca.var.per <- data.frame(pc = 1:length(pca.var), val = round(pca.var/sum(pca.var)*100,1))
  if(nrow(pca.var.per)> 20){
    pca.var.per <- pca.var.per[1:20,]
  }
  b <- ggplot(pca.var.per, aes(x = pc, y = val))+
    geom_line()+
    geom_point(pch = 15)+
    ggtitle(layers_names[i])+
    ylab("Variance in %")+
    xlab("PC")
  pca.data <- data.frame(Sample = rownames(pca$x),
                         X = pca$x[,PC_a],
                         Y = pca$x[,PC_b],
                         Group = data[[paste0("set",i, "_anno")]][[global_settings$voi]])
  p <- ggplot(pca.data, aes(x = X, y = Y, col = Group, label= Sample))+
    geom_point(size = 4)+
    ylab(paste0("PC ", PC_b, " (", pca.var.per[PC_b, 2], "%)"))+
    xlab(paste0("PC ", PC_a, " (", pca.var.per[PC_a, 2], "%)"))+
    ggtitle(layers_names[i])+
    scale_color_jama()
  
  load_a <- sort(abs(pca$rotation[,PC_a]), decreasing = T)
  load_b <- sort(abs(pca$rotation[,PC_b]), decreasing = T)
  top10_a <- data.frame(gene = names(load_a)[1:10], loading_score = pca$rotation[names(load_a)[1:10],PC_a])
  top10_b <- data.frame(gene = names(load_b)[1:10], loading_score = pca$rotation[names(load_b)[1:10],PC_b])
  
  
  if(es){
    ellbow <- ellbowplot(pca)
    output[["kmeans_data"]] <- ellbow[["dat"]]
    silhouette <- silhouette_plot(pca)
    cp2 <- cowplot::plot_grid(ellbow[["p"]], silhouette, ncol = 2, nrow = 1)
  }
  
  output[["pca"]] <- pca
  return(output)
}



plot_PCA_topvar <- function(){
  plotlist <- list()
  pca_list <- list()
  for(i in 1:length(layers)){
    
    pca <- prcomp(t(layer_specific_outputs[[paste0("set",i)]][["part1"]][["topvar"]]), scale = T)
    pca.var <- pca$sdev^2
    pca.var.per <- data.frame(pc = 1:length(pca.var), val = round(pca.var/sum(pca.var)*100,1))
    
    pca.data <- data.frame(Sample = rownames(pca$x),
                           X = pca$x[,1],
                           Y = pca$x[,2],
                           Group = data[[paste0("set",i, "_anno")]][[global_settings$voi]])
    p <- ggplot(pca.data, aes(x = X, y = Y, col = Group, label= Sample))+
      geom_point(size = 4)+
      ylab(paste0("PC 2", " (", pca.var.per[2, 2], "%)"))+
      xlab(paste0("PC 1", " (", pca.var.per[1, 2], "%)"))+
      ggtitle(paste0(layers_names[i], " by topvar"))+
      scale_color_jama()
    
    plotlist[[i]] <- p
    pca_list[[i]] <- pca
  }
  cp <- cowplot::plot_grid(plotlist = plotlist, ncol = 2)
  plot(cp)
  return(pca_list)
}


ellbowplot <- function(PCA, title = NULL){
  pca <- PCA
  pca.var <- pca$sdev^2
  pca.var.per <- data.frame(pc = 1:length(pca.var), val = round(pca.var/sum(pca.var)*100,1))
  i <- find_top85(pca.var.per)
  print(paste0("The first ", i, " PCs accounts for at lest 85% of the variance."))
  dat <- pca$x[,1:i]
  wss <- (nrow(dat)-1)*sum(apply(dat,2,var))
  set.seed(123)
  for (j in 2:15) wss[j] <- sum(kmeans(dat,
                                       centers=j,
                                       nstart = 25, 
                                       iter.max = 1000)$withinss)
  if(is.null(title)){
    p <- plot(2:15, wss[2:15], type="b", xlab="Number of Clusters",
              ylab="Within groups sum of squares", main = "ellbowplot")
  }else{
    p <- plot(2:15, wss[2:15], type="b", xlab="Number of Clusters",
              ylab="Within groups sum of squares", main = paste0("ellbowplot - ", title))
  }
  
  output <- list()
  output[["dat"]] <- dat
  output[["plot"]] <- p
  return(output)
}

silhouette_plot <- function(PCA, title = NULL){
  pca <- PCA
  pca.var <- pca$sdev^2
  pca.var.per <- data.frame(pc = 1:length(pca.var), val = round(pca.var/sum(pca.var)*100,1))
  i <- find_top85(pca.var.per)
  dat <- pca$x[,1:i]
  set.seed(123)
  df <- NULL
  for(k in 2:15){
    km <- kmeans(dat, centers = k, nstart=25, iter.max = 1000)
    ss <- silhouette(km$cluster, dist(dat))
    df <- rbind(df, data.frame(k = k, score = mean(ss[, 3])))
  }
  if(is.null(title)){
    p <- plot(df$k, type='b', df$score, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE,
              main = "silhouette plot")
  }else{
    p <- plot(df$k, type='b', df$score, xlab='Number of clusters', ylab='Average Silhouette Scores', frame=FALSE,
              main = paste0("silhouette plot - ", title))
  }
  
  return(p)
}



find_top85 <- function(var){
  perc <- 0
  for (i in 1:nrow(var)){
    perc <- perc + var[i,2]
    if(perc>= 85){
      return(i)
    }
  }
}


regroup_data_topvar <- function(dat, k, x){
  if(k > 0){
    set.seed(123)
    fit <- kmeans(x = dat, centers = k, nstart = 25, iter.max = 1000) # 5 cluster solution
    # get cluster means
    aggregate(dat,by=list(fit$cluster),FUN=mean)

    dat <- data.frame(sample = rownames(dat), cluster = paste0(layers_names[x], "_", fit$cluster))
    dat$cluster <- as.factor(dat$cluster)
    colnames(dat) <- c("Sample", global_settings$voi)
    anno <- data[[paste0("set", x, "_anno")]]
    
    anno[paste0(global_settings$voi, "_old")] <- anno[global_settings$voi]
    anno[global_settings$voi] <- NULL
    anno$Sample <- rownames(anno)
    
    anno <- merge(anno, dat, by = "Sample")
    rownames(anno) <- anno$Sample
    anno$Row.names <- NULL
    if(!global_settings$control == "none"){
      data[[paste0("set", x, "_anno")]] <<- determine_new_controls(anno, i = x)
    }else{
      data[[paste0("set", x, "_anno")]] <<- anno
    }
  }
  }

reset_regrouping <- function(reset = "no", sets_to_reset, algo, gtc = NULL){
  if(reset == "yes"){
    for(x in sets_to_reset){
      anno <- data[[paste0("set", x, "_anno")]]
      anno[[global_settings$voi]] <- anno[[paste0(global_settings$voi, "_old")]]
      anno[[paste0(global_settings$voi, "_old")]] <- NULL
      data[[paste0("set", x, "_anno")]] <<- anno
    }
  }else{
    for(x in 1:length(layers)){
      layer_specific_outputs[[paste0("set",x)]][["part2"]] <<- run_expression_analysis_2(x, grouping_v = NULL, 
                                                                                             plot_HM = T)
    }
    integrated_output$GFC_all_layers <<- merge_GFCs(GFC_when_not_expressed = -global_settings$range_GFC)
    if(!is.null(gtc)){
      integrated_output$cluster_calc$cluster_information <<- update_clustering_algorithm(alluvials = integrated_output$alluvials, new_algo = algo, gtc = gtc)
    }else{
      integrated_output$cluster_calc$cluster_information <<- update_clustering_algorithm(alluvials = integrated_output$alluvials, new_algo = algo)
    }
    
  }
}


determine_new_controls <- function(anno, i){
  ctrl_samples <- anno[grepl(global_settings$control, anno[[paste0(global_settings$voi, "_old")]], ignore.case = T),] %>% rownames()
  maxvec <- NULL
  anno[global_settings$voi] <- lapply(anno[global_settings$voi], as.character)
  anno$sample_name <- rownames(anno)
  for(x in unique(anno[global_settings$voi] %>% dplyr::pull(., global_settings$voi))){
    tmp <- anno[anno[global_settings$voi] == x, "sample_name"]
    l <- intersect(ctrl_samples, tmp) %>% length()
    maxvec <- c(maxvec, l)
  }
  names(maxvec) <- unique(anno[global_settings$voi] %>% dplyr::pull(., global_settings$voi))
  print(maxvec)
  new_ctrl <- names(which(maxvec == max(maxvec)))
  print(paste0("new control: ",new_ctrl))[1]
  anno[global_settings$voi][anno[global_settings$voi] == new_ctrl] <- paste0(layers_names[i], "_control")
  return(anno)
}


PCA_plot <- function(gtc = NULL, algo = NULL){
  PCAs <- list()
  PCAs[["topvar"]] <- plot_PCA_topvar()
  if(is.null(gtc)){
    if(is.null(algo)){
      for(a in names(integrated_output$alluvials)){
        PCAs[[paste0("module ", a)]] <- plot_PCA_cluster(gtc = integrated_output$alluvials[[a]], algo = a)
      }
    }else{
      PCAs[["module"]] <- plot_PCA_cluster(gtc = integrated_output$alluvials[[algo]], algo = algo)
    }
    
  }else{
    PCAs[[paste0("module my clusters")]]< - plot_PCA_cluster(gtc = gtc, algo = "my clusters")
  }
  
  return(PCAs)
}



