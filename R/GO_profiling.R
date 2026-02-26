go_profiling <- function(cluster_info, level, top, ont){
  clusters <- unique(cluster_info$color)
  clusters <- clusters[!clusters == "white"]
  
  top_GO <- list()
  for(c in clusters){
    
    genes <- dplyr::filter(cluster_info, color == c)%>%
      dplyr::pull(., "gene_n")%>%
      base::strsplit(., split = ",")%>%
      BiocGenerics::unlist(.)
    enrich <- clusterProfiler::enrichGO(genes,
                                        OrgDb = "org.Hs.eg.db",
                                        keyType = "SYMBOL",
                                        ont = ont,
                                        pvalueCutoff = 0.05)
    
    if(!is.null(enrich)){
      tmp <- gofilter(enrich, level = level)@result
      tmp <- tmp[order(tmp$Count, decreasing = T),]
      top_GO[[c]] <- tmp$Description[1:top]
    }
  }
  
  top_GO <- dplyr::bind_rows(top_GO)
  if(nrow(top_GO) == 1){
    top_GO <- t(top_GO)
  }

  test <- data.frame(terms = as.vector(as.matrix(top_GO)), cluster = rep(colnames(top_GO), each = nrow(top_GO)),
                     val = factor(rep(colnames(top_GO), each = nrow(top_GO))) %>% as.numeric())
  p <- ggplot(data = test, aes(x = val, y = terms))+
    geom_point(aes(col = cluster),size = 5, pch = 15)+
    scale_color_manual(values = as.character(unique(test[order(test$val), "cluster"])))+
    scale_x_discrete(limits = as.character(unique(test[order(test$val), "cluster"])))+
    theme(axis.text.x = element_text(angle=90))+
    geom_vline(xintercept = c(1:length(unique(test$cluster))), color = as.character(unique(test[order(test$val), "cluster"])), size = 1.5)+
    xlab("")+
    ylab("")+
    ggtitle("GO enrichment")
  
  
  Cairo::CairoPDF(file = paste0(working_directory, global_settings$save_folder, "/DJplot_GO_level_", level, "_top_",
                                top, ".pdf"),
                  width = 15, height = 15)
  plot(p)
  dev.off()
  plot(p)
  output <- list()
  output[["p"]] <- p
  output[["result"]] <- top_GO
  return(output)
}


KEGG_profiling <- function(cluster_info, top){
  clusters <- unique(cluster_info$color)
  clusters <- clusters[!clusters == "white"]
  
  top_GO <- list()
  for(c in clusters){
    
    genes <- dplyr::filter(cluster_info, color == c)%>%
      dplyr::pull(., "gene_n")%>%
      base::strsplit(., split = ",")%>%
      BiocGenerics::unlist(.)
    entrez <- clusterProfiler:: bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = F)$ENTREZID
    enrich <- clusterProfiler::enrichKEGG(entrez,
                                        organism = "hsa",
                                        pvalueCutoff = 0.05)
    
    if(!is.null(enrich)){
      tmp <- enrich@result
      tmp <- tmp[order(tmp$Count, decreasing = T),]
      top_GO[[c]] <- tmp$Description[1:top]
    }
    
  }
  top_GO <- dplyr::bind_rows(top_GO)
  if(nrow(top_GO) == 1){
    top_GO <- t(top_GO)
  }

  test <- data.frame(terms = as.vector(as.matrix(top_GO)), cluster = rep(colnames(top_GO), each = nrow(top_GO)),
                     val = factor(rep(colnames(top_GO), each = nrow(top_GO))) %>% as.numeric())
  p <- ggplot(data = test, aes(x = val, y = terms))+
    geom_point(aes(col = cluster),size = 5, pch = 15)+
    scale_color_manual(values = as.character(unique(test[order(test$val), "cluster"])))+
    scale_x_discrete(limits = as.character(unique(test[order(test$val), "cluster"])))+
    theme(axis.text.x = element_text(angle=90))+
    geom_vline(xintercept = c(1:length(unique(test$cluster))), color = as.character(unique(test[order(test$val), "cluster"])), size = 1.5)+
    xlab("")+
    ylab("")+
    ggtitle("KEGG enrichment")
  
  Cairo::CairoPDF(file = paste0(working_directory, global_settings$save_folder, "/DJplot_KEGG_top_",
                                top, ".pdf"),
                  width = 15, height = 15)
  plot(p)
  dev.off()
  plot(p)
  output <- list()
  output[["p"]] <- p
  output[["result"]] <- top_GO
  return(output)
}


HALLMARK_profiling <- function(cluster_info, top, hm){
  clusters <- unique(cluster_info$color)
  clusters <- clusters[!clusters == "white"]
  
  top_GO <- list()
  for(c in clusters){
    
    genes <- dplyr::filter(cluster_info, color == c)%>%
      dplyr::pull(., "gene_n")%>%
      base::strsplit(., split = ",")%>%
      BiocGenerics::unlist(.)
    enrich <- clusterProfiler::enricher(genes,
                                        TERM2GENE=hm,
                                        pvalueCutoff = 0.05,
                                        pAdjustMethod = "none",
                                        qvalueCutoff = 1.0)
    
    if(!is.null(enrich)){
      tmp <- enrich@result
      tmp <- tmp[order(tmp$Count, decreasing = T),]
      top_GO[[c]] <- tmp$Description[1:top]
    }
    
  }
  
  top_GO <- dplyr::bind_rows(top_GO)
  if(nrow(top_GO) == 1){
    top_GO <- t(top_GO)
  }
  test <- data.frame(terms = as.vector(as.matrix(top_GO)), cluster = rep(colnames(top_GO), each = nrow(top_GO)),
                     val = factor(rep(colnames(top_GO), each = nrow(top_GO))) %>% as.numeric())
  p <- ggplot(data = test, aes(x = val, y = terms))+
    geom_point(aes(col = cluster),size = 5, pch = 15)+
    scale_color_manual(values = as.character(unique(test[order(test$val), "cluster"])))+
    scale_x_discrete(limits = as.character(unique(test[order(test$val), "cluster"])))+
    theme(axis.text.x = element_text(angle=90))+
    geom_vline(xintercept = c(1:length(unique(test$cluster))), color = as.character(unique(test[order(test$val), "cluster"])), size = 1.5)+
    xlab("")+
    ylab("")+
    ggtitle("HALLMARK enrichment")
  
  Cairo::CairoPDF(file = paste0(working_directory, global_settings$save_folder, "/DJplot_HALLMARK_top_",
                                top, ".pdf"),
                  width = 15, height = 15)
  plot(p)
  dev.off()
  plot(p)
  output <- list()
  output[["p"]] <- p
  output[["result"]] <- top_GO
  return(output)
}



