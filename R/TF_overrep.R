makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
                                              blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}


TF_overrep <- function(topTF = 5, topTarget = 5){
  output <- list()
  gtc <- GeneToCluster() 
  colnames(gtc) <- c("gene", "cluster")
  for(c in unique(integrated_output$cluster_calc$cluster_information$color)){
    if(!c=="white"){
      genes = dplyr::filter(integrated_output$cluster_calc$cluster_information, color == c)%>%
        dplyr::pull(., "gene_n")%>%
        base::strsplit(., split = ",") %>%
        unlist(.)
      
      url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
      encode = "json"
      payload = list(query_name = "myQuery", gene_set = genes)
      
      #POST to ChEA3 server
      response = POST(url = url, body = payload, encode = encode)
      json = httr::content(response, as = "text")
      
      #results as list of R dataframes
      results = fromJSON(json)
      results <- results$`Integrated--meanRank`
      gtc_not_white <- gtc[!gtc$cluster == "white",]
      results <- dplyr::filter(results, TF %in% gtc_not_white$gene)
      # extract those from meanRank since meanRank scored as best method:
      resultlist <- list()
      for(i in 1:topTF){
        tf <- results$TF[i]
        overlapping_genes <- results$Overlapping_Genes[i]%>%
          base::strsplit(., split = ",")%>%
          unlist(.)
        resultlist[[tf]] <- list(TF = tf, targets = overlapping_genes[1:topTarget])
      }
      
      # filter for only those that are present in our network:
      output[[c]] <- resultlist
    }
  }

  return(output)
}


plot_TF_enrichment <- function(overrep_results){
  gtc <- GeneToCluster() 
  colnames(gtc) <- c("gene", "cluster")
  dflist <- list()
  TFs <- NULL
  for(c in names(overrep_results)){
    tmp <- NULL
    for(t in names(overrep_results[[c]])){
      TFs <- c(TFs, t)
      clt <- dplyr::filter(gtc, gene == t)%>%
        dplyr::pull(., "cluster")
      tmp_df <- data.frame(TF = rep(t, length(overrep_results[[c]][[t]][["targets"]])), 
                           Target = overrep_results[[c]][[t]][["targets"]],
                           ClusterTF = rep(clt, length(overrep_results[[c]][[t]][["targets"]])))
      #tmp_df <- merge(tmp_df, gtc, by.x = "Target", by.y = "gene")
      colnames(tmp_df) <- c("TF", "Target", "ClusterTF")
      tmp <- rbind(tmp, tmp_df)
    }
    dflist[[c]] <- tmp
  }
  
  
  TFs <- unique(as.character(TFs))
  edgelist <- integrated_output$combined_edgelist
  edgelist$merged <- paste0(as.character(edgelist$V1), as.character(edgelist$V2))
  edgelist$merged2 <- paste0(as.character(edgelist$V2), as.character(edgelist$V1))
  
  Cairo::CairoPDF(file = paste0(working_directory, global_settings$save_folder, "/circos_plots.pdf"), width = 12, 
                    height = 7, onefile = T)
  # loop in pairs of two to create douple plots: 
  for(n in 1:round((length(dflist)+1/2))){
    i1 <- n*2-1
    i2 <- n*2
    # set layout to plot two plots each in horizontal arrangement
    layout(matrix(1:2, 1, 2)) 
    
    for(j in c(i1, i2)){
      # catching 'out-of-bounds':
      if(j > length(dflist)){
        break
      }
      # create link data frame:
      fromto <- dflist[[j]]
      # only consider no TF targets:
      fromto <- dplyr::filter(fromto, !Target %in% TF)
      fromto <- fromto[complete.cases(fromto),]
      fromto <- unique(fromto)
      
      
      # dataframe that associates each gene with its clsuter colour:
      NodeToColor <- rbind(data.frame(gene = fromto$TF, color = fromto$ClusterTF),
                           data.frame(gene = fromto$Target, color = rep(names(dflist)[j], nrow(fromto))))%>%
        unique()
      
      # create plot factors:
      factors <- unique(as.character(NodeToColor$gene))
      circos.par(points.overflow.warning=FALSE)
      circos.initialize(factors, xlim = c(0, 1)) 
      circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, bg.col = as.character(NodeToColor$color),
                             bg.border = NA ) 
      # add sector labels:
      g <- circos.trackPlotRegion(track.index = 1, panel.fun = function(x,y){
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        sector.name = get.cell.meta.data("sector.index")
        if(sector.name %in% TFs){
          circos.text(mean(xlim), mean(ylim)+2.5, sector.name, facing = "inside", niceFacing = T, cex = .9, 
                      col = "turquoise3", font = 2)
        }else{
          circos.text(mean(xlim), mean(ylim)+2.5, sector.name, facing = "inside", niceFacing = T, cex = .9)
        }
        
      })
      # add links
      for(i in 1:nrow(fromto)) {
        merged <- paste0(as.character(fromto[i,1]), as.character(fromto[i,2]))
        if(merged %in% edgelist$merged | merged %in% edgelist $merged2){
          g <- circos.link(sector.index1 =  as.character(fromto[i,1]), c(0.45, 0.55),
                           sector.index2 =  as.character(fromto[i,2]), c(.92), 
                           col = as.character(fromto[i,3]),
                           directional = 1,
                           arr.width = .1,
                           arr.length = .1)
        }else{
          g <- circos.link(sector.index1 =  as.character(fromto[i,1]), c(0.45, 0.55),
                           sector.index2 =  as.character(fromto[i,2]), c(.92), 
                           col = makeTransparent(as.character(fromto[i,3]), alpha = 30),
                           directional = 1,
                           arr.width = .1,
                           arr.length = .1)
        }
        
        
      }
      title(names(dflist)[j])
      circos.clear()
    }
  }
  dev.off()
}
  



TF_enrich_all <- function(topTF = 100, topTarget = 30){
  gtc <- GeneToCluster() 
  colnames(gtc) <- c("gene", "cluster")
  genes <- gtc$gene
  
  
  url = "https://amp.pharm.mssm.edu/chea3/api/enrich/"
  encode = "json"
  payload = list(query_name = "myQuery", gene_set = genes)
  
  #POST to ChEA3 server
  response = POST(url = url, body = payload, encode = encode)
  json = httr::content(response, as = "text")
  
  #results as list of R dataframes
  results = fromJSON(json)
  results <- results$`Integrated--meanRank`
  gtc_not_white <- gtc[!gtc$cluster == "white",]
  results <- dplyr::filter(results, TF %in% gtc_not_white$gene)
  # extract those from meanRank since meanRank scored as best method:
  resultlist <- list()
  for(i in 1:topTF){
    tf <- results$TF[i]
    overlapping_genes <- results$Overlapping_Genes[i]%>%
      base::strsplit(., split = ",")%>%
      unlist(.)
    overlapping_genes <- overlapping_genes[overlapping_genes %in% gtc_not_white$gene]
    resultlist[[tf]] <- list(TF = tf, targets = overlapping_genes[1:topTarget])
  }
  return(resultlist)
}


check_TF <- function(TF, TF_Enrichment_Network){
  edgelist <- integrated_output$combined_edgelist[,c(1,2)]
  edgelist[] <- lapply(edgelist, as.character)
  edgelist$merged <- paste0(edgelist$V1, edgelist$V2)
  edgelist$merged2 <- paste0(edgelist$V2, edgelist$V1)
  
  gtc <- GeneToCluster()
  colnames(gtc) <- c("gene", "cluster")
  
  targets <- TF_Enrichment_Network[[TF]][["targets"]]
  
  edges <- data.frame(from = rep(TF, length(targets)), to = targets) %>% 
    unique()
  edges[] <- lapply(edges, as.character)
  edges <- edges[!edges$to == TF,]
  merged <- paste0(edges$from, edges$to)
  edges$color <- lapply(merged, function(x){
    if(x %in% edgelist$merged | x %in% edgelist$merged2){
      "black"
    }else{
      "grey"
    }
  })%>% unlist()
  
  nodes <- data.frame(name = unique(c(edges$from, edges$to)))
  nodes <- merge(nodes, gtc, by.x = "name", by.y = "gene") 
  colnames(nodes) <- c("name", "color")
  nodes <- nodes[order(nodes$color),]
  
  g <- igraph::graph_from_data_frame(d = edges, vertices = nodes$name)
  V(g)$color <- as.character(nodes$color)
  l <- igraph::layout.star(g, center = V(g)[TF])
  Cairo::CairoPDF(file = paste0(working_directory, global_settings$save_folder, "/testplot.pdf"), width = 10,
                                height = 10)
  igraph::plot.igraph(g, layout = l, edge.arrow.size = 0.5, vertex.label.color = "black", edge.color = edges$color,
                      vertex.label.cex = 0.7, vertex.label.font = 2, edge.width = 2, 
                      vertex.frame.color = as.character(nodes$color))
  dev.off()
  }


