# plot_cutoffs_internal <- function(cutoff_stats, 
#                                   hline = list("R.squared" = NULL, "no_edges" = NULL, "no_nodes" = NULL, "no_networks" = NULL),
#                                   x){
#   cutoff_stats$corr <- rownames(cutoff_stats) %>% as.numeric()
#   
#   
#   p1 <- ggplot(cutoff_stats, aes(x = corr))+
#     geom_vline(xintercept = c(seq(from = min(cutoff_stats$corr), to = 1, by = 0.01)), color = "darkgrey")+
#     geom_hline(yintercept = hline[[1]], color = "darkgrey", size = 1.7)+
#     geom_line(aes(y = R.squared), color = "#ff00b3")+
#     geom_point(aes(y = R.squared), color = "#ff00b3", size = 2)+
#     theme_light()+
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank())+
#     scale_x_continuous(breaks = seq(from = min(cutoff_stats$corr), to = 1, by = 0.005))+
#     scale_y_continuous(labels = comma)+
#     ggtitle(layers_names[x])
#   
#   
#   p2 <- ggplot(cutoff_stats, aes(x = corr))+
#     geom_vline(xintercept = c(seq(from = min(cutoff_stats$corr), to = 1, by = 0.01)), color = "darkgrey")+
#     geom_hline(yintercept = hline[[2]], color = "darkgrey", size = 1.7)+
#     geom_line(aes(y = no_edges), color = "#00e5ff")+ 
#     geom_point(aes(y = no_edges), color = "#00e5ff", size = 2)+
#     theme_light()+ 
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank())+
#     scale_x_continuous(breaks = seq(from = min(cutoff_stats$corr), to = 1, by = 0.005))+
#     scale_y_continuous(labels = comma)
#   
#   
#   p3 <- ggplot(cutoff_stats, aes(x = corr))+
#     geom_vline(xintercept = c(seq(from = min(cutoff_stats$corr), to = 1, by = 0.01)), color = "darkgrey")+
#     geom_hline(yintercept = hline[[3]], color = "darkgrey", size = 1.7)+
#     geom_line(aes(y = no_nodes), color = "#20d404")+ 
#     geom_point(aes(y = no_nodes), color = "#20d404", size = 2)+
#     theme_light()+ 
#     theme(axis.title.x = element_blank(),
#           axis.text.x = element_blank(),
#           axis.ticks.x = element_blank())+
#     scale_x_continuous(breaks = seq(from = min(cutoff_stats$corr), to = 1, by = 0.005))+
#     scale_y_continuous(labels = comma)
#   
# 
#   
#   
#   p5 <- ggplot(cutoff_stats, aes(x = corr))+
#     geom_vline(xintercept = c(seq(from = min(cutoff_stats$corr), to = 1, by = 0.01)), color = "darkgrey")+
#     geom_hline(yintercept = hline[[4]], color = "darkgrey", size = 1.7)+
#     geom_point(aes(y = no_of_networks), color = "#db0202", size = 2)+
#     scale_x_continuous(breaks = seq(from = min(cutoff_stats$corr), to = 1, by = 0.005))+
#     theme_light()+
#     theme(axis.text.x = element_text(angle=90))
#   
#   
#   p <- cowplot::plot_grid(p1, p2, p3, p5, ncol = 1, align = "v")
#   return(p)
# }


plot_cutoffs_internal <- function(cutoff_stats, 
                                  hline = list("R.squared" = NULL, "no_edges" = NULL, "no_nodes" = NULL, "no_networks" = NULL),
                                  x){
  cutoff_stats$corr <- rownames(cutoff_stats) %>% as.numeric()
  
  
  p1 <- plot_ly(cutoff_stats, x = ~corr, y = ~R.squared, type = 'scatter', 
                mode = 'lines+markers', name = "R", line = list(color = "lightblue"), marker = list(color = "lightblue")) 
  p2 <- plot_ly(cutoff_stats, x = ~corr, y = ~no_edges, type = 'scatter', 
                mode = 'lines+markers', name = "no. edges", line = list(color = "orange"), marker = list(color = "orange"))
  p3 <- plot_ly(cutoff_stats, x = ~corr, y = ~no_nodes, type = 'scatter', 
                mode = 'lines+markers', name = "no. nodes", line = list(color = "lightgreen"), marker = list(color = "lightgreen"))
  p4 <- plot_ly(cutoff_stats, x = ~corr, y = ~no_of_networks, type = 'scatter', 
                mode = "markers", name = "no. networks", marker = list(color = "yellow"))
  p <- plotly::subplot(p1, p2, p3, p4, nrows = 4, shareX = T)
  
  steps <- list()
  for(i in 1:length(cutoff_stats$corr)){
    
    step <- list(args = list("marker.color",list(rep("lightblue", length(cutoff_stats$corr)),
                                                 rep("orange", length(cutoff_stats$corr)),
                                                 rep("lightgreen", length(cutoff_stats$corr)),
                                                 rep("yellow", length(cutoff_stats$corr)))), 
                 label = paste0(as.character(cutoff_stats$corr[i]), ", R?: ", 
                                round(cutoff_stats$R.squared[i], 3),
                                "; no. edges: ", cutoff_stats$no_edges[i], "; no. nodes: ", 
                                cutoff_stats$no_nodes[i], "; no. networks: ", cutoff_stats$no_of_networks[i]), 
                 method = "restyle"
                 
    )
    # print(step)
    for(j in 1:4){
      step[["args"]][[2]][[j]][i] <- "red"
    }
    
    steps[[i]] <- step
  }
  
  
  p <- p %>% layout(hovermode = "x unified")%>%
    layout(title = paste0("Cut-off selection guide: ", layers_names[x]),
           sliders = list(
             list(pad = list(t=60),
                  active = 2, 
                  currentvalue = list(prefix = "Cut-off: ", font = list(color = "black", size = 14)), 
                  steps = steps,
                  font = list(color = "white", size = 0))),
           shapes = list(line))
  return(p)
}

plot_cutoffs <- function(
  hline = list(
    "R.squared"   = NULL,
    "no_edges"    = NULL,
    "no_nodes"    = NULL,
    "no_networks" = NULL
  )
) {
  p <- NULL
  for (x in seq_along(layers)) {
    cutoff_df <- layer_specific_outputs[[paste0("set", x)]]$part1$cutoff_calc_out$cutoff_stats_concise
    p <- plot_cutoffs_internal(
      cutoff_stats = cutoff_df,
      hline        = hline,
      x            = x
    )

    # Save as PDF if save folder is defined
    if (!is.null(global_settings$save_folder)) {
      pdf_path <- file.path(
        working_directory,
        global_settings$save_folder,
        paste0("cutoff_plot_", layers_names[x], ".pdf")
      )
      Cairo::CairoPDF(file = pdf_path, width = 15, height = 8)
      print(p)
      dev.off()
    }
  }
  invisible(p)
}
