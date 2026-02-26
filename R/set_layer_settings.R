set_layer_settings <- function(num_layers = length(layers), 
                               top_var, 
                               min_corr, 
                               range_cutoff_length, 
                               print_distribution_plots){
  
  for(x in 1:num_layers){
    layer_settings[[paste0("set", x)]] <- list()
    layer_settings[[paste0("set", x)]][["top_var"]] <- top_var[x]
    layer_settings[[paste0("set", x)]][["min_corr"]] <- min_corr[x]
    layer_settings[[paste0("set", x)]][["range_cutoff_length"]] <- range_cutoff_length[x]
    layer_settings[[paste0("set", x)]][["print_distribution_plots"]] <- print_distribution_plots[x]
  }
  
  return(layer_settings)
}