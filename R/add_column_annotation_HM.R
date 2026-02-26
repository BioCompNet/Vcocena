get_anno_matrix <- function(variables){
  all_mats <- list()
  for(v in 1:length(variables)){
    if(is.na(variables[v])){
      names <- dplyr::pull(data[[paste0("set", v, "_anno")]], global_settings$voi) %>% unique()
      layer_mat <- matrix(rep(0,length(names[names %in% colnames(integrated_output$GFC_all_layers)])), ncol = 1)
      rownames(layer_mat) <- names[names %in% colnames(integrated_output$GFC_all_layers)]
    }else{
      anno <- data[[paste0("set", v, "_anno")]]
      layer_mat <- NULL
      for(i in colnames(integrated_output$GFC_all_layers)[!colnames(integrated_output$GFC_all_layers) == "Gene"]){
        
        if(i %in% dplyr::pull(anno, global_settings$voi)){
          tmp <- dplyr::filter(anno, anno[global_settings$voi] == i)%>%
            dplyr::pull(., variables[v])
          
          tmp <- table(tmp)%>%
            t()%>%
            as.data.frame()
          
          tmp$Var1 <- NULL
          colnames(tmp) <- c("var", "freq")
          mat <- matrix(tmp$freq, nrow = 1, byrow = T)
          colnames(mat) <- tmp$var
          rownames(mat) <- i
          layer_mat <- bind_rows(as.data.frame(layer_mat),  as.data.frame(mat))
        }
      }
    }
    layer_mat[is.na(layer_mat)] <- 0
    layer_mat <- as.matrix(layer_mat)
    all_mats[[v]] <- layer_mat
  }
  return(all_mats)
}


unify_mats <- function(mat_list){
  merged_mat <- NULL
  new_rn <- NULL
  for(x in mat_list){
    if(!is.null(merged_mat)){
      if(all(colnames(merged_mat) %in% colnames(x))){
        new_rn <- c(new_rn, rownames(x))
        x <- x[, colnames(merged_mat)]
        merged_mat <- rbind(merged_mat,x)
      }else{
        new_rn <- c(new_rn, rownames(x))
        merged_mat <- merge(merged_mat, x, by = "row.names", all = T)
        merged_mat$Row.names <- NULL
      }
    }else{
      new_rn <- c(new_rn, rownames(x))
      merged_mat <- merge(merged_mat, x, by = "row.names", all = T)
      merged_mat$Row.names <- NULL
    }
    
  }
  merged_mat[is.na(merged_mat)] <- 0
  merged_mat <- merged_mat[, colSums(merged_mat)>0]
  rownames(merged_mat) <- new_rn
  return(merged_mat)
}

col_anno_categorical <- function(variables){
  ml <- get_anno_matrix(variables = variables)
  m <- unify_mats(ml)
  m <- apply(m, 1, function(x){
    if(sum(x) == 0){
      x
    }else{
      x/sum(x)*100
    }
  })%>%
    t()%>%
    as.data.frame()
  
  m[rowSums(m) == 0,] <- NA
  return(m)
}





col_anno_numerical <- function(variables){
  vals <- list()
  for(i in 1:length(variables)){
    anno <- data[[paste0("set", i, "_anno")]]
    if(!global_settings$control == "none"){
      ctrl <- unique(dplyr::pull(anno, global_settings$voi))[grepl(global_settings$control, 
                                                                   unique(dplyr::pull(anno, global_settings$voi)), 
                                                                   ignore.case = T) ==T]
    }
    for(j in unique(dplyr::pull(anno, global_settings$voi))){
      if(j %in% ctrl){
        next
      }
      if(is.na(variables[i])){
        vals[[j]] <- NA
      }else{
        tmp <- anno[anno[global_settings$voi] == j, ]
        vals[[j]] <- as.numeric(as.character(dplyr::pull(tmp, variables[i])))
      }
    }
  }
  return(vals)
}

