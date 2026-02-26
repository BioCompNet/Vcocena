antilog <- function(lx , base) {
  
  lbx <- lx/log(exp(1) , base = base)
  result <- exp(lbx)
  result
}

gfc_calc <- function(grp, trans_norm, group_means){
  df1= trans_norm[,grp]
  df2=gtools::foldchange(df1, group_means) %>% ifelse(.>global_settings$range_GFC, global_settings$range_GFC,.) %>% 
    ifelse(.< (-global_settings$range_GFC), -global_settings$range_GFC,.) %>%
    as.data.frame()
  colnames(df2) = paste0("",grp)
  return(df2)
}

GFC_calculation <- function(info_dataset, global_set, layer_set,input1, input2, grouping_v) {
  
  print(paste0("...calculate Group-Fold-Changes..."))
  
  if(!is.null(grouping_v)){
    info_dataset$grpvar <- info_dataset[,c(grouping_v)]
    print(paste0("User-defined variable ", grouping_v, " will be used for grouping the data."))
  }
  else if(intersect(global_set$voi,colnames(info_dataset))%>% length()>0){
    
    print(paste("variables:",paste0(intersect(global_set$voi,colnames(info_dataset)), collapse = ","),
                "will be used as grouping variables,",
                "if these are different in # from voi_id then please check if",
                "the variables mentioned in voi_id are present in the metadata"))
    
    info_dataset$grpvar =purrr::pmap(input2$heatmap_out$corresp_info[intersect(global_set$voi,colnames(info_dataset))],
                                                        paste, sep="-") %>% unlist()
  } else {
    
    print(paste("the first column in the metadata will be used as the grouping variable",
                "since the voi_id is not present in the metadata"))
    info_dataset$grpvar = info_dataset[,1]
  }

  # Safety fallback: ensure a grouping column exists.
  # If grpvar was not created above (e.g., due to mismatched column names),
  # default to the first variable of interest when available.
  if (!"grpvar" %in% colnames(info_dataset)) {
    if (!is.null(global_set$voi) && global_set$voi[1] %in% colnames(info_dataset)) {
      info_dataset$grpvar <- info_dataset[[global_set$voi[1]]]
    } else {
      # As a last resort, use the first column.
      info_dataset$grpvar <- info_dataset[[1]]
    }
  }
  
  
  if(global_settings$control == "none"){
    # if(length(layers) > 1){
    #   stop("You provided no controls but more than one dataset. When analysing serveral datasets, each must have a control to guarantee
    #         comparibility between datasets.")
    # }else{
      print("GFC calculation with foldchange from mean")
      # GFC calculation with foldchange from mean
      norm_data_anno = merge(t(input1$topvar), subset(info_dataset, select = c("grpvar")), by="row.names", all.x=T)
      
      norm_data_anno = norm_data_anno[,-1]
      
      norm_data_anno <- norm_data_anno[ , c(ncol(norm_data_anno) , 1:(ncol(norm_data_anno)-1))]
      
      trans_norm <- setNames(data.frame(t(norm_data_anno[ , -1])) , norm_data_anno[,1])

      if (isTRUE(global_set$data_in_log)) {
        trans_norm <- antilog(trans_norm , 2)
      }
      
      trans_norm <- t(apply(trans_norm , 1 , function(x) tapply(x , colnames(trans_norm) , mean)))
      trans_norm <- cbind(trans_norm , rowMeans(trans_norm))
      colnames(trans_norm)[ncol(trans_norm)] <- "group_mean"
      grplist=colnames(trans_norm)[-(ncol(trans_norm))]
      
      
      GFC_all_genes <- do.call(
        "cbind",
        lapply(grplist, gfc_calc, trans_norm = trans_norm, group_means = trans_norm[, "group_mean"])
      )
      GFC_all_genes <- as.data.frame(GFC_all_genes)
      GFC_all_genes[] <- lapply(GFC_all_genes, function(col) {
        if (is.numeric(col)) round(col, 3) else col
      })
      GFC_all_genes$Gene <- rownames(GFC_all_genes)
      
      
      return(GFC_all_genes)
      
    # }
  }else{
    print("GFC calculation with foldchange from control")
    # GFC calculation with foldchange from control
    norm_data_anno = merge(t(input1$topvar), subset(info_dataset, select = c("grpvar")), by="row.names", all.x=T)
    #contains a column called Row.names (1st col)
    norm_data_anno = norm_data_anno[,-1]

    norm_data_anno <- norm_data_anno[ , c(ncol(norm_data_anno) , 1:(ncol(norm_data_anno)-1))]
    
    trans_norm <- setNames(data.frame(t(norm_data_anno[ , -1])) , norm_data_anno[,1])

    if (isTRUE(global_set$data_in_log)) {
      trans_norm <- antilog(trans_norm , 2)
    }
    
    trans_norm <- t(apply(trans_norm , 1 , function(x) tapply(x , colnames(trans_norm) , mean)))
    
    
    
    trans_norm_no_ctrl <- trans_norm[, grepl(global_settings$control, colnames(trans_norm), ignore.case = T) == F]
    if(is.vector(trans_norm_no_ctrl)){
      
      trans_norm_no_ctrl <- data.frame(trans_norm_no_ctrl = trans_norm_no_ctrl)
      
      colnames(trans_norm_no_ctrl) <- colnames(trans_norm)[grepl(global_settings$control, 
                                                                 colnames(trans_norm), ignore.case = T)==F]
      rownames(trans_norm_no_ctrl) <- rownames(trans_norm)
    }
    trans_norm <- cbind(trans_norm_no_ctrl, trans_norm[, grepl(global_settings$control, 
                                                               colnames(trans_norm), ignore.case = T) == T])
    rownames(trans_norm) <- rownames(trans_norm_no_ctrl)
    colnames(trans_norm)[ncol(trans_norm)] <- "group_mean"
    grplist=colnames(trans_norm)[-(ncol(trans_norm))]
    
    
    zero_rownames <- trans_norm[trans_norm[, ncol(trans_norm)]==0,] %>% rownames()
    # print(head(zero_rownames))
    
    # apply(trans_norm,2, function(x){
    #   print(0 %in% x)
    # })
    # print("stop")
    GFC_all_genes <- do.call(
      "cbind",
      lapply(grplist, gfc_calc, trans_norm = trans_norm, group_means = trans_norm[, "group_mean"])
    )
    # apply(GFC_all_genes,2, function(x){
    #   print(NA %in% x)
    # })
    # test <- cbind(GFC_all_genes, trans_norm[, ncol(trans_norm)])
    # GFC_NA <- test[is.na(test),] 
    GFC_all_genes <- as.data.frame(GFC_all_genes)
    GFC_all_genes[] <- lapply(GFC_all_genes, function(col) {
      if (is.numeric(col)) round(col, 3) else col
    })
    GFC_all_genes[is.na(GFC_all_genes)] <- 0
    # GFC_NA2 <- GFC_all_genes[is.na(GFC_all_genes),] %>% rownames()
    # apply(GFC_all_genes,2, function(x){
    #   print(NA %in% x)
    # })
    rownames(GFC_all_genes) <- rownames(trans_norm)
    GFC_all_genes$Gene = rownames(GFC_all_genes)
    tmp_col_names <- colnames(GFC_all_genes)[grepl(global_settings$control, 
                                                   colnames(GFC_all_genes), ignore.case = T)==F]
    GFC_all_genes <- GFC_all_genes[, colnames(GFC_all_genes) %in% tmp_col_names]
    
    return(GFC_all_genes)
  }
  
  
}
