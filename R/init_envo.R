
#source(paste0(working_directory,"scripts/", "make_rownames_unique.R"))

# function to initialize the save folder
init_save_folder <- function(name){
  
  save_folder <- name
  
  if(!save_folder %in% list.dirs(working_directory)) {
    
    dir.create(paste0(working_directory,save_folder))
    
  }
  
  return(name)
}



make_rownames_unique <- function(counts_new, gene_symbol_col){
  counts_new <- counts_new[!duplicated(counts_new[gene_symbol_col]),] %>%
    tibble::remove_rownames() %>%
    tibble::column_to_rownames(., gene_symbol_col)
  
  # remove all non-numeric columns (description, gene id, etc.):
  for (x in colnames(counts_new)){
    
    if(!is.numeric(counts_new[[x]])){
      
      counts_new[[x]] <- NULL
      
    }
    
  }
  
  if(ncol(counts_new) == 0){
    
    print("All columns were deleted when removing non-numeric columns. Please check the data type of your expression values.")
    
  }
  
  return(counts_new)
}



# fucntion to load expression data:
read_expression_data <- function(file, rown = T, sep = "\t", gene_symbol_col){
  
  if(rown){
    
    expression_data <- read.table(file = file, row.names = 1,
                                  stringsAsFactors = F, sep = sep, check.names = F, header = T)
    
  }else{
    
    expression_data <- read.table(file = file, 
                                  stringsAsFactors = F, sep = sep, check.names = F, header = T)
  }
  expression_data <- make_rownames_unique(expression_data, gene_symbol_col = gene_symbol_col)
  
  return(expression_data)
}



# fucntion to load annotation data:
read_anno <- function(file, rown = T, sep = "\t", sample_col){
  
  if(rown){
    
    anno <- read.table(file = file, row.names = 1,
                       stringsAsFactors = F, sep = sep, check.names = T, header = T) 
    
  }else{
    
    anno <- read.table(file = file, 
                       stringsAsFactors = F, sep = sep, check.names = T, header = T)
    
  }
  rownames(anno) <- dplyr::pull(anno, sample_col) 
  
  anno[] <- lapply(anno, base::factor)
  
  return(anno)
}



loading_data <- function(wd, layers, sep_counts , sep_anno, gene_symbol_col, sample_col, count_has_rn, anno_has_rn){
  
  data <- list()
  
  for (x in 1:length(layers)){
    if(exists(layers[[paste0("set",x)]][1])){
      
      data[[paste0("set", x, "_counts")]] <- get(layers[[paste0("set",x)]][1])
      
    }else{
      
      data[[paste0("set", x, "_counts")]] <- read_expression_data(file = paste0(wd, "data/", layers[[paste0("set",x)]][1]), rown = count_has_rn,
                                                                  sep = sep_counts, gene_symbol_col = gene_symbol_col)
    }
    if(exists(layers[[paste0("set",x)]][2])){
      anno <- get(layers[[paste0("set",x)]][2])
      
      anno[] <- lapply(anno, base::factor)
      
      data[[paste0("set", x, "_anno")]] <- anno
      
    }else{
      
      data[[paste0("set", x, "_anno")]] <- read_anno(file = paste0(wd,"sample_info/", layers[[paste0("set",x)]][2]), rown = anno_has_rn,
                                                     sep = sep_anno, sample_col = sample_col)
      
    }
    
    if(!ncol(data[[paste0("set", x, "_counts")]]) == nrow(data[[paste0("set", x, "_anno")]])){
      stop(paste0("The count table has ",  ncol(data[[paste0("set", x, "_counts")]]), " columns but the annotation has ", 
                  nrow(data[[paste0("set", x, "_anno")]]), " rows. These values are required to be the same sice they
                     should correspond to the number of samples. THE LOADING OF THE DATA WILL BE TERMINATED."))
    }
    
  }
  return(data)
}



initialize_environment <- function(sep_counts = ",", sep_anno = ",", gene_symbol_col, sample_col, count_has_rn, anno_has_rn){
  
  knitr::opts_knit$set(root.dir = working_directory)
  options(dplyr.summarise.inform = F)
  data <- loading_data(wd = working_directory, layers = layers, sep_counts = sep_counts, sep_anno = sep_counts, 
                       gene_symbol_col = gene_symbol_col, sample_col = sample_col, count_has_rn = count_has_rn, anno_has_rn = anno_has_rn)
  
  return(data)
}



load_supplementary_data <- function(){
  
  data <- list()
  
  data[["TF"]] <- read.delim(paste0(working_directory, "data/reference_files/", supplement[1]),
                             header=TRUE,
                             check.names=F)
  
  data[["hallmark"]] <- clusterProfiler::read.gmt(paste0(working_directory, "data/reference_files/", supplement[2]))
  
  data[["go"]] <- clusterProfiler::read.gmt(paste0(working_directory, "data/reference_files/", supplement[3]))
  
  data[["kegg"]] <- clusterProfiler::read.gmt(paste0(working_directory, "data/reference_files/", supplement[4]))
  
  data[["reactome"]] <- clusterProfiler::read.gmt(paste0(working_directory, "data/reference_files/", supplement[5]))
  
  return(data)
}
