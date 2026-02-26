# BiocManager::install("EnsDb.Hsapiens.v86")
# BiocManager::install("ensembldb")
# BiocManager::install("AnnotationHub")
# library(magrittr)
# library(ensembldb)
# library(EnsDb.Hsapiens.v86)
# edb <- EnsDb.Hsapiens.v86
# Dna <- ensembldb::getGenomeTwoBitFile(edb)
# dna <- ensembldb::getGenomeFaFile(edb)
# transcripts <- extractTranscriptSeqs(Dna, edb)
# seq <- transcripts %>% as.data.frame()

# library('GenomicFeatures')
# library(AnnotationDbi)
################################ BiocManager::install("BSgenome")
# library(BSgenome)
################################ BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
# BiocManager::install("motifRG")

# genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")

# download.file("https://www.encodeproject.org/files/gencode.v24.primary_assembly.annotation/@@download/gencode.v24.primary_assembly.annotation.gtf.gz", "gencode.v24.primary_assembly.annotation.gtf.gz")

get_gencode_annotaiton <- function(download = F, version = 35, file = NULL){
  if(download){
    print(paste0("gencode version ", version, " is being downloaded. This might take a few minutes."))
    print(paste0("The file will be saved at: ", working_directory, "reference_files/gencode.v", version, ".primary_assembly.annotation.gtf.gz"))
    download.file(url = paste0("ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_", 
                               version, "/gencode.v", version, ".primary_assembly.annotation.gtf.gz"), destfile = paste0(working_directory, "reference_files/gencode.v", version, ".primary_assembly.annotation.gtf.gz"))
    txdb = makeTxDbFromGFF(paste0(working_directory, "reference_files/gencode.v", version, ".primary_assembly.annotation.gtf.gz"))
    saveDb(txdb, paste0(working_directory, 'reference_files/txdb.gencode', version, '.sqlite'))
    txdb = loadDb(file = paste0(working_directory, 'reference_files/txdb.gencode', version, '.sqlite'))
    
  }else{
    if(is.null(file)){
      stop("No file path defined.")
    }else{
      txdb = loadDb(file = file)
    }
  }
  return(txdb)
}


get_GRanges_transcripts <- function(txdb){
  return(GenomicFeatures::transcriptsBy(txdb))
}


get_lnc_sequences <- function(txdb, id = NULL, symbol = F, ID2SYMBOL = NULL){
  if(!symbol == F){
    if(is.null(ID2SYMBOL)){
      stop("When providing a symbol instead of an ensembl ID, you must provide ID2SYMBOL_gencode_vX_transcript.txt.")
    }else{
      if(!"symbol" %in% colnames(ID2SYMBOL) | !"ensembl" %in% colnames(ID2SYMBOL)){
        stop("ID2SYMBOL must have columns named 'symbol' (containing gene symbols) and 'ensembl' (containing ensembl gene ids).")
      }
      I2S <- ID2SYMBOL[ID2SYMBOL[,"symbol"] == symbol, ]
      id  <- I2S[1, "ensembl"]
    }
  }
  
  transcripts <- get_GRanges_transcripts(txdb)
  
  if(global_settings$organism %in% c("human", "Human")){
    genome <- BSgenome::getBSgenome("BSgenome.Hsapiens.UCSC.hg38")
  }else if(global_settings$organism %in% c("mouse", "Mouse")){
    genome <- BSgenome::getBSgenome("BSgenome.Mmusculus.UCSC.mm10")
  }else{
    stop("Organism not supported.")
  }
  
  seq_to_name <- BSgenome::getSeq(genome, transcripts[[id]]) %>% as.data.frame()
  colnames(seq_to_name) <- c("genomic_seq")
  seq_to_name$ID <- transcripts[[id]]$tx_name
  seq_to_name$rna_seq <- lapply(seq_to_name$genomic_seq, DNA_to_RNA) %>% unlist()
  return(seq_to_name)
  
}

DNA_to_RNA <- function(seq){
  seq <- strsplit(seq, split = "") %>%  unlist()
  seq <- lapply(seq, function(x){
    if(x == "A"){
      "U"
    }else if(x == "T"){
      "A"
    }else if(x == "C"){
      "G"
    }else{
      "C"
    }
  }) %>% unlist()
  
  seq <- paste0(seq, collapse = "")
  return(seq)
}

run_vienna <- function(seqs){
  for(i in 1:nrow(seqs)){
    create_vienna_input(name = seqs[i, "ID"], seq = seqs[i, "rna_seq"])
    write_cmd(fasta_file = paste0("ViennaInput_", seqs[i, "ID"], ".txt"))
    shell(cmd = paste0(working_directory, global_settings$save_folder, "/runvienna.cmd"), translate = F)
    plot_vienna_output(path = paste0(working_directory, global_settings$save_folder, "/vienna_out_", seqs[i, "ID"],".txt"))
  }
}


free_seq <- function(vienna_out){
  rnastruct <- read_delim(vienna_out, 
                          "\t", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE) %>% t()
  seq <- rnastruct[1,2]
  names(seq) <- NULL
  dotbracket <- rnastruct[1,3] %>% base::strsplit(., split = "\\(") %>% magrittr::extract2(1)
  dotbracket <- dotbracket[1:(length(dotbracket) - 1)] %>% paste0(., collapse = "(") %>% base::strsplit(., split = " ") %>% magrittr::extract2(1)
  names(dotbracket) <- NULL
  
  ct <- RRNA::makeCt(struct = dotbracket, seq = seq)
  coord <- RRNA::ct2coord(ct)
  
  freeseqs <- apply(coord, 1, function(x){
    if(as.numeric(x[9]) == 0){
      x[7]
    }else{
      "-"
    }
  }) %>% paste0(., collapse = "") %>% base::strsplit(., split = "-") %>% unlist()
  
  freeseqs <- freeseqs[!freeseqs == ""]
  
  
  return(freeseqs)
}

# Sys.setenv(PATH = paste0(Sys.getenv("PATH"), ";C:\\Program Files (x86)\\ViennaRNA Package\\"))

# txdb <- get_gencode_annotaiton(download = F, file = "C:/Users/Marie/Documents/txdb.gencode24.sqlite")

# seqs <- get_lnc_sequences(txdb = txdb, id = "ENSG00000259928.1")

# run_vienna(seqs = seqs)

# test <- free_seq(paste0(working_directory, global_settings$save_folder, "/vienna_out.txt"))



GeneToCluster <- function(){
  gtc <- do.call(rbind, apply(integrated_output$cluster_calc$cluster_information,1, function(x){
    tmp <- x["gene_n"] %>%
      base::strsplit(., split = ",")%>%
      unlist(.)
    data.frame(gene = tmp, color = rep(x["color"], length(tmp)))
  }))
  return(gtc)
}


entity_layer_cluster <- function(){
  tmp <- NULL
  for(i in 1:length(layers)){
    tmp <- rbind(tmp, data.frame(entity = V(layer_specific_outputs[[i]]$network)$name, layer = layers_names[i]))
  }
  
  gtc <- GeneToCluster()
  out <- merge(tmp, gtc, by.x = "entity", by.y = "gene")
  return(out)
}

# etl <- entity_layer_cluster()