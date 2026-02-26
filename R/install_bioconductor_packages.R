install_bioconductor_packages <- function(install){
  if(install){
    BiocManager::install(c("biomaRt",
                           "BiocParallel",
                           "ComplexHeatmap",
                           "clusterProfiler",
                           "DOSE",
                           "enrichplot",
                           "org.Hs.eg.db",
                           "pcaGoPromoter", 
                           "pcaGoPromoter.Hs.hg19",
                           "RCy3",
                           "ReactomePA",
                           "STRINGdb",
                           "org.Mm.eg.db",
                           "pcaGoPromoter.Mm.mm9"))
  }
}
