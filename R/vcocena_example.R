#' Example data paths for vcocena
#'
#' Returns file paths for the small example dataset included with the package.
#' Use these paths to build `layers` and `supplement` inputs for
#' \\code{run_vcocena()}.
#'
#' @return A named list with paths to example count tables, annotations,
#'   and reference files.
#' @export
vcocena_example <- function() {
  base_dir <- system.file("extdata", "example_vertical", package = "vcocena")
  if (base_dir == "") {
    stop("Example data not found. Is the package installed correctly?")
  }

  list(
    working_directory = base_dir,
    counts_seq = "data_seq_processed.txt",
    counts_array = "data_array_processed.txt",
    anno_seq = "annotation_seq.txt",
    anno_array = "annotation_array.txt",
    reference_dir = file.path(base_dir, "data", "reference_files"),
    tf = "TFcat.txt",
    hallmark = "toy_hallmark.gmt",
    go = "toy_go.gmt",
    kegg = "toy_kegg.gmt",
    reactome = "toy_reactome.gmt"
  )
}
