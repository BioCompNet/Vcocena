## High-level wrapper for running the Vertical CoCena² pipeline

#' Run the Vertical CoCena² pipeline
#'
#' This function mirrors the main steps of `vCoCena.Rmd` and orchestrates
#' data loading, per-layer correlation / cutoff analysis, GFC calculation,
#' network construction, clustering, and vertical integration across layers.
#'
#' @param working_directory Path to the project directory that contains
#'   `data/`, `sample_info/`, and `reference_files/`.
#' @param layers Named list describing omics layers.
#'   Each entry should be a character vector of length 2:
#'   \code{c(count_file, annotation_file)}.
#' @param supplement Character vector of length 5 with file names (or paths
#'   relative to \code{reference_files/}) for TF list, Hallmark, GO, KEGG,
#'   and Reactome GMT files.
#' @param layers_names Character vector with a name for each layer in
#'   \code{layers} (used in plots and summaries).
#' @param save_folder Name of the output folder (created inside
#'   \code{working_directory}) where plots and intermediate files are saved.
#' @param layer_top_var Numeric vector (length = number of layers) with the
#'   number of top variable genes/entities to consider per layer, or "all".
#' @param layer_min_corr Numeric vector with the minimum correlation used to
#'   start the cutoff search per layer.
#' @param layer_range_cutoff_length Integer vector with the number of
#'   cutoffs to test between \code{layer_min_corr} and the maximum observed
#'   correlation per layer.
#' @param print_distribution_plots Logical; whether to save distribution
#'   plots for all tested cutoffs.
#' @param cutoff_vec Optional numeric vector of final correlation cutoffs
#'   (one per layer). If \code{NULL}, automatically chosen via
#'   \code{choose_auto_cutoff()}.
#' @param sep_counts Separator for count files.
#' @param sep_anno Separator for annotation files.
#' @param gene_symbol_col Column name in count tables that contains gene
#'   symbols.
#' @param sample_col Column name in annotation tables that contains sample
#'   IDs.
#' @param count_has_rn Logical; does the count file already have row names?
#' @param anno_has_rn Logical; does the annotation file already have row
#'   names?
#'
#' @return A list with the main objects produced by the pipeline:
#'   \itemize{
#'     \item \code{global_settings}
#'     \item \code{layer_settings}
#'     \item \code{data}
#'     \item \code{supplementary_data}
#'     \item \code{layer_specific_outputs}
#'     \item \code{integrated_output}
#'   }
#'
#' @export
run_vcocena <- function(
  working_directory,
  layers,
  supplement,
  layers_names,
  save_folder = "out_all_1",
  layer_top_var,
  layer_min_corr,
  layer_range_cutoff_length,
  print_distribution_plots = FALSE,
  cutoff_vec = NULL,
  sep_counts = "\t",
  sep_anno = "\t",
  gene_symbol_col = "SYMBOL",
  sample_col = "ID",
  count_has_rn = FALSE,
  anno_has_rn = FALSE
) {
  if (length(layers) != length(layers_names)) {
    stop("`layers` and `layers_names` must have the same length.")
  }

  n_layers <- length(layers)

  if (length(layer_top_var) != n_layers ||
      length(layer_min_corr) != n_layers ||
      length(layer_range_cutoff_length) != n_layers) {
    stop("`layer_top_var`, `layer_min_corr`, and `layer_range_cutoff_length` must each have length equal to `length(layers)`.")
  }

  # Ensure trailing slash on working_directory
  if (!grepl("[/\\\\]$", working_directory)) {
    working_directory <- paste0(working_directory, "/")
  }

  # Assign shared objects into the global environment so that the existing
  # helper functions (which were written assuming a script/Rmd workflow)
  # can find them via lexical scoping. This is not ideal for a CRAN-style
  # package, but is pragmatic for this pipeline-style workflow.
  pkg_env <- .GlobalEnv

  assign("working_directory", working_directory, envir = pkg_env)
  assign("layers", layers, envir = pkg_env)
  assign("supplement", supplement, envir = pkg_env)
  assign("layers_names", layers_names, envir = pkg_env)

  # Global settings (minimal initialization here; you can enrich this list
  # as needed, e.g. by adding organism, voi, control, etc.).
  global_settings <- list()
  global_settings[["save_folder"]] <- init_save_folder(name = save_folder)
  global_settings[["voi"]] <- "Condition"
  global_settings[["control"]] <- "none"
  global_settings[["organism"]] <- "human"
  global_settings[["min_nodes_number_for_cluster"]] <- 3
  global_settings[["data_in_log"]] <- FALSE
  global_settings[["range_GFC"]] <- 3
  assign("global_settings", global_settings, envir = pkg_env)

  # Layer-specific settings
  layer_settings <- list()
  assign("layer_settings", layer_settings, envir = pkg_env)

  layer_settings <- set_layer_settings(
    top_var = layer_top_var,
    min_corr = layer_min_corr,
    range_cutoff_length = layer_range_cutoff_length,
    print_distribution_plots = rep(print_distribution_plots, n_layers)
  )
  assign("layer_settings", layer_settings, envir = pkg_env)

  # Data import and supplementary data
  data <- initialize_environment(
    sep_counts = sep_counts,
    sep_anno = sep_anno,
    gene_symbol_col = gene_symbol_col,
    sample_col = sample_col,
    count_has_rn = count_has_rn,
    anno_has_rn = anno_has_rn
  )
  supplementary_data <- load_supplementary_data()

  assign("data", data, envir = pkg_env)
  assign("supplementary_data", supplementary_data, envir = pkg_env)

  # Per-layer outputs (part1: correlations and cutoff stats)
  layer_specific_outputs <- vector("list", n_layers)
  names(layer_specific_outputs) <- paste0("set", seq_len(n_layers))
  assign("layer_specific_outputs", layer_specific_outputs, envir = pkg_env)

  for (x in seq_len(n_layers)) {
    layer_specific_outputs[[x]][["part1"]] <- run_expression_analysis_1(x)
  }
  assign("layer_specific_outputs", layer_specific_outputs, envir = pkg_env)

  # Choose final cutoffs per layer (either automatic or user-provided)
  if (is.null(cutoff_vec)) {
    cutoff_vec <- choose_auto_cutoff()
  } else if (length(cutoff_vec) != n_layers) {
    stop("`cutoff_vec` must be NULL or have one cutoff value per layer.")
  }
  assign("cutoff_vec", cutoff_vec, envir = pkg_env)

  # Per-layer outputs (part2: heatmaps and GFCs)
  for (x in seq_len(n_layers)) {
    layer_specific_outputs[[x]][["part2"]] <- run_expression_analysis_2(
      x = x,
      grouping_v = NULL,
      plot_HM = TRUE
    )
  }
  assign("layer_specific_outputs", layer_specific_outputs, envir = pkg_env)

  # Build and cluster networks per layer
  create_network_per_layer()
  cluster_networks()
  cluster_GFC_df_per_layer()

  # Vertical integration
  integrated_output <- list()
  assign("integrated_output", integrated_output, envir = pkg_env)

  integrated_output[["all_cluster_corrs"]] <- find_all_cluster_corrs()
  integrated_output[["identicals"]] <- ids(min_corr = 0.8)

  integrated_output[["cross_corr_cutoff_stats"]] <- find_cross_corr_cutoff(
    n = 50,
    threshold = 0.8
  )

  # Choose a cross-layer correlation cutoff based on the highest mixture score
  best_idx <- which.max(integrated_output[["cross_corr_cutoff_stats"]]$mix_score)
  cross_threshold <- as.numeric(
    integrated_output[["cross_corr_cutoff_stats"]]$cutoff[best_idx]
  )

  integrated_output[["cross_correlations"]] <- cross_corrs(
    threshold = cross_threshold
  )

  integrated_output[["integrated_edgelist"]] <- create_integrated_edgelist(
    el = integrated_output[["cross_correlations"]]
  )

  merged_net <- igraph::graph_from_data_frame(
    integrated_output[["integrated_edgelist"]],
    directed = FALSE
  )
  merged_net <- igraph::simplify(
    merged_net,
    remove.multiple = TRUE,
    remove.loops = TRUE
  )

  integrated_output[["merged_net"]] <- merged_net
  integrated_output[["cluster_calc"]] <- list()
  integrated_output$cluster_calc[["cluster_information"]] <- cluster_integrated_network(
    network = merged_net
  )
  assign("integrated_output", integrated_output, envir = pkg_env)

  # Derive integrated GFC per cluster (stored inside integrated_output)
  cluster_GFC_df_integrated()

  # Return a concise summary object
  list(
    global_settings = global_settings,
    layer_settings = layer_settings,
    data = data,
    supplementary_data = supplementary_data,
    layer_specific_outputs = layer_specific_outputs,
    integrated_output = integrated_output
  )
}
