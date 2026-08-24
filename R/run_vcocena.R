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
#' @param voi Character; the "variable of interest" column in the annotation
#'   tables used to define groups for Group Fold Change calculation.
#' @param control Character; name of the control group, or "none".
#' @param organism Character; organism name used by downstream enrichment.
#' @param min_nodes_number_for_cluster Integer; minimum cluster size to retain.
#' @param min_nodes_number_for_network Integer; minimum connected-component
#'   size counted as a "network" when scoring correlation cutoffs. Previously
#'   unset, which made the network count always zero and silently degraded
#'   automatic cutoff selection; 40 matches the original pipeline default.
#' @param data_in_log Logical; is the input already log-transformed? Set TRUE
#'   for rlog/VST/log2 data, otherwise Group Fold Changes will be wrong.
#' @param range_GFC Numeric; Group Fold Changes are clamped to +/- this value.
#' @param identicals_min_corr Numeric; minimum cluster-to-cluster correlation
#'   for two clusters in different layers to be treated as "identical".
#' @param cross_corr_threshold Numeric or NULL. Cutoff applied to cross-layer
#'   entity correlations when building the integrated edge list. NULL (default)
#'   selects it automatically as the cutoff with the highest mixture score;
#'   supply a number to pin it, e.g. to reproduce a published analysis.
#' @param cross_corr_scan_n Integer; number of cutoffs scanned by
#'   \code{find_cross_corr_cutoff()}.
#' @param cross_corr_scan_min Numeric; lower bound of that scan.
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
  anno_has_rn = FALSE,
  voi = "Condition",
  control = "none",
  organism = "human",
  min_nodes_number_for_cluster = 3,
  min_nodes_number_for_network = 40,
  data_in_log = FALSE,
  range_GFC = 3,
  identicals_min_corr = 0.8,
  cross_corr_threshold = NULL,
  cross_corr_scan_n = 50,
  cross_corr_scan_min = 0.8
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

  # Global settings, driven by the corresponding function arguments.
  global_settings <- list()
  global_settings[["save_folder"]] <- init_save_folder(name = save_folder)
  global_settings[["voi"]] <- voi
  global_settings[["control"]] <- control
  global_settings[["organism"]] <- organism
  global_settings[["min_nodes_number_for_cluster"]] <- min_nodes_number_for_cluster
  # Read by cutoff_prep() via `global_set$min_nodes_number_for_network`; if it
  # is absent the component count collapses to 0 for every cutoff tested.
  global_settings[["min_nodes_number_for_network"]] <- min_nodes_number_for_network
  global_settings[["data_in_log"]] <- data_in_log
  global_settings[["range_GFC"]] <- range_GFC
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
  assign("integrated_output", integrated_output, envir = pkg_env)

  integrated_output[["identicals"]] <- ids(
    crosscorrlist = integrated_output[["all_cluster_corrs"]],
    min_corr = identicals_min_corr
  )
  assign("integrated_output", integrated_output, envir = pkg_env)

  if (length(integrated_output[["identicals"]]) == 0) {
    warning("No cross-layer cluster pairs above min_corr; proceeding without cross-layer edges.")
    integrated_output[["cross_corr_cutoff_stats"]] <- data.frame(
      mix_score = numeric(),
      cutoff = character(),
      num_clusters = integer(),
      num_cross_edges = integer(),
      stringsAsFactors = FALSE
    )
    integrated_output[["cross_correlations"]] <- data.frame(
      V1 = character(),
      V2 = character(),
      weight = numeric(),
      stringsAsFactors = FALSE
    )
    integrated_output[["integrated_edgelist"]] <- get_layer_edgelists()
  } else {
    integrated_output[["cross_corr_cutoff_stats"]] <- find_cross_corr_cutoff(
      n = cross_corr_scan_n,
      threshold = cross_corr_scan_min
    )
    assign("integrated_output", integrated_output, envir = pkg_env)

    if (nrow(integrated_output[["cross_corr_cutoff_stats"]]) == 0) {
      warning("Cross-correlation cutoff stats are empty; proceeding without cross-layer edges.")
      integrated_output[["cross_correlations"]] <- data.frame(
        V1 = character(),
        V2 = character(),
        weight = numeric(),
        stringsAsFactors = FALSE
      )
      integrated_output[["integrated_edgelist"]] <- get_layer_edgelists()
    } else {
      # Use the pinned cutoff when supplied, otherwise the highest mixture score.
      if (is.null(cross_corr_threshold)) {
        best_idx <- which.max(integrated_output[["cross_corr_cutoff_stats"]]$mix_score)
        cross_threshold <- as.numeric(
          integrated_output[["cross_corr_cutoff_stats"]]$cutoff[best_idx]
        )
      } else {
        cross_threshold <- as.numeric(cross_corr_threshold)
      }

      integrated_output[["cross_correlations"]] <- cross_corrs(
        threshold = cross_threshold
      )

      integrated_output[["integrated_edgelist"]] <- create_integrated_edgelist(
        el = integrated_output[["cross_correlations"]]
      )
    }
  }
  assign("integrated_output", integrated_output, envir = pkg_env)

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
