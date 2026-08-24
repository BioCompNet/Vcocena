#!/usr/bin/env Rscript
# Tier 2: end-to-end run on the real two-layer dataset.
#
#   docker run --rm \
#     -v "$PWD/seq_micro_array_example_vertical:/data:ro" \
#     -v "$PWD/out_docker:/out" \
#     vcocena:1.0.0 Rscript /opt/docker/run_full.R
#
# 11,796 genes x 9 samples x 2 layers (RNA-seq + Illumina BeadArray),
# IFNg / IL4 / baseline, balanced 3/3/3, both layers already intersected to a
# common gene universe. Expect 10-30 min and >= 8 GB: correlation_actions()
# materializes a data frame of all choose(5000,2) ~ 12.5M gene pairs per layer.

suppressPackageStartupMessages(library(vcocena))

data_mount <- Sys.getenv("VCOCENA_DATA", "/data")
out_mount  <- Sys.getenv("VCOCENA_OUT",  "/out")
if (!dir.exists(data_mount))
  stop("dataset mount not found at ", data_mount,
       " -- bind-mount the dataset directory there")

# /data is mounted read-only, but the pipeline writes its save_folder inside
# working_directory. Copy the (small, ~2.3 MB) inputs to writable scratch.
wd <- "/work/dataset"
unlink(wd, recursive = TRUE)
dir.create(file.path(wd, "data"),        recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(wd, "sample_info"), recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(file.path(data_mount, "data"), full.names = TRUE),
          file.path(wd, "data"), recursive = TRUE)
file.copy(list.files(file.path(data_mount, "sample_info"), full.names = TRUE),
          file.path(wd, "sample_info"), recursive = TRUE)

save_folder <- "out_docker_full"
t0 <- Sys.time()

res <- run_vcocena(
  working_directory         = wd,
  layers                    = list(
    set1 = c("data_seq_processed.txt",   "annotation_seq.txt"),
    set2 = c("data_array_processed.txt", "annotation_array.txt")
  ),
  supplement                = c("TFcat.txt", "toy_hallmark.gmt", "toy_go.gmt",
                                "toy_kegg.gmt", "toy_reactome.gmt"),
  layers_names              = c("sequencing", "microarray"),
  save_folder               = save_folder,

  # Parameters matching the historical run recorded in the tree. Note
  # range_cutoff_length = 50, not the 500 in vCoCena_stepwise_part1.R:61 -- the
  # 52 distinct cutoffs in dir_DegreeDistribution_5000/ show the run that
  # produced out_stepwise/ used 50, and 500 is ~10x slower for no benefit.
  layer_top_var             = c(5000, 5000),
  layer_min_corr            = c(0.6, 0.6),
  layer_range_cutoff_length = c(50, 50),
  print_distribution_plots  = FALSE,

  # This dataset is rlog/log2 with batch effects removed. Leaving the old
  # hardcoded data_in_log = FALSE would produce wrong Group Fold Changes.
  data_in_log                  = TRUE,
  voi                          = "Condition",
  min_nodes_number_for_cluster = 50,
  range_GFC                    = 3
)

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

# run_vcocena() returns a partially stale object: create_network_per_layer(),
# cluster_networks() and cluster_GFC_df_per_layer() write to the GLOBAL
# layer_specific_outputs via <<-, but run_vcocena.R returns its local copy.
# Read the globals for anything network- or cluster-related.
lso <- get("layer_specific_outputs", envir = .GlobalEnv)
int <- get("integrated_output",      envir = .GlobalEnv)

cat("\n================ FULL RUN ================\n")
cat("elapsed:            ", elapsed, "min\n")
cat("cutoffs chosen:     ", paste(get("cutoff_vec", envir = .GlobalEnv), collapse = ", "), "\n")

for (i in seq_along(lso)) {
  net <- lso[[i]]$network
  ci  <- lso[[i]]$cluster_calc$cluster_information
  cat(sprintf("layer %d (%-10s): %6d nodes %8d edges %3d clusters\n",
              i, c("sequencing", "microarray")[i],
              igraph::gorder(net), igraph::gsize(net),
              if (is.null(ci)) 0L else nrow(ci)))
}

mn     <- int$merged_net
ci_int <- int$cluster_calc$cluster_information
kept   <- if (is.null(ci_int)) 0L else sum(ci_int$cluster_included == "yes")

cat(sprintf("integrated network:  %6d nodes %8d edges\n",
            igraph::gorder(mn), igraph::gsize(mn)))
cat("integrated modules: ", if (is.null(ci_int)) 0L else nrow(ci_int),
    "(", kept, "retained )\n")
cat("cross-layer edges:  ", nrow(int$cross_correlations), "\n")

# Publish results to the output mount.
outdir <- file.path(wd, save_folder)
files  <- list.files(outdir, recursive = TRUE)
if (dir.exists(out_mount)) {
  file.copy(list.files(outdir, full.names = TRUE), out_mount, recursive = TRUE)
  # Freeze the numeric results as TSV. The pipeline has no tabular export path
  # of its own, so this is what a future regression baseline compares against.
  if (!is.null(ci_int))
    write.table(ci_int, file.path(out_mount, "integrated_cluster_information.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(int$cluster_calc$GFC_per_cluster))
    write.table(int$cluster_calc$GFC_per_cluster,
                file.path(out_mount, "integrated_GFC_per_cluster.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(int$integrated_edgelist))
    write.table(int$integrated_edgelist,
                file.path(out_mount, "integrated_edgelist.tsv"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  cat("results copied to: ", out_mount, "\n")
}
cat("files written:      ", length(files), "\n")
for (f in files) cat("  -", f, "\n")

stopifnot(
  igraph::gorder(mn) > 0,
  igraph::gsize(mn)  > 0,
  !is.null(ci_int), nrow(ci_int) >= 1
)
cat("\nFULL RUN PASSED\n")
