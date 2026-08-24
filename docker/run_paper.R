#!/usr/bin/env Rscript
# Reproduction run against Carraro et al., eLife 2022 (doi:10.7554/eLife.78012),
# the paper that introduced vCoCena.
#
# Inputs and ground truth both come from the published R environment
# `vCoCena_3-CePs_GH` (Zenodo 10.5281/zenodo.6984701, CC-BY-4.0), so this
# compares our package against the authors' own recorded output rather than
# against a re-derivation of it.
#
# Every parameter below is read from the published `global_settings` /
# `layer_settings`, not from the manuscript text.

suppressPackageStartupMessages(library(vcocena))

data_mount <- Sys.getenv("VCOCENA_DATA", "/data")
out_mount  <- Sys.getenv("VCOCENA_OUT",  "/out")

wd <- "/work/paper"
unlink(wd, recursive = TRUE)
dir.create(file.path(wd, "data"),        recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(wd, "sample_info"), recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(file.path(data_mount, "data"), full.names = TRUE),
          file.path(wd, "data"), recursive = TRUE)
file.copy(list.files(file.path(data_mount, "sample_info"), full.names = TRUE),
          file.path(wd, "sample_info"), recursive = TRUE)

t0 <- Sys.time()
res <- run_vcocena(
  working_directory = wd,
  layers = list(set1 = c("counts_RNA.txt",  "anno_RNA.txt"),
                set2 = c("counts_ATAC.txt", "anno_ATAC.txt")),
  supplement   = c("TFcat.txt", "toy_hallmark.gmt", "toy_go.gmt",
                   "toy_kegg.gmt", "toy_reactome.gmt"),
  layers_names = c("RNA", "ATAC"),
  save_folder  = "out_paper",

  # The published run did not use automatic cutoff selection -- `cutoff_vec` is
  # NULL in its saved environment, and the retained edge weights bottom out at
  # 0.707 (RNA) / 0.693 (ATAC), i.e. cutoffs chosen by hand from the diagnostic
  # plots, which is the normal CoCena workflow. Pin them so this run compares
  # the integration and clustering rather than the cutoff heuristic.
  cutoff_vec = as.numeric(strsplit(Sys.getenv("VCOCENA_CUTOFFS", "0.707,0.693"), ",")[[1]]),

  # --- published layer_settings ---
  layer_top_var             = c(1794, 1624),
  layer_min_corr            = c(0.5, 0.5),
  layer_range_cutoff_length = c(100, 100),

  # --- published global_settings ---
  voi                          = "Group",
  control                      = "none",
  organism                     = "human",
  min_nodes_number_for_cluster = 15,
  min_nodes_number_for_network = 15,
  range_GFC                    = 2,
  data_in_log                  = TRUE,

  # --- published integration thresholds (from vCoCena_3-CePs_GH.Rmd) ---
  identicals_min_corr  = 0.65,
  cross_corr_threshold = 0.665
)
mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)

lso <- get("layer_specific_outputs", envir = .GlobalEnv)
int <- get("integrated_output",      envir = .GlobalEnv)
ci  <- int$cluster_calc$cluster_information

cat("\n================ PAPER REPRODUCTION ================\n")
cat("elapsed:", mins, "min\n")
cat("cutoffs chosen:", paste(get("cutoff_vec", envir = .GlobalEnv), collapse = ", "), "\n")
for (i in seq_along(lso))
  cat(sprintf("layer %d (%-4s): %5d nodes %7d edges %3d clusters\n", i,
              c("RNA","ATAC")[i], igraph::gorder(lso[[i]]$network),
              igraph::gsize(lso[[i]]$network),
              nrow(lso[[i]]$cluster_calc$cluster_information)))
cat(sprintf("integrated: %d nodes, %d edges, %d modules (%d included)\n",
            igraph::gorder(int$merged_net), igraph::gsize(int$merged_net),
            nrow(ci), sum(ci$cluster_included == "yes")))
cat("cross-layer edges:", nrow(int$cross_correlations), "\n")

if (dir.exists(out_mount)) {
  saveRDS(list(cluster_information = ci,
               GFC_per_cluster = int$cluster_calc$GFC_per_cluster,
               n_nodes = igraph::gorder(int$merged_net),
               n_edges = igraph::gsize(int$merged_net),
               n_edgelist = nrow(int$integrated_edgelist),
               n_cross = nrow(int$cross_correlations),
               cutoffs = get("cutoff_vec", envir = .GlobalEnv)),
          file.path(out_mount, "ours_paper_run.rds"))
  write.table(ci, file.path(out_mount, "paper_cluster_information.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  file.copy(list.files(file.path(wd, "out_paper"), full.names = TRUE),
            out_mount, recursive = TRUE)
  cat("results written to", out_mount, "\n")
}
cat("PAPER RUN COMPLETE\n")
