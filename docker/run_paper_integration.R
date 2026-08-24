#!/usr/bin/env Rscript
# Integration-only reproduction test.
#
# The published vCoCena run consumed hCoCena-combined layer networks, not count
# matrices, so run_vcocena() cannot reproduce it end-to-end. This harness skips
# correlation/cutoff/per-layer clustering entirely and injects the PUBLISHED
# layer networks and per-layer modules straight into the vertical-integration
# stage -- the part vcocena actually owns -- then compares the integrated
# modules against the published ones.
suppressPackageStartupMessages(library(vcocena))

rd <- Sys.getenv("VCOCENA_RDATA", "/paper/vCoCena_3-CePs_GH.RData")
out <- Sys.getenv("VCOCENA_OUT", "/out")
e <- new.env(); cat("loading published environment...\n"); load(rd, envir = e)

# The pipeline's helpers resolve their inputs through .GlobalEnv, so seed it
# with the published state rather than anything we recomputed.
G <- .GlobalEnv
assign("layers",       e$layers,       envir = G)
assign("layers_names", e$layers_names, envir = G)
assign("data",         e$data,         envir = G)
assign("layer_settings", e$layer_settings, envir = G)
gs <- e$global_settings; gs$save_folder <- "."
assign("global_settings", gs, envir = G)
assign("working_directory", "/work/", envir = G)

# Published layer networks + published per-layer modules, verbatim.
lso <- e$layer_specific_outputs
assign("layer_specific_outputs", lso, envir = G)
for (s in names(lso))
  cat(sprintf("  %-5s network %d nodes / %d edges, %d modules\n", s,
              igraph::gorder(lso[[s]]$network), igraph::gsize(lso[[s]]$network),
              nrow(lso[[s]]$cluster_calc$cluster_information)))

io <- list(); assign("integrated_output", io, envir = G)
t0 <- Sys.time()

io[["all_cluster_corrs"]] <- vcocena:::find_all_cluster_corrs()
assign("integrated_output", io, envir = G)
cat("cluster-pair correlations:", nrow(io$all_cluster_corrs[[1]]), "pairs\n")

io[["identicals"]] <- vcocena:::ids(crosscorrlist = io$all_cluster_corrs, min_corr = 0.65)
assign("integrated_output", io, envir = G)
cat("identical cluster pairs (r > 0.65):", nrow(io$identicals[[1]]), "\n")

io[["cross_correlations"]] <- vcocena:::cross_corrs(threshold = 0.665)
assign("integrated_output", io, envir = G)
cat("cross-layer edges (r >= 0.665):", nrow(io$cross_correlations), "\n")

io[["integrated_edgelist"]] <- vcocena:::create_integrated_edgelist(el = io$cross_correlations)
assign("integrated_output", io, envir = G)

mn <- igraph::graph_from_data_frame(io$integrated_edgelist, directed = FALSE)
mn <- igraph::simplify(mn, remove.multiple = TRUE, remove.loops = TRUE)
io[["merged_net"]] <- mn
io[["cluster_calc"]] <- list()
assign("integrated_output", io, envir = G)

io$cluster_calc[["cluster_information"]] <- vcocena:::cluster_integrated_network(network = mn)
assign("integrated_output", io, envir = G)

mins <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
ci <- io$cluster_calc$cluster_information

cat("\n=========== INTEGRATION-ONLY REPRODUCTION ===========\n")
cat("elapsed:", mins, "min\n")
cat(sprintf("integrated: %d nodes, %d edges, %d modules (%d included)\n",
            igraph::gorder(mn), igraph::gsize(mn), nrow(ci),
            sum(ci$cluster_included == "yes")))
saveRDS(list(cluster_information = ci,
             GFC_per_cluster = NULL,
             n_nodes = igraph::gorder(mn), n_edges = igraph::gsize(mn),
             n_edgelist = nrow(io$integrated_edgelist),
             n_cross = nrow(io$cross_correlations)),
        file.path(out, "ours_integration_only.rds"))
cat("saved to", file.path(out, "ours_integration_only.rds"), "\n")
