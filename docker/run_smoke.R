#!/usr/bin/env Rscript
# Tier 1: fast self-contained smoke test on the packaged 100-gene example.
#
# This fixture is `head -101` of the real matrices, so its row variances top out
# around 0.02/0.16 versus 9.7/11.0 in the full data. It is statistically
# degenerate by construction. Assert only that the pipeline COMPLETES and emits
# files -- never on cluster identity or membership.

suppressPackageStartupMessages(library(vcocena))

ex <- vcocena_example()

# vcocena_example() points at the installed package library, which is read-only.
# init_save_folder() (init_envo.R:5-15) calls dir.create() into
# working_directory and the pipeline ggsave()s into it, so copy to /work first.
wd <- "/work/example_vertical"
unlink(wd, recursive = TRUE)
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
file.copy(
  list.files(ex$working_directory, full.names = TRUE),
  wd, recursive = TRUE
)

save_folder <- "out_smoke"
t0 <- Sys.time()

res <- run_vcocena(
  working_directory         = wd,
  layers                    = list(
    set1 = c(ex$counts_seq,   ex$anno_seq),
    set2 = c(ex$counts_array, ex$anno_array)
  ),
  supplement                = c(ex$tf, ex$hallmark, ex$go, ex$kegg, ex$reactome),
  layers_names              = c("sequencing", "microarray"),
  save_folder               = save_folder,
  layer_top_var             = c(100, 100),
  layer_min_corr            = c(0.2, 0.2),
  layer_range_cutoff_length = c(10, 10),
  data_in_log               = TRUE
)

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
outdir  <- file.path(wd, save_folder)
files   <- list.files(outdir, recursive = TRUE)

cat("\n================ SMOKE TEST ================\n")
cat("elapsed:      ", elapsed, "s\n")
cat("output dir:   ", outdir, "\n")
cat("files written:", length(files), "\n")
for (f in files) cat("  -", f, "\n")

stopifnot(is.list(res), length(files) > 0)
cat("\nSMOKE TEST PASSED\n")
