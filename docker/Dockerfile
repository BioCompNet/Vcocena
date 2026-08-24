# vcocena -- Vertical CoCena2 multi-omics integration pipeline
#
# Build from the repository root (this file lives in docker/):
#   docker build -f docker/Dockerfile -t vcocena:1.0.0 .
#
# The Bioconductor base carries R 4.5.2 plus a matched BiocManager and most of
# the system libraries the Bioconductor dependency set needs. R 4.5.x also
# matches the R 4.5.0 that produced the existing outputs in out_stepwise/.
FROM bioconductor/bioconductor_docker:RELEASE_3_22-R-4.5.2

LABEL org.opencontainers.image.title="vcocena" \
      org.opencontainers.image.description="Vertical CoCena2 multi-omics integration and clustering" \
      org.opencontainers.image.source="https://github.com/BioCompNet/Vcocena" \
      org.opencontainers.image.licenses="MIT"

# System libraries the Bioconductor base does not already provide.
#   libcairo2-dev / libxt-dev -> the Cairo package AND R's own cairo_pdf device,
#     which heatmap_top_var_script.R:87 uses on the main execution path
#   libglpk-dev               -> igraph community detection
#   libgmp-dev / libmpfr-dev  -> transitive (bnstruct, qgraph)
RUN apt-get update && apt-get install -y --no-install-recommends \
        libcairo2-dev \
        libxt-dev \
        libglpk-dev \
        libgmp-dev \
        libmpfr-dev \
    && rm -rf /var/lib/apt/lists/*

# --- Layer 1: dependencies (slow, ~30 min). Only DESCRIPTION is copied here, so
# --- editing package source below does not invalidate this layer.
COPY DESCRIPTION /opt/deps/DESCRIPTION
COPY docker/install_deps.R /opt/docker/install_deps.R
RUN Rscript /opt/docker/install_deps.R

# --- Layer 2: the package itself (fast). Datasets are NOT baked in -- they are
# --- bind-mounted at runtime.
COPY . /opt/Vcocena/
COPY docker/ /opt/docker/
RUN Rscript -e 'remotes::install_local("/opt/Vcocena", dependencies = FALSE, upgrade = "never")'

# Build-time assertions: fail here rather than three steps into a pipeline run.
RUN Rscript -e ' \
      stopifnot(capabilities("cairo")); \
      req <- c("Hmisc","gtools","tibble","grid","knitr","MCDA","combinat", \
               "lineup","clusterProfiler","igraph","pheatmap","ggplot2"); \
      missing <- req[!vapply(req, requireNamespace, logical(1), quietly = TRUE)]; \
      if (length(missing)) stop("missing packages: ", paste(missing, collapse=", ")); \
      library(vcocena); \
      stopifnot(is.function(vcocena::run_vcocena)); \
      cat("build assertions passed; vcocena", \
          as.character(packageVersion("vcocena")), "\n") \
    '

# Writable scratch. run_vcocena() calls dir.create() into working_directory and
# ggsave()s into it, so the analysis directory can never be a read-only mount
# nor the installed package library.
RUN mkdir -p /work /out
ENV HOME=/work
WORKDIR /work

CMD ["Rscript", "/opt/docker/run_smoke.R"]
