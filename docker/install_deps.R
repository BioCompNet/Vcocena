#!/usr/bin/env Rscript
# Dependency installer for the vcocena image.
#
# Order matters: MCDA was removed from CRAN on 2025-11-20 but is still declared
# in DESCRIPTION and used on the main path (optimal_cutoff_MultiOmics.R:33,36).
# It must be present BEFORE install_deps() resolves dependencies, otherwise
# resolution fails on a package no repository can serve.

options(
  Ncpus       = max(1L, parallel::detectCores()),
  warn        = 1,
  timeout     = 1800,
  repos       = BiocManager::repositories()   # CRAN + all Bioconductor repos
)

message("== R ", getRversion(), " / Bioconductor ", BiocManager::version(), " ==")

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

install_tarball <- function(url, name) {
  if (requireNamespace(name, quietly = TRUE)) {
    message("-- ", name, " already present, skipping")
    return(invisible(TRUE))
  }
  message("-- installing ", name, " from ", url)
  remotes::install_url(url, upgrade = "never")
  if (!requireNamespace(name, quietly = TRUE))
    stop("failed to install ", name, " from ", url)
}

# 1a. glpkAPI: archived from CRAN on 2026-03-20, and an Imports of every MCDA
#     release. Needs the GLPK C library, which the apt layer installs.
install_tarball(
  "https://cran.r-project.org/src/contrib/Archive/glpkAPI/glpkAPI_1.3.4.1.tar.gz",
  "glpkAPI"
)

# 1b. MCDA: archived on CRAN 2025-11-20, still imported wholesale by the package
#     NAMESPACE and used unqualified in optimal_cutoff_MultiOmics.R:33,36.
#     Its remaining Imports (Rglpk, triangle, plyr, combinat, ggplot2,
#     RColorBrewer) are still live on CRAN and resolve normally.
install_tarball(
  "https://cran.r-project.org/src/contrib/Archive/MCDA/MCDA_0.1.0.tar.gz",
  "MCDA"
)

# 2. combinat: a 2012 package with an open CRAN check deadline. The tarball is
#    currently live; pinning it here means an archival does not break the build.
install_tarball(
  "https://cran.r-project.org/src/contrib/combinat_0.0-8.tar.gz",
  "combinat"
)

# 3. Everything else, resolved from DESCRIPTION alone.
#    Only /opt/deps/DESCRIPTION is present at this layer -- the package source is
#    copied later, so editing R/ does not invalidate this (slow) layer.
#    NOTE: R/install_bioconductor_packages.R is deliberately NOT sourced -- it
#    installs pcaGoPromoter{,.Hs.hg19,.Mm.mm9}, all removed from Bioconductor
#    as of release 3.12, and would abort the build.
message("-- installing declared dependencies from DESCRIPTION")
#    Hard dependencies only. Suggests (devtools, pkgdown, rcmdcheck, covr,
#    rmarkdown, testthat) are build-time tooling, not runtime requirements, and
#    pulling them roughly doubles the install for no benefit in this image.
remotes::install_deps(
  "/opt/deps",
  dependencies = c("Depends", "Imports", "LinkingTo"),
  upgrade      = "never"
)

# 4. org.Hs.eg.db is undeclared in DESCRIPTION but needed by the GO profiling
#    helpers (GO_profiling.R:13,67). Not on the run_vcocena() path, but cheap
#    insurance so the image can do enrichment without a rebuild.
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db", ask = FALSE, update = FALSE)

message("== dependency install complete ==")
