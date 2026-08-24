# vcocena

[![R-CMD-check](https://github.com/BioCompNet/Vcocena/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BioCompNet/Vcocena/actions/workflows/R-CMD-check.yaml)
[![Codecov](https://codecov.io/gh/BioCompNet/Vcocena/branch/main/graph/badge.svg)](https://codecov.io/gh/BioCompNet/Vcocena)
[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://biocompnet.github.io/Vcocena/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`vcocena` provides a packaged, runnable version of the Vertical CoCena² multi‑omics
pipeline, plus plotting helpers.

Where *horizontal* integration ([hCoCena](https://github.com/MarieOestreich/hCoCena))
combines several datasets measured on one platform, **vertical** integration combines
different omics layers — RNA‑seq and microarray, or transcriptome and chromatin
accessibility. The pipeline builds a co‑expression network per layer, clusters each,
correlates clusters across layers, and merges everything into one integrated network
of modules.

## Install

Two dependencies are **no longer on CRAN** — `MCDA` (archived 2025-11-20) and its
dependency `glpkAPI` (archived 2026-03-20). Both sit on the main execution path, so
they cannot be skipped. Add a dated repository snapshot from before the archivals and
they resolve as ordinary packages:

```r
options(repos = c(
  CRAN     = "https://cloud.r-project.org",
  snapshot = "https://packagemanager.posit.co/cran/2025-11-01"
))

# install.packages("remotes")
remotes::install_github("BioCompNet/Vcocena")
```

The snapshot is a *secondary* repository: R takes the highest available version of
each package, so everything else still comes from current CRAN and only the two
archived packages come from the snapshot.

`glpkAPI` compiles against the GLPK C library — `libglpk-dev` on Debian/Ubuntu,
`brew install glpk` on macOS. The remaining Bioconductor dependencies install normally
via `BiocManager::install()`.

Do **not** call `install_bioconductor_packages()` from this package — it still lists
`pcaGoPromoter`, which was removed from Bioconductor in release 3.12.

### Docker (recommended)

The container pins every dependency, including the two archived ones, and needs no
manual setup:

```bash
git clone https://github.com/BioCompNet/Vcocena.git && cd Vcocena
docker build -f docker/Dockerfile -t vcocena:1.0.0 .
docker run --rm vcocena:1.0.0 Rscript /opt/docker/run_smoke.R
```

See [`docker/README.md`](docker/README.md) for running your own data.

## Quickstart

```r
library(vcocena)

ex <- vcocena_example()

# vcocena_example() points into the installed package library, which is usually
# read-only. run_vcocena() creates its save_folder inside working_directory and
# writes plots there, so copy the example somewhere writable first.
wd <- file.path(tempdir(), "example_vertical")
dir.create(wd, recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(ex$working_directory, full.names = TRUE), wd, recursive = TRUE)

out <- run_vcocena(
  working_directory = wd,
  layers = list(
    set1 = c(ex$counts_seq,   ex$anno_seq),
    set2 = c(ex$counts_array, ex$anno_array)
  ),
  supplement   = c(ex$tf, ex$hallmark, ex$go, ex$kegg, ex$reactome),
  layers_names = c("sequencing", "microarray"),
  save_folder  = "out_v1",
  layer_top_var             = c(100, 100),
  layer_min_corr            = c(0.2, 0.2),
  layer_range_cutoff_length = c(10, 10),
  data_in_log               = TRUE
)
```

The bundled example is a 100‑gene subset intended only as a smoke test; its variance
range is too narrow for meaningful modules. Use it to check the install, not to judge
results.

### Reading the results

`run_vcocena()` returns a list, but the network and cluster objects are written to the
global environment rather than into that return value. Read them from there:

```r
lso <- get("layer_specific_outputs", envir = .GlobalEnv)  # per-layer networks, clusters
int <- get("integrated_output",      envir = .GlobalEnv)  # integrated network + modules

int$cluster_calc$cluster_information   # module assignments
int$cluster_calc$GFC_per_cluster       # Group Fold Change per module per condition
```

There is no tabular export; results live as in‑memory objects and PDFs. `docker/run_full.R`
shows how to write them out as TSV.

## Important parameters

| Argument | Why it matters |
| --- | --- |
| `data_in_log` | Set `TRUE` for rlog/VST/log2 input. Leaving it `FALSE` silently produces wrong Group Fold Changes. |
| `layer_top_var` | Genes/features per layer. `suggest_topvar()` proposes values from variance inflection points. |
| `cutoff_vec` | Correlation cutoff per layer. If `NULL` these are chosen automatically — see the caveat below. |
| `voi` | Annotation column defining groups for GFC calculation. |
| `identicals_min_corr`, `cross_corr_threshold` | Control which cross‑layer cluster pairs are linked. |

## Reproducibility notes

**Automatic cutoff selection is stringent.** On a published dataset, `choose_auto_cutoff()`
selected 0.984/0.991 where the original authors had chosen 0.707/0.693 by hand from the
diagnostic plots — retaining 878 entities instead of 3,405. Inspect `cutoff_plot()` and
consider passing `cutoff_vec` explicitly rather than relying on the automatic choice.

**Module counts will not match analyses run before ~2024.** igraph 2.0 rewrote community
detection. On an identical network, `cluster_louvain` returns 6 communities where igraph
1.2.6 returned 147. Re‑running an older CoCena analysis on current igraph yields coarser
modules; the entities do not scatter, they merge.

**`run_vcocena()` takes count matrices, one per layer.** If your layers are themselves
hCoCena integrations across several datasets, this entry point does not accept them —
that was the design in Carraro et al. 2022, and it needs the hCoCena stage first.

`docker/run_paper.R` and `docker/run_paper_integration.R` check this package against
[Carraro et al., eLife 2022](https://doi.org/10.7554/eLife.78012) using the authors'
own saved R environment as ground truth. The vertical‑integration stage reproduces the
published network exactly: 3,405 nodes, 186,723 edges, 35,453 cross‑layer edges.

## Citation

See [`CITATION.cff`](CITATION.cff). The method is described in
[Carraro et al., eLife 2022](https://doi.org/10.7554/eLife.78012).
