# vcocena

[![R-CMD-check](https://github.com/BioCompNet/Vcocena/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/BioCompNet/Vcocena/actions/workflows/R-CMD-check.yaml)
[![Codecov](https://codecov.io/gh/BioCompNet/Vcocena/branch/main/graph/badge.svg)](https://codecov.io/gh/BioCompNet/Vcocena)
[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://biocompnet.github.io/Vcocena/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

`vcocena` provides a packaged, runnable version of the Vertical CoCena2 multi‑omics pipeline, plus plotting helpers.

## Install

From a local checkout:

```r
# install.packages("devtools")
devtools::install("vcocena")
```

## Quickstart

```r
library(vcocena)

ex <- vcocena_example()

layers <- list(
  set1 = c(ex$counts_seq, ex$anno_seq),
  set2 = c(ex$counts_array, ex$anno_array)
)

supplement <- c(
  ex$tf,
  ex$hallmark,
  ex$go,
  ex$kegg,
  ex$reactome
)

layers_names <- c("sequencing", "microarray")

out <- run_vcocena(
  working_directory = ex$working_directory,
  layers = layers,
  supplement = supplement,
  layers_names = layers_names,
  save_folder = "out_v1",
  layer_top_var = c(100, 100),
  layer_min_corr = c(0.2, 0.2),
  layer_range_cutoff_length = c(10, 10)
)
```

## Notes

This package depends on several CRAN and Bioconductor packages. If you are missing Bioconductor dependencies, install them with `BiocManager::install()`.
