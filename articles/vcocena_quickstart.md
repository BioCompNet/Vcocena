# vcocena Quickstart

## Overview

This vignette shows how to run a minimal vCoCena pipeline using the
bundled example dataset.

## Example data

``` r
library(vcocena)
```

    #> 

    #> Warning: replacing previous import 'biomaRt::select' by
    #> 'clusterProfiler::select' when loading 'vcocena'

    #> Registered S3 method overwritten by 'ggnetwork':
    #>   method         from  
    #>   fortify.igraph ggtree

    #> Warning: replacing previous import 'dplyr::as_data_frame' by
    #> 'igraph::as_data_frame' when loading 'vcocena'

    #> Warning: replacing previous import 'dplyr::groups' by 'igraph::groups' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'clusterProfiler::simplify' by
    #> 'igraph::simplify' when loading 'vcocena'

    #> Warning: replacing previous import 'dplyr::union' by 'igraph::union' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'circlize::degree' by 'igraph::degree' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'ComplexHeatmap::pheatmap' by
    #> 'pheatmap::pheatmap' when loading 'vcocena'

    #> Warning: replacing previous import 'ComplexHeatmap::add_heatmap' by
    #> 'plotly::add_heatmap' when loading 'vcocena'

    #> Warning: replacing previous import 'igraph::groups' by 'plotly::groups' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'BSgenome::export' by 'plotly::export' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'httr::config' by 'plotly::config' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'ggplot2::last_plot' by 'plotly::last_plot'
    #> when loading 'vcocena'

    #> Warning: replacing previous import 'magrittr::set_names' by 'purrr::set_names'
    #> when loading 'vcocena'

    #> Warning: replacing previous import 'igraph::simplify' by 'purrr::simplify' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'jsonlite::flatten' by 'purrr::flatten' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'igraph::compose' by 'purrr::compose' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'readr::col_factor' by 'scales::col_factor'
    #> when loading 'vcocena'

    #> Warning: replacing previous import 'purrr::discard' by 'scales::discard' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'igraph::crossing' by 'tidyr::crossing' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'magrittr::extract' by 'tidyr::extract' when
    #> loading 'vcocena'

    #> Warning: replacing previous import 'bnstruct::complete' by 'tidyr::complete'
    #> when loading 'vcocena'

``` r
ex <- vcocena_example()
```

## Define layers and supplements

``` r
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
```

## Run the pipeline (lightweight settings)

``` r
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

The full pipeline can take time on real data. Start with small values
and increase once the workflow is confirmed.
