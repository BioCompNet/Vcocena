# Run the Vertical CoCena² pipeline

This function mirrors the main steps of `vCoCena.Rmd` and orchestrates
data loading, per-layer correlation / cutoff analysis, GFC calculation,
network construction, clustering, and vertical integration across
layers.

## Usage

``` r
run_vcocena(
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
)
```

## Arguments

- working_directory:

  Path to the project directory that contains `data/`, `sample_info/`,
  and `reference_files/`.

- layers:

  Named list describing omics layers. Each entry should be a character
  vector of length 2: `c(count_file, annotation_file)`.

- supplement:

  Character vector of length 5 with file names (or paths relative to
  `reference_files/`) for TF list, Hallmark, GO, KEGG, and Reactome GMT
  files.

- layers_names:

  Character vector with a name for each layer in `layers` (used in plots
  and summaries).

- save_folder:

  Name of the output folder (created inside `working_directory`) where
  plots and intermediate files are saved.

- layer_top_var:

  Numeric vector (length = number of layers) with the number of top
  variable genes/entities to consider per layer, or "all".

- layer_min_corr:

  Numeric vector with the minimum correlation used to start the cutoff
  search per layer.

- layer_range_cutoff_length:

  Integer vector with the number of cutoffs to test between
  `layer_min_corr` and the maximum observed correlation per layer.

- print_distribution_plots:

  Logical; whether to save distribution plots for all tested cutoffs.

- cutoff_vec:

  Optional numeric vector of final correlation cutoffs (one per layer).
  If `NULL`, automatically chosen via `choose_auto_cutoff()`.

- sep_counts:

  Separator for count files.

- sep_anno:

  Separator for annotation files.

- gene_symbol_col:

  Column name in count tables that contains gene symbols.

- sample_col:

  Column name in annotation tables that contains sample IDs.

- count_has_rn:

  Logical; does the count file already have row names?

- anno_has_rn:

  Logical; does the annotation file already have row names?

- voi:

  Character; the "variable of interest" column in the annotation tables
  used to define groups for Group Fold Change calculation.

- control:

  Character; name of the control group, or "none".

- organism:

  Character; organism name used by downstream enrichment.

- min_nodes_number_for_cluster:

  Integer; minimum cluster size to retain.

- min_nodes_number_for_network:

  Integer; minimum connected-component size counted as a "network" when
  scoring correlation cutoffs. Previously unset, which made the network
  count always zero and silently degraded automatic cutoff selection; 40
  matches the original pipeline default.

- data_in_log:

  Logical; is the input already log-transformed? Set TRUE for
  rlog/VST/log2 data, otherwise Group Fold Changes will be wrong.

- range_GFC:

  Numeric; Group Fold Changes are clamped to +/- this value.

- identicals_min_corr:

  Numeric; minimum cluster-to-cluster correlation for two clusters in
  different layers to be treated as "identical".

- cross_corr_threshold:

  Numeric or NULL. Cutoff applied to cross-layer entity correlations
  when building the integrated edge list. NULL (default) selects it
  automatically as the cutoff with the highest mixture score; supply a
  number to pin it, e.g. to reproduce a published analysis.

- cross_corr_scan_n:

  Integer; number of cutoffs scanned by `find_cross_corr_cutoff()`.

- cross_corr_scan_min:

  Numeric; lower bound of that scan.

## Value

A list with the main objects produced by the pipeline:

- `global_settings`

- `layer_settings`

- `data`

- `supplementary_data`

- `layer_specific_outputs`

- `integrated_output`
