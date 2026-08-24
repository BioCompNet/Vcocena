# Circular module heatmap

Draw a circular heatmap for one or more modules, optionally splitting
tracks by omics layer and aggregating samples into condition means.

## Usage

``` r
plot_circular_module_heatmap(
  modules = NULL,
  cluster_info = NULL,
  layer_outputs = NULL,
  layers_names = NULL,
  layer_indices = NULL,
  data_list = NULL,
  condition_col = NULL,
  condition_order = NULL,
  condition_mode = c("intersect", "union"),
  condition_clean_regex = NULL,
  max_genes = 200,
  order_by = c("mean_abs", "none"),
  scale_rows = TRUE,
  show_gene_labels = FALSE,
  label_cex = 0.4,
  show_legend = TRUE,
  legend_title = "Scaled expr",
  legend_at = NULL,
  legend_labels = NULL,
  legend_x = NULL,
  legend_y = NULL,
  legend_just = c("right", "center"),
  show_condition_labels = TRUE,
  condition_label_prefix = "Track (inner→outer): ",
  condition_label_text = NULL,
  condition_label_gp = NULL,
  output_file = NULL,
  width_in = 8,
  height_in = 8,
  dpi = 300,
  na_col = "grey90",
  low_col = "navy",
  mid_col = "white",
  high_col = "firebrick",
  zlim = NULL,
  start_degree = 90,
  gap_degree = 10,
  track_height = 0.3,
  cell_padding = c(0, 0, 0, 0)
)
```

## Arguments

- modules:

  Module identifiers (e.g., `"M1"`) or module colors.

- cluster_info:

  Cluster information table.

- layer_outputs:

  Layer-specific outputs (from the pipeline).

- layers_names:

  Character vector of layer names.

- layer_indices:

  Integer vector of layer indices to include.

- data_list:

  Data list containing `set*_counts` and `set*_anno`.

- condition_col:

  Column name in annotations to define conditions.

- condition_order:

  Optional vector to order condition columns.

- condition_mode:

  How to reconcile condition columns across layers.

- condition_clean_regex:

  Regex to clean condition names (e.g. suffix removal).

- max_genes:

  Maximum number of genes to plot.

- order_by:

  Gene ordering method.

- scale_rows:

  Logical; scale per-gene values.

- show_gene_labels:

  Logical; show gene labels around the circle.

- output_file:

  Optional output PDF path.

## Value

Invisibly returns `NULL`; draws plot to device.
