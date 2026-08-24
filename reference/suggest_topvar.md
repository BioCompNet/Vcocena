# Suggest top variable features per layer

Computes a suggested number of top variable genes/features for each
layer.

## Usage

``` r
suggest_topvar(
  data = data,
  layers_names = layers_names,
  layers = layers,
  inflection_plot_file = NULL,
  summary_plot_file = NULL
)
```

## Arguments

- data:

  Data list with `set*_counts` entries.

- layers_names:

  Character vector of layer names.

- layers:

  List of layers (used for iteration).

- inflection_plot_file:

  Optional PDF path to save inflection plots.

- summary_plot_file:

  Optional PDF path to save a summary barplot.

## Value

A numeric vector of suggested top variable feature counts.
