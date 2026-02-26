rank_variance <- function(expr) {
  var_df <- base::data.frame(variance = base::apply(expr, 1, stats::var))
  var_df[["gene"]] <- base::rownames(expr)
  var_df <- var_df[base::order(var_df[["variance"]], decreasing = TRUE), ]
  var_df[["rank"]] <- base::seq_len(base::nrow(var_df))
  var_df
}

find_span <- function(var_df) {
  for (i in base::rev(base::seq(0.5, 1, 0.01))) {
    lo <- stats::loess(base::log(var_df$variance) ~ var_df$rank, span = i)
    xl <- base::seq(base::min(var_df$rank), base::max(var_df$rank),
                    (base::max(var_df$rank) - base::min(var_df$rank)) / 1000)
    out <- stats::predict(lo, xl)
    infl <- base::as.logical(base::diff(base::sign(base::diff(base::diff(out, differences = 1)))))
    if (base::length(infl[infl == TRUE]) > 0) {
      break
    }
  }
  i
}

plot_inflections <- function(var_df, i, setname = "my data") {
  lo <- stats::loess(base::log(var_df$variance) ~ var_df$rank, span = i)
  xl <- base::seq(base::min(var_df$rank), base::max(var_df$rank),
                  (base::max(var_df$rank) - base::min(var_df$rank)) / 1000)
  out <- stats::predict(lo, xl)
  infl <- base::as.logical(base::diff(base::sign(base::diff(base::diff(out, differences = 1)))))
  graphics::plot(var_df$rank, base::log(var_df$variance), type = "l",
                 xlab = "Rank", ylab = "log(Variance)", main = setname)
  graphics::points(xl[infl], out[infl], col = "red")
  base::cat(setname, ": Inflection points at the following #genes:\n")
  base::print(base::ceiling(xl[infl]))
  base::ceiling(xl[infl])
}

#' Suggest top variable features per layer
#'
#' Computes a suggested number of top variable genes/features for each layer.
#'
#' @param data Data list with \code{set*_counts} entries.
#' @param layers_names Character vector of layer names.
#' @param layers List of layers (used for iteration).
#' @param inflection_plot_file Optional PDF path to save inflection plots.
#' @param summary_plot_file Optional PDF path to save a summary barplot.
#' @return A numeric vector of suggested top variable feature counts.
#' @export
suggest_topvar <- function(data = data,
                           layers_names = layers_names,
                           layers = layers,
                           inflection_plot_file = NULL,
                           summary_plot_file = NULL) {
  if (!base::is.null(inflection_plot_file)) {
    grDevices::pdf(inflection_plot_file, width = 8, height = 6)
  }
  suggested <- base::numeric(base::length(layers_names))
  for (l in base::seq_along(layers)) {
    var_df <- rank_variance(data[[base::paste0("set", l, "_counts")]])
    i <- find_span(var_df = var_df)
    inflections <- plot_inflections(var_df = var_df, i = i, setname = layers_names[[l]])
    if (base::length(inflections) == 0) {
      suggested[l] <- NA_real_
    } else {
      suggested[l] <- inflections[1]
    }
  }
  if (!base::is.null(inflection_plot_file)) {
    grDevices::dev.off()
  }
  if (!base::is.null(summary_plot_file)) {
    grDevices::pdf(summary_plot_file, width = 8, height = 4)
    graphics::barplot(suggested,
                      names.arg = layers_names,
                      las = 2,
                      ylab = "Suggested top_var (# genes)",
                      main = "Suggested top_var per layer")
    grDevices::dev.off()
  }
  suggested
}
