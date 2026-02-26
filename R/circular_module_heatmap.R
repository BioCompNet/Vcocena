#' Circular module heatmap
#'
#' Draw a circular heatmap for one or more modules, optionally splitting tracks
#' by omics layer and aggregating samples into condition means.
#'
#' @param modules Module identifiers (e.g., \code{"M1"}) or module colors.
#' @param cluster_info Cluster information table.
#' @param layer_outputs Layer-specific outputs (from the pipeline).
#' @param layers_names Character vector of layer names.
#' @param layer_indices Integer vector of layer indices to include.
#' @param data_list Data list containing \code{set*_counts} and \code{set*_anno}.
#' @param condition_col Column name in annotations to define conditions.
#' @param condition_order Optional vector to order condition columns.
#' @param condition_mode How to reconcile condition columns across layers.
#' @param condition_clean_regex Regex to clean condition names (e.g. suffix removal).
#' @param max_genes Maximum number of genes to plot.
#' @param order_by Gene ordering method.
#' @param scale_rows Logical; scale per-gene values.
#' @param show_gene_labels Logical; show gene labels around the circle.
#' @param output_file Optional output PDF path.
#' @return Invisibly returns \code{NULL}; draws plot to device.
#' @export
plot_circular_module_heatmap <- function(
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
) {
  order_by <- match.arg(order_by)
  condition_mode <- match.arg(condition_mode)
  if (is.null(cluster_info)) {
    if (exists("integrated_output", envir = .GlobalEnv)) {
      cluster_info <- get("integrated_output", envir = .GlobalEnv)$cluster_calc$cluster_information
    } else {
      stop("cluster_info is NULL and integrated_output not found.")
    }
  }
  if (is.null(layer_outputs)) {
    if (exists("layer_specific_outputs", envir = .GlobalEnv)) {
      layer_outputs <- get("layer_specific_outputs", envir = .GlobalEnv)
    } else {
      stop("layer_outputs is NULL and layer_specific_outputs not found.")
    }
  }
  if (is.null(layers_names)) {
    if (exists("layers_names", envir = .GlobalEnv)) {
      layers_names <- get("layers_names", envir = .GlobalEnv)
    } else {
      stop("layers_names is NULL and layers_names not found.")
    }
  }
  if (is.null(layer_indices)) {
    layer_indices <- seq_along(layer_outputs)
  }
  if (is.null(data_list)) {
    if (exists("data", envir = .GlobalEnv)) {
      data_list <- get("data", envir = .GlobalEnv)
    }
  }

  included_colors <- unique(cluster_info$color[cluster_info$cluster_included == "yes"])

  resolve_module_colors <- function(mods) {
    if (is.null(mods)) {
      return(included_colors)
    }
    out <- character(0)
    for (m in mods) {
      if (is.character(m) && grepl("^M\\d+$", m)) {
        idx <- as.integer(sub("^M", "", m))
        if (!is.na(idx) && idx >= 1 && idx <= length(included_colors)) {
          out <- c(out, included_colors[idx])
        }
      } else {
        out <- c(out, m)
      }
    }
    unique(out)
  }

  pick_condition_col <- function(anno_df) {
    if (is.null(anno_df) || ncol(anno_df) == 0) return(NULL)
    candidates <- setdiff(colnames(anno_df), c("Sample", "ID"))
    if (length(candidates) == 0) return(NULL)
    preferred <- intersect(c("Group", "group", "Condition", "condition", "merged", "Merged",
                             "treatment", "Treatment", "time", "Time", "TIME"), candidates)
    if (length(preferred) > 0) return(preferred[1])
    uniq_counts <- sapply(candidates, function(nm) length(unique(anno_df[[nm]])))
    uniq_counts <- uniq_counts[uniq_counts > 1]
    if (length(uniq_counts) == 0) return(NULL)
    names(sort(uniq_counts))[1]
  }

  prep_counts_matrix <- function(counts_df) {
    if (is.null(counts_df)) return(NULL)
    counts_df <- as.data.frame(counts_df)
    if (is.null(rownames(counts_df)) || any(rownames(counts_df) == "")) {
      if ("Gene" %in% colnames(counts_df)) {
        rownames(counts_df) <- counts_df$Gene
        counts_df$Gene <- NULL
      }
    }
    num_cols <- vapply(counts_df, is.numeric, logical(1))
    counts_df <- counts_df[, num_cols, drop = FALSE]
    as.matrix(counts_df)
  }

  get_layer_matrix <- function(layer_idx, genes) {
    counts <- NULL
    anno <- NULL
    if (!is.null(data_list)) {
      counts_name <- paste0("set", layer_idx, "_counts")
      anno_name <- paste0("set", layer_idx, "_anno")
      if (counts_name %in% names(data_list)) {
        counts <- data_list[[counts_name]]
      }
      if (anno_name %in% names(data_list)) {
        anno <- data_list[[anno_name]]
      }
    }

    if (!is.null(counts) && !is.null(anno)) {
      mat <- prep_counts_matrix(counts)
      if (is.null(mat)) return(NULL)

      if (is.null(rownames(anno)) || any(rownames(anno) == "")) {
        if ("Sample" %in% colnames(anno)) {
          rownames(anno) <- anno$Sample
        } else if ("ID" %in% colnames(anno)) {
          rownames(anno) <- anno$ID
        }
      }

      common_samples <- intersect(colnames(mat), rownames(anno))
      if (length(common_samples) == 0) {
        return(NULL)
      }
      mat <- mat[genes, common_samples, drop = FALSE]
      anno <- anno[common_samples, , drop = FALSE]

      cond_col <- condition_col
      if (is.null(cond_col)) {
        cond_col <- pick_condition_col(anno)
      }
      if (is.null(cond_col) || !cond_col %in% colnames(anno)) {
        stop("Could not determine condition column for layer ", layer_idx,
             ". Provide `condition_col`.")
      }
      groups <- as.character(anno[[cond_col]])
      if (!is.null(condition_order)) {
        groups <- factor(groups, levels = condition_order)
      }
      group_levels <- if (is.factor(groups)) levels(groups) else unique(groups)

      if (scale_rows && nrow(mat) > 1) {
        mat <- t(scale(t(mat)))
      }
      mat[is.na(mat)] <- 0

      agg <- sapply(group_levels, function(g) {
        rowMeans(mat[, groups == g, drop = FALSE], na.rm = TRUE)
      })
      if (is.null(dim(agg))) {
        agg <- matrix(agg, ncol = 1, dimnames = list(rownames(mat), group_levels))
      }
      colnames(agg) <- group_levels
      if (!is.null(condition_clean_regex)) {
        colnames(agg) <- gsub(condition_clean_regex, "", colnames(agg))
      }
      if (!is.null(condition_order)) {
        keep <- condition_order[condition_order %in% colnames(agg)]
        agg <- agg[, keep, drop = FALSE]
      }
      return(agg)
    }

    # Fallback to GFC matrix if counts/anno not available
    gfc <- layer_outputs[[paste0("set", layer_idx)]]$part2$GFC_all_genes
    rownames(gfc) <- gfc$Gene
    gfc$Gene <- NULL
    mat <- gfc[genes, , drop = FALSE]
    if (scale_rows && nrow(mat) > 1) {
      mat <- t(scale(t(mat)))
    }
    mat[is.na(mat)] <- 0
    mat
  }

  module_colors <- resolve_module_colors(modules)
  if (length(module_colors) == 0) {
    stop("No modules resolved. Check `modules` or `cluster_info`.")
  }

  if (!is.null(output_file)) {
    if (requireNamespace("Cairo", quietly = TRUE)) {
      Cairo::Cairo(file = output_file, width = width_in, height = height_in,
                   dpi = dpi, type = "pdf", units = "in")
    } else {
      pdf(file = output_file, width = width_in, height = height_in)
    }
  }

  for (mod_col in module_colors) {
    genes_str <- cluster_info$gene_n[cluster_info$color == mod_col][1]
    if (is.na(genes_str) || !nzchar(genes_str)) {
      next
    }
    genes <- unlist(strsplit(genes_str, ","))
    if (length(genes) == 0) next

    # Build per-layer matrices (rows = genes, cols = conditions)
    mats <- lapply(layer_indices, function(i) get_layer_matrix(i, genes))
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (length(mats) == 0) next

    # Harmonize condition columns across layers
    if (condition_mode == "intersect") {
      common_cols <- Reduce(intersect, lapply(mats, colnames))
      mats <- lapply(mats, function(m) m[, common_cols, drop = FALSE])
    } else {
      all_cols <- Reduce(union, lapply(mats, colnames))
      mats <- lapply(mats, function(m) {
        out <- matrix(NA_real_, nrow = nrow(m), ncol = length(all_cols),
                      dimnames = list(rownames(m), all_cols))
        out[, colnames(m)] <- m
        out
      })
    }

    all_genes <- genes
    # Optional ordering by mean absolute signal across all layers
    if (!is.null(max_genes) && length(all_genes) > max_genes && order_by == "mean_abs") {
      scores <- sapply(all_genes, function(g) {
        vals <- unlist(lapply(mats, function(m) if (g %in% rownames(m)) m[g, ] else rep(NA, ncol(m))))
        mean(abs(vals), na.rm = TRUE)
      })
      ord <- order(scores, decreasing = TRUE)
      all_genes <- all_genes[ord][seq_len(max_genes)]
    }

    mats <- lapply(seq_along(mats), function(idx) {
      m <- mats[[idx]]
      m <- m[intersect(rownames(m), all_genes), , drop = FALSE]
      suffix <- layers_names[layer_indices[idx]]
      rownames(m) <- paste0(rownames(m), "_", suffix)
      m
    })

    combined <- do.call(rbind, mats)
    combined[is.na(combined)] <- 0
    rownames_side <- "none"

    split_vec <- rep(layers_names[layer_indices], times = vapply(mats, nrow, integer(1)))

    if (is.null(zlim)) {
      if (scale_rows) {
        zlim <- c(-2, 0, 2)
      } else {
        max_abs <- max(abs(combined), na.rm = TRUE)
        if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
        zlim <- c(-max_abs, 0, max_abs)
      }
    }
    if (length(zlim) == 2) {
      zlim <- c(zlim[1], 0, zlim[2])
    }
    col_fun <- circlize::colorRamp2(zlim, c(low_col, mid_col, high_col))

    circlize::circos.clear()
    # Reduce gap if too many sectors
    if (nrow(combined) > 1) {
      max_gap <- 360 / nrow(combined) * 0.2
      if (gap_degree > max_gap) {
        gap_degree <- max_gap
      }
    }
    if (nrow(combined) > 200 && track_height < 0.6) {
      track_height <- 0.6
    }
    draw_heatmap <- function(th) {
      circlize::circos.par(start.degree = start_degree, gap.degree = gap_degree,
                           points.overflow.warning = FALSE, cell.padding = cell_padding,
                           track.margin = c(0, 0))
      circlize::circos.heatmap(combined, col = col_fun, na.col = na_col,
                               split = split_vec, track.height = th,
                               dend.side = "inside", rownames.side = rownames_side)
    }
    th <- track_height
    repeat {
      circlize::circos.clear()
      ok <- tryCatch({
        draw_heatmap(th)
        TRUE
      }, error = function(e) {
        if (grepl("not enough space for cells", e$message)) {
          FALSE
        } else {
          stop(e)
        }
      })
      if (ok) break
      th <- th * 0.7
      if (th < 0.05) {
        stop("Not enough space to draw circos heatmap. Reduce conditions or increase device size.")
      }
    }

    if (show_condition_labels) {
      if (is.null(condition_label_text)) {
        condition_label_text <- paste0(
          condition_label_prefix,
          paste(colnames(combined), collapse = ", ")
        )
      }
      if (is.null(condition_label_gp)) {
        condition_label_gp <- grid::gpar(cex = 0.7)
      }
      grid::pushViewport(grid::viewport())
      grid::grid.text(
        condition_label_text,
        x = if (is.null(legend_x)) unit(1, "npc") - unit(2, "mm") else legend_x,
        y = unit(0.95, "npc"),
        just = c("right", "top"),
        gp = condition_label_gp
      )
      grid::popViewport()
    }

    if (show_legend && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      if (is.null(legend_at)) {
        legend_at <- zlim
      }
      if (is.null(legend_labels)) {
        legend_labels <- legend_at
      }
      lgd <- ComplexHeatmap::Legend(
        col_fun = col_fun,
        title = legend_title,
        at = legend_at,
        labels = legend_labels
      )
      grid::pushViewport(grid::viewport())
      ComplexHeatmap::draw(
        lgd,
        x = if (is.null(legend_x)) unit(1, "npc") - unit(2, "mm") else legend_x,
        y = if (is.null(legend_y)) unit(0.5, "npc") else legend_y,
        just = legend_just
      )
      grid::popViewport()
    }

    if (show_gene_labels && nrow(combined) > 0) {
      circlize::circos.trackPlotRegion(
        track.index = 1,
        panel.fun = function(x, y) {
          sector.index <- circlize::get.cell.meta.data("sector.index")
          xcenter <- circlize::get.cell.meta.data("xcenter")
          ytop <- circlize::get.cell.meta.data("ylim")[2]
          circlize::circos.text(
            xcenter, ytop, sector.index,
            facing = "clockwise", niceFacing = TRUE, adj = c(0, 0),
            cex = label_cex
          )
        },
        bg.border = NA
      )
    }

    title(main = paste0(mod_col, " module"))
  }

  if (!is.null(output_file)) {
    dev.off()
  }
}
