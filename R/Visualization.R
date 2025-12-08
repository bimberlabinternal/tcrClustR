#' Create a ComplexHeatmap for a single TCR distance assay
#'
#' @description
#' Internal utility to generate a ComplexHeatmap for a single TCR distance assay.
#'
#' @param seuratObj_TCR A Seurat object containing TCR distance assays.
#' @param assay Character string specifying the assay name to plot.
#' @param cluster_info Factor vector of cluster assignments for each cell in the assay.
#' @param cluster_colors Named character vector of colors corresponding to cluster levels.
#' @param annotate_clusters Boolean specifying whether to display clustering information.
#' @param verbose Boolean controlling whether to display processing steps. Default is FALSE.
#'
#' @return A ComplexHeatmap object ready for drawing.
#' @importFrom ComplexHeatmap Heatmap draw HeatmapAnnotation rowAnnotation
#' @keywords internal
#' @examples
#' \dontrun{
#' .TCRDistanceHeatmap(
#'   seuratObj_TCR = seuratObj,
#'   assay = "TCRAssay",
#'   cluster_info = info_df,
#'   cluster_colors = color_list,
#'   annotate_clusters = TRUE,
#'   verbose = TRUE
#' )
#' }
.TCRDistanceHeatmap <- function(
  seuratObj_TCR = NULL,
  assay = NULL,
  cluster_info = NULL,
  cluster_colors = NULL,
  annotate_clusters = TRUE,
  verbose = FALSE
) {
  if (is.null(seuratObj_TCR)) {
    stop("seuratObj_TCR must not be NULL.")
  }
  if (is.null(assay)) {
    stop("assay must not be NULL.")
  }
  if (is.null(cluster_info)) {
    stop("cluster_info must not be NULL.")
  }
  if (is.null(cluster_colors)) {
    stop("cluster_colors must not be NULL.")
  }

  if (verbose) message("Creating heatmap for assay: ", assay)

  m <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))

  if (verbose) {
    message("Distance matrix dimensions: ", nrow(m), " x ", ncol(m))
    message("Number of cluster levels: ", length(levels(factor(cluster_info))))
    message("Cluster info length: ", length(cluster_info))
  }

  # Validate dimensions
  if (nrow(m) == 0 || ncol(m) == 0) {
    stop("Distance matrix has zero dimensions for assay ", assay, ": ", nrow(m), " x ", ncol(m))
  }

  if (length(cluster_info) != ncol(m)) {
    stop("Cluster info length (", length(cluster_info), ") does not match matrix columns (", ncol(m), ") for assay ", assay)
  }

  cluster_info <- factor(cluster_info)

  if (annotate_clusters) {
    col_annotation <- ComplexHeatmap::HeatmapAnnotation(
      cluster = cluster_info,
      col = list(cluster = cluster_colors),
      show_annotation_name = FALSE,
      annotation_name_side = "left",
      which = "column",
      show_legend = TRUE
    )

    row_annotation <- ComplexHeatmap::rowAnnotation(
      cluster = cluster_info,
      col     = list(cluster = cluster_colors),
      show_annotation_name = FALSE,
      show_legend          = FALSE
    )

    heatmap <- ComplexHeatmap::Heatmap(
      m,
      name               = assay,
      column_title       = assay,
      border_gp          = grid::gpar(col = "black", lty = 2),
      top_annotation     = col_annotation,
      left_annotation    = row_annotation,
      use_raster         = TRUE,
      cluster_columns    = TRUE,
      cluster_rows       = TRUE,
      row_dend_side      = "left",
      column_dend_side   = "top",
      column_split       = cluster_info,
      row_split          = cluster_info,
      show_heatmap_legend = TRUE,
      show_column_names    = FALSE,
      show_row_names       = FALSE
    )
  } else {
    heatmap <- ComplexHeatmap::Heatmap(
      m,
      name               = assay,
      column_title       = assay,
      border_gp          = grid::gpar(col = "black", lty = 2),
      use_raster         = TRUE,
      cluster_columns    = FALSE,
      cluster_rows       = FALSE,
      show_heatmap_legend = TRUE,
      show_column_names    = FALSE,
      show_row_names       = FALSE
    )
  }
  return(heatmap)
}

#' Generate and display heatmaps of TCR similarity for each assay
#'
#' @description
#' Generates and displays ComplexHeatmaps of TCR similarity for each assay in a Seurat object.
#'
#' @param seuratObj_TCR A Seurat object containing one or more TCR distance assays.
#' @param assayList Character vector of assay names to include. Default is NULL (all assays).
#' @param resolution Numeric clustering resolution parameter matching metadata column suffix.
#' @param annotate_clusters Boolean specifying whether to display clustering information.
#' @param verbose Boolean controlling whether to display processing steps. Default is FALSE.
#'
#' @return Returns the patchwork object containing ComplexHeatmaps of the TCR distance data.
#' @export
#' @examples
#' \dontrun{
#' TCRDistanceHeatmaps(
#'   seuratObj_TCR = seuratObj,
#'   assayList = NULL,
#'   resolution = 0.1,
#'   annotate_clusters = TRUE,
#'   verbose = TRUE
#' )
#' }
TCRDistanceHeatmaps <- function(
  seuratObj_TCR = NULL,
  assayList = NULL,
  resolution = 0.1,
  annotate_clusters = TRUE,
  verbose = FALSE
) {
  if (is.null(seuratObj_TCR)) {
    stop("Please provide a Seurat Object with TCR distance assays.")
  }

  if (verbose) {
    message("Starting TCRDistanceHeatmaps with resolution: ", resolution)
    message("Available assays in Seurat object: ", paste(SeuratObject::Assays(seuratObj_TCR), collapse = ", "))
    message("Available metadata columns: ", paste(colnames(seuratObj_TCR@meta.data), collapse = ", "))
  }

  assays_to_use <- if (is.null(assayList)) {
    SeuratObject::Assays(seuratObj_TCR)
  } else {
    assayList
  }

  if (verbose) message("Assays to process: ", paste(assays_to_use, collapse = ", "))

  heatmaps <- list()

  for (assay in assays_to_use) {
    if (verbose) message("Processing assay: ", assay)

    distance_matrix <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))

    if (verbose) message("Distance matrix dimensions for ", assay, ": ", nrow(distance_matrix), " x ", ncol(distance_matrix))

    # Find metadata column from clustering pipeline
    cluster_column <- paste0("TcrClustR_", assay, "_", resolution)
    if (verbose) message("Looking for cluster column: ", cluster_column)

    if (!(cluster_column %in% colnames(seuratObj_TCR@meta.data))) {
      if (verbose) message("Skipping ", assay, " - cluster column ", cluster_column, " not found")
      message(paste("Skipping", assay, "- no", cluster_column))
      next
    }

    # Cluster metadata is an array of assay_size x n_assay, so iterate and slice
    # the correct portion of that array
    full_cluster_info <- seuratObj_TCR@meta.data[[cluster_column]]
    n_cells_in_assay <- ncol(distance_matrix)

    if (verbose) {
      message("Full cluster info length: ", length(full_cluster_info))
      message("Cells in current assay (", assay, "): ", n_cells_in_assay)
    }

    assay_start_index <- 1
    for (a in SeuratObject::Assays(seuratObj_TCR)) {
      if (a == assay) break
      assay_cells <- ncol(Seurat::GetAssayData(seuratObj_TCR, assay = a, layer = "counts"))
      if (verbose) message("Assay ", a, " has ", assay_cells, " cells, start index now: ", assay_start_index + assay_cells)
      assay_start_index <- assay_start_index + assay_cells
    }
    assay_end_index <- assay_start_index + n_cells_in_assay - 1

    if (verbose) message("Slicing cluster info from index ", assay_start_index, " to ", assay_end_index)

    if (assay_end_index > length(full_cluster_info)) {
      if (verbose) message("ERROR: End index (", assay_end_index, ") exceeds cluster info length (", length(full_cluster_info), ")")
      message("Skipping ", assay, " - cluster indexing error")
      next
    }

    cluster_info <- full_cluster_info[assay_start_index:assay_end_index]
    cluster_info <- as.factor(cluster_info)
    cluster_levels <- levels(cluster_info)

    if (verbose) {
      message("Cluster levels for ", assay, ": ", paste(cluster_levels, collapse = ", "))
      message("Number of cluster levels: ", length(cluster_levels))
    }
    cluster_colors <- stats::setNames(
      if (length(cluster_levels) <= 8) {
        RColorBrewer::brewer.pal(length(cluster_levels), "Set2")
      } else if (length(cluster_levels) <= 12) {
        RColorBrewer::brewer.pal(length(cluster_levels), "Set3")
      } else {
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(cluster_levels))
      },
      cluster_levels
    )

    # Get a ComplexHeatmap
    if (verbose) message("Creating ComplexHeatmap for ", assay)

    heatmap_obj <- .TCRDistanceHeatmap(seuratObj_TCR,
                                       assay,
                                       cluster_info,
                                       cluster_colors,
                                       annotate_clusters,
                                       verbose = verbose)
    drawn_heatmap <- draw(heatmap_obj,
                          merge_legend = FALSE,
                          heatmap_legend_side = "right",
                          annotation_legend_side = "right",
                          newpage = FALSE)

    if (!is.null(drawn_heatmap)) {
      if (verbose) message("Successfully created heatmap for ", assay)
      heatmaps[[assay]] <- drawn_heatmap
    } else {
      if (verbose) message("Failed to create heatmap for ", assay)
    }
  }

  if (verbose) message("Created ", length(heatmaps), " heatmaps total")

  # Composite with patchwork
  combined_heatmaps <- patchwork::wrap_plots(lapply(heatmaps, function(hm) grid::grid.grabExpr(draw(hm))), ncol = 1)

  final_plot <- combined_heatmaps +
    patchwork::plot_annotation(title = "TCR Similarity",
                               theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)))

  #print(final_plot)

  return(final_plot)
}

#' Plot histograms of summed TCR distances by cluster for each assay
#'
#' @description
#' Plots histograms of summed TCR distances by cluster for each assay in a Seurat object.
#'
#' @param seuratObj_TCR A Seurat object with TCR distance assay data.
#' @param assayList Character vector of assays to plot. Default is NULL (all assays).
#' @param resolution Numeric clustering resolution matching metadata column suffix.
#'
#' @return Returns the patchwork object containing the TCR distance histograms.
#' @importFrom stats setNames
#' @export
#' @examples
#' \dontrun{
#' TCRDistanceHistograms(
#'   seuratObj_TCR = seuratObj,
#'   assayList = NULL,
#'   resolution = 0.1
#' )
#' }
TCRDistanceHistograms <- function(
  seuratObj_TCR = NULL,
  assayList  = NULL,
  resolution = 0.1
) {
  if (is.null(seuratObj_TCR)) {
    stop("Please provide a Seurat object with TCR distance assays.")
  }

  assays <- if (is.null(assayList)) {
    SeuratObject::Assays(seuratObj_TCR)
  } else {
    assayList
  }

  # precompute how many cells each assay has
  cell_counts <- stats::setNames(
    vapply(assays, function(a) {
      ncol(Seurat::GetAssayData(seuratObj_TCR, assay = a, layer = "counts"))
    }, integer(1)),
    assays
  )

  # compute start/end indices for slicing the metadata vector
  starts <- cumsum(c(1, head(cell_counts, -1)))
  ends   <- cumsum(cell_counts)

  plots <- list()
  for (i in seq_along(assays)) {
    assay <- assays[i]
    cluster_col <- paste0("TcrClustR_", assay, "_", resolution)

    if (!(cluster_col %in% colnames(seuratObj_TCR@meta.data))) {
      message("Skipping ", assay, ": no metadata column ", cluster_col)
      next
    }

    dist_mat  <- as.matrix(Seurat::GetAssayData(
      seuratObj_TCR, assay = assay, layer = "counts"
    ))

    # slice the full clustering vector down to this assay
    full_info    <- seuratObj_TCR@meta.data[[cluster_col]]
    cluster_info <- factor(full_info[starts[i]:ends[i]])

    # Build palette
    n_clust <- length(levels(cluster_info))
    pal     <- RColorBrewer::brewer.pal(min(n_clust, 8), "Set2")
    cl_cols <- stats::setNames(pal, levels(cluster_info))

    df <- data.frame(
      DistanceSum = rowSums(dist_mat),
      Cluster     = cluster_info
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = DistanceSum, fill = Cluster)) +
      ggplot2::geom_histogram(bins = 50, color = "black") +
      ggplot2::scale_fill_manual(values = cl_cols) +
      ggplot2::facet_wrap(~ Cluster, nrow = 1, scales = "free_y") +
      ggplot2::labs(
        title = assay,
        x     = "Summed Distances",
        y     = "# of cells"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(size = 14, face = "bold", hjust = 0),
        axis.title      = ggplot2::element_text(size = 12),
        axis.text       = ggplot2::element_text(size = 10),
        strip.text      = ggplot2::element_text(face = "bold"),
        legend.position = "none"
      )

    plots[[assay]] <- p
  }

  if (length(plots) == 0) {
    stop("No assays with the requested clustering metadata were found.")
  }

  combined <- patchwork::wrap_plots(plots, ncol = 1)
  print(combined)

  return(combined)
}

#' @title Get Example Markdown
#'
#' @description Save a template R markdown file, showing usage of this package
#' @param dest The destination filepath, where the file will be saved
#' @export
GetExampleMarkdown <- function(dest) {
  source <- system.file("rmd/tcrClustR.rmd", package = "tcrClustR")
  if (!file.exists(source)) {
    stop(paste0('Unable to find file: ', source))
  }

  dest <- normalizePath(dest, mustWork = F)
  success <- file.copy(source, dest, overwrite = TRUE)

  if (!success) {
    stop(paste0('Unable to copy file to: ', dest))
  }
}