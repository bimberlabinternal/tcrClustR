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

# internal palette generator shared by TCR network functions
.generate_tcrnetwork_palette <- function(n) {
  if (n == 0L) return(character(0L))
  if (n <= 8L) {
    RColorBrewer::brewer.pal(max(3L, n), "Set2")[seq_len(n)]
  } else if (n <= 12L) {
    RColorBrewer::brewer.pal(n, "Set3")
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(8L, "Set3"))(n)
  }
}

#' Build a TCR distance network graph
#'
#' @description
#' Constructs a \code{\link[tidygraph]{tbl_graph}} from a TCR pairwise distance matrix
#' stored in a tcrClustR Seurat object. Each node represents a unique clone; edges connect
#' pairs of clones whose pairwise distance is at or below \code{distanceThreshold}.
#' Clone metadata from the Seurat object is attached to each node, enabling downstream
#' visualisation with \code{\link{TCRDistanceNetwork}} or custom igraph analyses.
#'
#' @param seuratObj_TCR A Seurat object produced by \code{\link{CalculateTcrDistances}},
#'   containing distance assays in \code{@misc$TCR_Distances}.
#' @param chains Character. A single chain (\code{"TRA"}, \code{"TRB"}, \code{"TRG"},
#'   \code{"TRD"}) or a two-element vector (e.g. \code{c("TRA", "TRB")} for paired
#'   distances).
#' @param cdr3Only Logical. If \code{TRUE}, use CDR3-only distances; otherwise use
#'   full-length (CDR3 + V + J) distances. Default is \code{FALSE}.
#' @param distanceThreshold Numeric. Maximum pairwise distance for an edge to be drawn.
#'   Default is \code{50}.
#' @param edgeType Character. \code{"binary"} gives every retained edge weight 1.
#'   \code{"continuous"} adds \code{weight = 1/(distance+1)} and a \code{norm_weight}
#'   column normalised to the range 0--1. Default is \code{"binary"}.
#' @param verbose Logical. Print progress messages. Default is \code{FALSE}.
#'
#' @return A \code{\link[tidygraph]{tbl_graph}} (which is also a valid \code{igraph})
#'   with one node per clone. Node attributes include all metadata columns joinable via
#'   the chain's \code{CloneIdx} column. Edge attributes include \code{distance},
#'   \code{weight}, and (for \code{edgeType = "continuous"}) \code{norm_weight}.
#'
#' @seealso \code{\link{TCRDistanceNetwork}}, \code{\link{GetDistanceMatrix}}
#' @importFrom tidygraph tbl_graph
#' @importFrom rlang .data
#' @export
#' @examples
#' \dontrun{
#' seuratObj_TCR <- CalculateTcrDistances(seuratObj, chains = "TRB")
#' tcrGraph <- BuildTCRDistanceGraph(seuratObj_TCR, chains = "TRB", distanceThreshold = 50)
#' igraph::vcount(tcrGraph)
#' }
BuildTCRDistanceGraph <- function(
    seuratObj_TCR,
    chains,
    cdr3Only          = FALSE,
    distanceThreshold = 50,
    edgeType          = "binary",
    verbose           = FALSE
) {
  if (!edgeType %in% c("binary", "continuous")) {
    stop("edgeType must be 'binary' or 'continuous', got: ", edgeType)
  }
  if (!is.numeric(distanceThreshold) || length(distanceThreshold) != 1L || distanceThreshold <= 0) {
    stop("distanceThreshold must be a single positive number.")
  }

  chainPrefix <- if (length(chains) > 1L) .get_chain_field_prefix(chains) else chains
  if (verbose) message("Building TCR distance graph for chain prefix: ", chainPrefix)

  dist_m <- as.matrix(GetDistanceMatrix(seuratObj_TCR, chains, cdr3Only = cdr3Only))
  if (verbose) message("Distance matrix: ", nrow(dist_m), " x ", ncol(dist_m), " clones")

  # build upper-triangle edge list
  idx     <- which(upper.tri(dist_m), arr.ind = TRUE)
  edge_df <- data.frame(
    from     = rownames(dist_m)[idx[, 1L]],
    to       = colnames(dist_m)[idx[, 2L]],
    distance = dist_m[idx],
    stringsAsFactors = FALSE
  )
  edge_df <- edge_df[edge_df$distance <= distanceThreshold, , drop = FALSE]

  if (edgeType == "continuous") {
    edge_df$weight      <- 1 / (edge_df$distance + 1)
    edge_df$norm_weight <- if (nrow(edge_df) > 0L) {
      edge_df$weight / max(edge_df$weight)
    } else {
      numeric(0L)
    }
  } else {
    edge_df$weight <- 1L
  }

  if (verbose) message("Edges after thresholding: ", nrow(edge_df))

  # build node table, joining Seurat metadata (one row per clone)
  cloneIdxCol <- paste0(chainPrefix, "_CloneIdx")
  node_df     <- data.frame(name = rownames(dist_m), stringsAsFactors = FALSE)

  if (cloneIdxCol %in% colnames(seuratObj_TCR@meta.data)) {
    meta_src <- seuratObj_TCR@meta.data
    meta_src <- meta_src[!duplicated(meta_src[[cloneIdxCol]]), , drop = FALSE]
    node_df  <- dplyr::left_join(
      node_df,
      meta_src,
      by = stats::setNames(cloneIdxCol, "name")
    )
  } else if (verbose) {
    message("Note: ", cloneIdxCol, " not found in metadata — ",
            "nodes will carry no metadata annotations.")
  }

  tidygraph::tbl_graph(nodes = node_df, edges = edge_df, directed = FALSE)
}

#' Plot a TCR distance network graph
#'
#' @description
#' Visualises pairwise TCR distances as a network. Nodes represent unique clones; edges
#' connect clones whose pairwise distance is within \code{distanceThreshold}. Communities
#' can be detected via Louvain clustering (\code{communityMethod = "threshold"}) or read
#' from pre-computed DIANA results (\code{communityMethod = "DIANA"}, requires
#' \code{\link{RunTcrClustering}} to have been run first).
#'
#' @param seuratObj_TCR A Seurat object produced by \code{\link{CalculateTcrDistances}},
#'   containing distance assays in \code{@misc$TCR_Distances}.
#' @param chains Character. A single chain (\code{"TRA"}, \code{"TRB"}, \code{"TRG"},
#'   \code{"TRD"}) or a two-element vector (e.g. \code{c("TRA", "TRB")}).
#' @param cdr3Only Logical. If \code{TRUE}, use CDR3-only distances. Default \code{FALSE}.
#' @param communityMethod Character. How to assign community membership:
#'   \code{"threshold"} runs Louvain community detection on the thresholded graph;
#'   \code{"DIANA"} reads pre-computed DIANA cluster assignments from
#'   \code{RunTcrClustering}. Default \code{"threshold"}.
#' @param distanceThreshold Numeric. Maximum pairwise distance for drawing an edge.
#'   Default is \code{50}.
#' @param edgeType Character. \code{"binary"} draws uniform edges; \code{"continuous"}
#'   maps edge width and opacity to inverse distance. Default is \code{"binary"}.
#' @param colorBy Character or \code{NULL}. Metadata column name used for node fill
#'   colour (e.g. \code{"outcome"}, \code{"SubjectId"}). Ignored when
#'   \code{colorByCommunity = TRUE}. Default is \code{NULL}.
#' @param colorByCommunity Logical. If \code{TRUE}, node fill reflects community
#'   assignment rather than \code{colorBy}. Default is \code{FALSE}.
#' @param colorPalette Named character vector passed through as the colour scale. When
#'   \code{NULL} (default), an automatic RColorBrewer palette is generated.
#' @param labelBy Character or \code{NULL}. Metadata column name drawn as repelled text
#'   labels on each node. Default is \code{NULL} (no labels).
#' @param nodeSizeBy Character or \code{NULL}. Numeric metadata column used to scale node
#'   size (e.g. \code{"TRB_CloneSize"}). When \code{NULL}, all nodes use a fixed size.
#' @param layout Character. ggraph/igraph layout algorithm: \code{"fr"}
#'   (Fruchterman-Reingold, default), \code{"kk"} (Kamada-Kawai), \code{"circle"}, or
#'   \code{"nicely"}.
#' @param showIsolated Logical. When \code{FALSE} (default), clones with no edges within
#'   \code{distanceThreshold} are dropped. When \code{TRUE}, they are placed on a
#'   peripheral ring outside the connected component.
#' @param title Character or \code{NULL}. Plot title; auto-generated from assay name when
#'   \code{NULL}.
#' @param verbose Logical. Print progress messages. Default is \code{FALSE}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{\code{plot}}{A \code{ggplot}/\code{ggraph} object that can be further
#'       customised with standard ggplot2 layers.}
#'     \item{\code{layout}}{The \code{create_layout} data frame used to position nodes,
#'       containing x/y coordinates and all node attributes. Useful for reusing the same
#'       node positions across multiple plots or for direct inspection.}
#'   }
#'
#' @seealso \code{\link{BuildTCRDistanceGraph}}, \code{\link{RunTcrClustering}}
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_text
#'   create_layout theme_graph scale_edge_width_continuous scale_edge_alpha_continuous
#' @importFrom tidygraph activate as_tbl_graph
#' @importFrom igraph cluster_louvain membership degree E V delete_vertices
#'   ecount vcount induced_subgraph layout_with_fr
#' @importFrom rlang .data
#' @export
#' @examples
#' \dontrun{
#' seuratObj_TCR <- CalculateTcrDistances(seuratObj, chains = "TRB")
#' result <- TCRDistanceNetwork(seuratObj_TCR, chains = "TRB", distanceThreshold = 50,
#'                              colorBy = "outcome")
#' result$plot
#' head(result$layout[, c("x", "y", "name")])
#'
#' # Use DIANA communities (requires RunTcrClustering first)
#' seuratObj_TCR <- RunTcrClustering(seuratObj_TCR, chainsToCluster = "TRB")
#' result <- TCRDistanceNetwork(seuratObj_TCR, chains = "TRB", communityMethod = "DIANA",
#'                              colorByCommunity = TRUE)
#' result$plot
#' }
TCRDistanceNetwork <- function(
    seuratObj_TCR,
    chains,
    cdr3Only          = FALSE,
    communityMethod   = "DIANA",
    distanceThreshold = 50,
    edgeType          = "binary",
    colorBy           = NULL,
    colorByCommunity  = FALSE,
    colorPalette      = NULL,
    labelBy           = NULL,
    nodeSizeBy        = NULL,
    layout            = "fr",
    showIsolated      = FALSE,
    title             = NULL,
    verbose           = FALSE
) {
  # validate arguments
  if (!communityMethod %in% c("threshold", "DIANA")) {
    stop("communityMethod must be 'threshold' or 'DIANA', got: ", communityMethod)
  }
  if (!edgeType %in% c("binary", "continuous")) {
    stop("edgeType must be 'binary' or 'continuous', got: ", edgeType)
  }

  chainPrefix <- if (length(chains) > 1L) .get_chain_field_prefix(chains) else chains
  assayName   <- paste0(chainPrefix, "_", ifelse(cdr3Only, "cdr3", "fl"))

  if (communityMethod == "DIANA") {
    clusterIdxCol <- paste0(assayName, "_ClusterIdx")
    if (!clusterIdxCol %in% colnames(seuratObj_TCR@meta.data)) {
      stop(
        "communityMethod = 'DIANA' requires pre-computed cluster assignments. ",
        "Column '", clusterIdxCol, "' not found in metadata. ",
        "Run RunTcrClustering() first."
      )
    }
  }

  for (param_name in c("colorBy", "labelBy", "nodeSizeBy")) {
    param_val <- get(param_name)
    if (!is.null(param_val) && !param_val %in% colnames(seuratObj_TCR@meta.data)) {
      stop(param_name, " = '", param_val, "' is not a column in seuratObj_TCR@meta.data.")
    }
  }

  # build the network graph
  tcrGraph <- BuildTCRDistanceGraph(
    seuratObj_TCR     = seuratObj_TCR,
    chains            = chains,
    cdr3Only          = cdr3Only,
    distanceThreshold = distanceThreshold,
    edgeType          = edgeType,
    verbose           = verbose
  )

  # assign community membership
  g_full  <- tidygraph::as.igraph(tcrGraph)
  n_nodes <- igraph::vcount(g_full)

  if (communityMethod == "threshold") {
    if (igraph::ecount(g_full) == 0L) {
      if (verbose) message("No edges found within distanceThreshold; community set to NA for all nodes.")
      community_vec <- rep(NA_character_, n_nodes)
    } else {
      set.seed(42L)
      comms         <- igraph::cluster_louvain(g_full, weights = igraph::E(g_full)$weight)
      community_vec <- as.character(igraph::membership(comms))
    }
  } else {
    # DIANA: cluster column was already joined into node data by BuildTCRDistanceGraph
    clusterIdxCol <- paste0(assayName, "_ClusterIdx")
    node_data_df  <- as.data.frame(tcrGraph)
    community_vec <- as.character(node_data_df[[clusterIdxCol]])
    # noise clones (value 0) are treated as unassigned
    community_vec[!is.na(community_vec) & community_vec == "0"] <- NA_character_
  }

  tcrGraph <- tcrGraph %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(.community = community_vec)

  # handle isolated nodes (no edges within threshold)
  deg_full    <- igraph::degree(tidygraph::as.igraph(tcrGraph))
  is_isolated <- deg_full == 0L

  if (!showIsolated) {
    isolated_flags <- is_isolated
    tcrGraph <- tcrGraph %>%
      tidygraph::activate(nodes) %>%
      dplyr::mutate(.isolated_flag = isolated_flags) %>%
      dplyr::filter(!.data[[".isolated_flag"]]) %>%
      dplyr::select(-.data[[".isolated_flag"]])

    set.seed(42L)
    layout_obj          <- ggraph::create_layout(tcrGraph, layout = layout)
    layout_obj$isolated <- FALSE

  } else {
    conn_idx <- which(!is_isolated)
    iso_idx  <- which(is_isolated)
    n_iso    <- length(iso_idx)

    if (length(conn_idx) > 0L) {
      g_conn  <- igraph::induced_subgraph(tidygraph::as.igraph(tcrGraph), conn_idx)
      set.seed(42L)
      conn_xy <- igraph::layout_with_fr(
        g_conn,
        weights = if (edgeType == "continuous") igraph::E(g_conn)$weight else NULL
      )
      max_r <- max(sqrt(rowSums(conn_xy ^ 2))) * 1.5
    } else {
      conn_xy <- matrix(nrow = 0L, ncol = 2L)
      max_r   <- 5
    }

    layout_mat <- matrix(0, nrow = n_nodes, ncol = 2L)
    if (length(conn_idx) > 0L) layout_mat[conn_idx, ] <- conn_xy
    if (n_iso > 0L) {
      angles                <- seq(0, 2 * pi, length.out = n_iso + 1L)[seq_len(n_iso)]
      layout_mat[iso_idx, ] <- cbind(max_r * cos(angles), max_r * sin(angles))
    }

    layout_obj          <- ggraph::create_layout(
      tcrGraph, layout = "manual",
      x = layout_mat[, 1L], y = layout_mat[, 2L]
    )
    layout_obj$isolated <- is_isolated
  }

  # resolve fill variable and colour palette
  if (colorByCommunity) {
    fill_col   <- ".community"
    fill_label <- "Community"
  } else if (!is.null(colorBy)) {
    fill_col   <- colorBy
    fill_label <- colorBy
  } else {
    fill_col   <- NULL
    fill_label <- NULL
  }

  resolved_palette <- NULL
  if (!is.null(fill_col)) {
    fill_vals <- layout_obj[[fill_col]]
    if (!is.numeric(fill_vals)) {
      color_levels     <- unique(stats::na.omit(as.character(fill_vals)))
      n_colors         <- length(color_levels)
      resolved_palette <- if (is.null(colorPalette)) {
        stats::setNames(.generate_tcrnetwork_palette(n_colors), color_levels)
      } else {
        colorPalette
      }
    }
  }

  # assemble ggraph plot
  if (is.null(title)) {
    title <- paste0("TCR Distance Network: ", assayName)
  }
  subtitle <- paste0(
    "Threshold: ", distanceThreshold,
    "  |  Edge type: ", edgeType,
    "  |  Community: ", communityMethod
  )

  p <- ggraph::ggraph(layout_obj)

  # edges
  if (edgeType == "binary") {
    p <- p + ggraph::geom_edge_link(colour = "grey70", alpha = 0.5)
  } else {
    p <- p +
      ggraph::geom_edge_link(
        ggplot2::aes(width = norm_weight, alpha = norm_weight),
        colour  = "grey60",
        lineend = "round"
      ) +
      ggraph::scale_edge_width_continuous(
        range = c(0.3, 3),
        name  = "Closeness (1/dist)",
        guide = ggplot2::guide_legend(override.aes = list(alpha = 1))
      ) +
      ggraph::scale_edge_alpha_continuous(range = c(0.15, 0.9), guide = "none")
  }

  # nodes - four combinations of fill and size
  if (!is.null(fill_col) && !is.null(nodeSizeBy)) {
    p <- p + ggraph::geom_node_point(
      ggplot2::aes(fill = .data[[fill_col]], size = .data[[nodeSizeBy]]),
      shape = 21L, colour = "white", stroke = 0.7
    )
  } else if (!is.null(fill_col)) {
    p <- p + ggraph::geom_node_point(
      ggplot2::aes(fill = .data[[fill_col]]),
      shape = 21L, colour = "white", stroke = 0.7, size = 3
    )
  } else if (!is.null(nodeSizeBy)) {
    p <- p + ggraph::geom_node_point(
      ggplot2::aes(size = .data[[nodeSizeBy]]),
      shape = 21L, fill = "steelblue", colour = "white", stroke = 0.7
    )
  } else {
    p <- p + ggraph::geom_node_point(
      shape = 21L, fill = "steelblue", colour = "white", stroke = 0.7, size = 3
    )
  }

  # fill scale
  if (!is.null(fill_col)) {
    if (!is.null(resolved_palette)) {
      p <- p + ggplot2::scale_fill_manual(
        values   = resolved_palette,
        name     = fill_label,
        na.value = "grey70"
      )
    } else {
      # numeric fill - use a continuous colour scale
      p <- p + ggplot2::scale_fill_distiller(
        palette   = "Blues",
        direction = 1L,
        name      = fill_label,
        na.value  = "grey70"
      )
    }
  }

  # optional node size scale
  if (!is.null(nodeSizeBy)) {
    p <- p + ggplot2::scale_size_continuous(name = nodeSizeBy)
  }

  # optional text labels
  if (!is.null(labelBy)) {
    p <- p + ggraph::geom_node_text(
      ggplot2::aes(label = .data[[labelBy]]),
      repel    = TRUE,
      size     = 2.5,
      fontface = "bold"
    )
  }

  p <- p +
    ggplot2::labs(title = title, subtitle = subtitle) +
    ggraph::theme_graph(base_family = "sans", base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle   = ggplot2::element_text(size = 9, colour = "grey40"),
      legend.position = "right"
    )

  return(list(plot = p, layout = layout_obj))
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