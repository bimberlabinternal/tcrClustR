#' @include Utils.R

utils::globalVariables(
  names = c('group', 'distance', 'Cluster', 'Mu', 'Sigma',
            'MixingProportion', 'PointsPerCluster'),
  package = 'tcrClustR',
  add = TRUE
)

########################
### Helper Functions ###
########################

#' Derive chainId from an assay name
#' @description Strips the type suffix (fl/cdr3) from an assay name to get the chain prefix.
#' @param assayName Character string (e.g. "TRA_fl", "TRA_TRB_cdr3")
#' @return Character chain prefix (e.g. "TRA", "TRA_TRB")
#' @keywords internal
.assay_to_chain_id <- function(assayName) {
  assay_name_parts <- unlist(strsplit(assayName, split = '_'))
  paste0(head(assay_name_parts, -1), collapse = '_')
}


#' Extract Within-Group TCR Distance Vectors
#'
#' @description For each level of \code{splitField} in the Seurat metadata, maps
#' cells to their clone indices, subsets the clone-level distance matrix, and
#' extracts the upper triangle as a flat numeric vector.
#'
#' @param seuratObj A Seurat object produced by \code{CalculateTcrDistances}
#'   (must have \code{@misc$TCR_Distances}).
#' @param assayName Character. The distance assay to use (e.g., \code{"TRA_fl"},
#'   \code{"TRB_cdr3"}, \code{"TRA_TRB_fl"}).
#' @param splitField Character. A metadata column name whose levels define the
#'   groups.
#' @param minClonesPerGroup Integer. Minimum number of unique clones required in
#'   a group to produce a distance vector. Default 2.
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return A named list of numeric vectors, one per group level.
#' @export
ExtractGroupDistanceVectors <- function(seuratObj,
                                        assayName,
                                        splitField,
                                        minClonesPerGroup = 2,
                                        verbose = TRUE) {
  if (!'TCR_Distances' %in% names(seuratObj@misc)) {
    stop('seuratObj does not contain TCR_Distances in @misc. Run CalculateTcrDistances() first.')
  }

  if (!assayName %in% names(seuratObj@misc$TCR_Distances)) {
    stop(paste0(
      'Assay "', assayName, '" not found in @misc$TCR_Distances. Available: ',
      paste0(names(seuratObj@misc$TCR_Distances), collapse = ', ')
    ))
  }

  if (!splitField %in% names(seuratObj@meta.data)) {
    stop(paste0(
      'splitField "', splitField, '" not found in metadata. Available: ',
      paste0(sort(names(seuratObj@meta.data)), collapse = ', ')
    ))
  }

  clone_distance_matrix <- Seurat::GetAssayData(seuratObj@misc$TCR_Distances[[assayName]], layer = 'counts')
  clone_distance_matrix <- as.matrix(clone_distance_matrix)

  chain_id <- .assay_to_chain_id(assayName)
  clone_index_column <- paste0(chain_id, '_CloneIdx')

  if (!clone_index_column %in% names(seuratObj@meta.data)) {
    stop(paste0('Clone index column "', clone_index_column, '" not found in metadata.'))
  }

  group_levels <- unique(seuratObj@meta.data[[splitField]])
  group_levels <- group_levels[!is.na(group_levels)]

  distance_vectors_by_group <- list()
  for (current_group in group_levels) {
    cells_in_group_mask <- seuratObj@meta.data[[splitField]] == current_group &
      !is.na(seuratObj@meta.data[[splitField]])
    clones_for_group_cells <- seuratObj@meta.data[cells_in_group_mask, clone_index_column]
    unique_clones_in_group <- unique(stats::na.omit(clones_for_group_cells))
    unique_clones_in_group <- intersect(as.character(unique_clones_in_group), colnames(clone_distance_matrix))

    if (length(unique_clones_in_group) < minClonesPerGroup) {
      if (verbose) {
        message(paste0(
          'Group "', current_group, '": only ', length(unique_clones_in_group),
          ' clones (min ', minClonesPerGroup, '). Skipping.'
        ))
      }
      next
    }

    within_group_distance_matrix <- clone_distance_matrix[
      unique_clones_in_group, unique_clones_in_group, drop = FALSE
    ]
    pairwise_distance_vector <- within_group_distance_matrix[
      upper.tri(within_group_distance_matrix, diag = FALSE)
    ]

    if (verbose) {
      message(paste0(
        'Group "', current_group, '": ', length(unique_clones_in_group), ' clones, ',
        length(pairwise_distance_vector), ' pairwise distances.'
      ))
    }

    distance_vectors_by_group[[as.character(current_group)]] <- pairwise_distance_vector
  }

  attr(distance_vectors_by_group, 'assayName') <- assayName
  return(distance_vectors_by_group)
}


#' Downsample a distance vector with decile-stratified sampling
#'
#' @description Draws up to \code{maxSamples} values. Divides distances into
#' \code{nBins} equal-frequency quantile bins and draws up to
#' \code{samplesPerBin} from each, preserving the distributional shape.
#'
#' @param distances Numeric vector of pairwise distances.
#' @param maxSamples Integer. Hard cap on total returned samples.
#' @param nBins Integer. Number of equal-frequency bins. Default \code{10}.
#' @param samplesPerBin Integer. Samples per bin. Default \code{10}.
#' @param seed Integer. RNG seed. Default \code{42}.
#'
#' @return Numeric vector of length \eqn{\leq} \code{maxSamples}.
#' @keywords internal
.SampleDistanceVector <- function(distances,
                                   maxSamples = 100,
                                   nBins = 10,
                                   samplesPerBin = 10,
                                   seed = 42) {
  set.seed(seed)

  if (length(distances) <= maxSamples) {
    return(distances)
  }

  bin_assignments <- dplyr::ntile(distances, nBins)
  sampled_values <- c()
  for (bin_idx in seq_len(nBins)) {
    bin_values <- distances[bin_assignments == bin_idx]
    n_draw <- min(length(bin_values), samplesPerBin)
    sampled_values <- c(sampled_values, sample(bin_values, size = n_draw))
  }

  if (length(sampled_values) > maxSamples) {
    sampled_values <- sample(sampled_values, size = maxSamples)
  }

  return(sampled_values)
}


#' Extract cluster parameters from a fitted dirichletprocess model
#'
#' @description Pulls the DirichletProcess-discovered cluster means, standard deviations,
#' mixing proportions, and point counts from a fitted \code{dirichletprocess}
#' object.
#'
#' @param dpModel A fitted \code{dirichletprocess} object.
#' @param groupName Character label for this group (added as a column).
#'
#' @return A \code{data.frame} with columns:
#'   \code{Cluster}, \code{Mu}, \code{Sigma}, \code{MixingProportion},
#'   \code{PointsPerCluster}, \code{Group}.
#' @keywords internal
.ExtractClusterParams <- function(dpModel, groupName) {
  n_clusters <- dpModel$numberClusters
  mus    <- dpModel$clusterParameters[[1]][1, , ]
  sigmas <- dpModel$clusterParameters[[2]][1, , ]

  # When there is exactly one cluster, the indexing returns scalars
  if (n_clusters == 1) {
    mus    <- as.numeric(mus)
    sigmas <- as.numeric(sigmas)
  }

  data.frame(
    Cluster            = seq_len(n_clusters),
    Mu                 = as.numeric(mus),
    Sigma              = sqrt(as.numeric(sigmas)),
    MixingProportion   = as.numeric(dpModel$weights),
    PointsPerCluster   = as.integer(dpModel$pointsPerCluster),
    Group              = groupName,
    stringsAsFactors   = FALSE
  )
}


########################
### Public functions ###
########################

#' Dirichlet Process Cluster Analysis of TCR Distances
#'
#' @description High-level entry point for exploring the natural clustering
#' structure of pairwise TCR distances. For each level of \code{splitField},
#' the function extracts within-group distances from \code{assayName},
#' downsamples to \code{maxSamples} via decile-stratified sampling, fits a
#' Gaussian Dirichlet process mixture, and collects the per-cluster parameters
#' (mean distance, spread, mixing weight) into a tidy data frame.
#'
#' The resulting cluster means tell you where the natural modes of the distance
#' distribution sit, which directly informs the \code{dianaHeight} parameter in
#' \code{\link{RunTcrClustering}}.
#'
#' @param seuratObj A Seurat object produced by \code{CalculateTcrDistances}
#'   (must have \code{@misc$TCR_Distances}).
#' @param assayName Character. Distance assay to analyse (e.g., \code{"TRA_fl"},
#'   \code{"TRB_cdr3"}, \code{"TRA_TRB_fl"}).
#' @param splitField Character. Metadata column whose levels define groups
#'   (e.g., \code{"cDNA_ID"}, \code{"Tissue"}).
#' @param maxSamples Integer. Maximum pairwise distances to sample per group
#'   before fitting the DirichletProcess. Controls compute cost. Default \code{100}.
#' @param nIterations Integer. MCMC iterations for each DirichletProcess fit. Default
#'   \code{1000}.
#' @param minClonesPerGroup Integer. Groups with fewer clones are skipped.
#'   Default \code{2}.
#' @param nBins Integer. Number of equal-frequency quantile bins used during
#'   decile-stratified downsampling. Increase to preserve rare modes at the
#'   tails of the distance distribution. Default \code{10}.
#' @param samplesPerBin Integer. Samples drawn from each bin. Defaults to
#'   \code{max(1, floor(maxSamples / nBins))}, which evenly distributes the
#'   budget across bins. Override to draw more observations from each stratum.
#' @param seed Integer. RNG seed for reproducible downsampling. Default
#'   \code{42}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A list of class \code{"tcrDirichletResult"} containing:
#' \describe{
#'   \item{cluster_summary}{A tidy \code{data.frame} with one row per
#'     group-cluster combination. Columns: \code{Cluster}, \code{Mu},
#'     \code{Sigma}, \code{MixingProportion}, \code{PointsPerCluster},
#'     \code{Group}.}
#'   \item{models}{Named list of raw \code{dirichletprocess} model objects.}
#'   \item{assayName}{The assay used.}
#'   \item{splitField}{The metadata split field used.}
#' }
#'
#' @examples
#' \dontrun{
#' # after CalculateTcrDistances() + RunTcrClustering():
#' dp <- DirichletClusterAnalysis(seuratObj, "TRA_fl", "metadata_variable")
#'
#' # diagnostic plots
#' PlotClusterMeans(dp)
#' PlotMixingProportions(dp)
#'
#' # use cluster means to pick a dianaHeight:
#' dp$cluster_summary
#' 
#' # additionally, you can use quantile sampling to sample rare modes
#' dp <- DirichletClusterAnalysis(seuratObj, "TRA_fl", "metadata_variable", 
#'                                nBins = 10, samplesPerBin = 100)

#' 
#' }
#'
#' @importFrom dirichletprocess DirichletProcessGaussian Fit
#' @export
DirichletClusterAnalysis <- function(seuratObj,
                                      assayName,
                                      splitField,
                                      maxSamples = 100,
                                      nIterations = 1000,
                                      minClonesPerGroup = 2,
                                      nBins = 10,
                                      samplesPerBin = NULL,
                                      seed = 42,
                                      verbose = TRUE) {
  if (is.null(samplesPerBin)) {
    samplesPerBin <- max(1L, floor(maxSamples / nBins))
  }
  # extract per-group distance vectors
  group_distance_vectors <- ExtractGroupDistanceVectors(
    seuratObj         = seuratObj,
    assayName         = assayName,
    splitField        = splitField,
    minClonesPerGroup = minClonesPerGroup,
    verbose           = verbose
  )

  if (length(group_distance_vectors) == 0) {
    stop("No groups with sufficient clones found. Please check splitField and minClonesPerGroup.")
  }

  models_by_group <- list()
  param_frames    <- list()

  for (group_name in names(group_distance_vectors)) {
    if (verbose) message(paste0("Fitting DirichletProcess for group: ", group_name))

    raw_distances <- group_distance_vectors[[group_name]]

    # downsample
    fit_distances <- .SampleDistanceVector(
      distances     = raw_distances,
      maxSamples    = maxSamples,
      nBins         = nBins,
      samplesPerBin = samplesPerBin,
      seed          = seed
    )

    if (verbose && length(fit_distances) < length(raw_distances)) {
      message(paste0("  Downsampled ", length(raw_distances), " -> ",
                     length(fit_distances), " distances."))
    }

    if (length(unique(fit_distances)) < 2) {
      warning(paste0("Group '", group_name, "': constant distances; skipping."))
      next
    }

    # fit DirichletProcess Gaussian mixture
    dp_model <- tryCatch({
      m <- dirichletprocess::DirichletProcessGaussian(fit_distances)
      dirichletprocess::Fit(m, nIterations)
    }, error = function(e) {
      warning(paste0("DirichletProcess fitting failed for group '", group_name, "': ", e$message))
      NULL
    })

    if (is.null(dp_model)) next

    models_by_group[[group_name]] <- dp_model

    # finally, extract cluster params
    param_frames[[group_name]] <- .ExtractClusterParams(dp_model, group_name)
  }

  if (length(models_by_group) == 0) {
    stop("No groups produced a valid Dirichlet process model.")
  }

  cluster_summary <- do.call(rbind, param_frames)
  rownames(cluster_summary) <- NULL

  result <- list(
    cluster_summary = cluster_summary,
    models          = models_by_group,
    assayName       = assayName,
    splitField      = splitField
  )
  class(result) <- c("tcrDirichletResult", "list")

  return(result)
}


#' Plot Cluster Mean Distances with Error Bars
#'
#' @description Produces the primary diagnostic plot: for each group, the
#' DirichletProcess-discovered cluster means (mu) are shown as points with +/- 1 sigma
#' error bars. The x-axis is cluster index, and groups are jitter-dodged so
#' multiple groups can be compared. The y-axis (\code{Mu}) directly corresponds
#' to average TCR distances and can be used to select a \code{dianaHeight}
#' cutoff.
#'
#' @param dpResult A \code{tcrDirichletResult} object from
#'   \code{\link{DirichletClusterAnalysis}}.
#' @param verbose Logical. Default \code{FALSE}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' dp <- DirichletClusterAnalysis(seuratObj, "TRA_fl", "cDNA_ID")
#' PlotClusterMeans(dp)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_errorbar
#'   position_jitterdodge scale_x_continuous scale_y_continuous
#'   labs theme_bw
#' @export
PlotClusterMeans <- function(dpResult, verbose = FALSE) {
  if (!inherits(dpResult, "tcrDirichletResult")) {
    stop("dpResult must be a tcrDirichletResult from DirichletClusterAnalysis().")
  }

  df <- dpResult$cluster_summary

  #natural-sort the group levels
  group_levels <- unique(df$Group)
  if (requireNamespace("naturalsort", quietly = TRUE)) {
    group_levels <- naturalsort::naturalsort(group_levels)
  } else {
    group_levels <- sort(group_levels)
  }
  df$Group <- factor(df$Group, levels = group_levels)

  jd_pos <- ggplot2::position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75)

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x     = Cluster,
    y     = Mu,
    color = Group
  )) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Mu - Sigma, ymax = Mu + Sigma),
      width    = 0.2,
      position = jd_pos,
      lineend  = "square"
    ) +
    ggplot2::geom_point(position = jd_pos) +
    ggplot2::scale_x_continuous(breaks = seq_len(max(df$Cluster))) +
    ggplot2::scale_y_continuous(
      breaks = seq(0, max(df$Mu + df$Sigma, na.rm = TRUE), by = 10)
    ) +
    ggplot2::labs(
      title = paste0("DirichletProcess Cluster Means: ", dpResult$assayName),
      subtitle = paste0("Split by: ", dpResult$splitField),
      x = "Cluster Index",
      y = "Average TCR Distance"
    ) +
    egg::theme_article()

  return(p)
}


#' Plot Mixing Proportions per Cluster
#'
#' @description Companion to \code{\link{PlotClusterMeans}}. Shows the
#' proportion of data points assigned to each DirichletProcess cluster as a dodged bar
#' chart, one bar per group. Clusters with high mixing weight and low mean
#' distance identify well-supported clonotype families.
#'
#' @param dpResult A \code{tcrDirichletResult} object from
#'   \code{\link{DirichletClusterAnalysis}}.
#' @param verbose Logical. Default \code{FALSE}.
#'
#' @return A \pkg{ggplot2} object.
#'
#' @examples
#' \dontrun{
#' dp <- DirichletClusterAnalysis(seuratObj, "TRA_fl", "cDNA_ID")
#' PlotMixingProportions(dp)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_col position_dodge2
#'   scale_x_continuous scale_y_continuous labs theme_bw
#' @export
PlotMixingProportions <- function(dpResult, verbose = FALSE) {
  if (!inherits(dpResult, "tcrDirichletResult")) {
    stop("dpResult must be a tcrDirichletResult from DirichletClusterAnalysis().")
  }

  df <- dpResult$cluster_summary

  group_levels <- unique(df$Group)
  if (requireNamespace("naturalsort", quietly = TRUE)) {
    group_levels <- naturalsort::naturalsort(group_levels)
  } else {
    group_levels <- sort(group_levels)
  }
  df$Group <- factor(df$Group, levels = group_levels)

  p <- ggplot2::ggplot(df, ggplot2::aes(
    x    = Cluster,
    y    = MixingProportion,
    fill = Group
  )) +
    ggplot2::geom_col(position = ggplot2::position_dodge2(preserve = "single")) +
    ggplot2::scale_x_continuous(breaks = seq_len(max(df$Cluster))) +
    ggplot2::scale_y_continuous(
      labels = scales::percent_format(),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      title = paste0(" Mixing Proportions: ", dpResult$assayName),
      subtitle = paste0("Split by: ", dpResult$splitField),
      x = "Cluster Index",
      y = "Proportional Cluster Assignment"
    ) +
    egg::theme_article()

  return(p)
}
