#' @include Visualization.R
NULL

# ---------------------------------------------------------------------------
# internal: fit a GLM or GLMM on dyad-level data
# ---------------------------------------------------------------------------
.FitDyadModel <- function(data, response, family, predictor_term,
                          approach, from_subj, to_subj, verbose = FALSE) {
  if (approach == "clustered") {
    data$Dyad_ID <- paste0(
      pmin(data[[from_subj]], data[[to_subj]]), "__",
      pmax(data[[from_subj]], data[[to_subj]])
    )
    fml_str <- paste0(response, " ~ ", predictor_term)
    if (verbose) message("Fitting GLM (clustered SEs): ", fml_str)
    fit <- stats::glm(stats::as.formula(fml_str), data = data, family = family)
  } else {
    rand_term <- paste0("(1|", from_subj, ") + (1|", to_subj, ")")
    fml_str   <- paste0(response, " ~ ", predictor_term, " + ", rand_term)
    if (verbose) message("Fitting GLMM: ", fml_str)
    fit <- lme4::glmer(stats::as.formula(fml_str), data = data, family = family)
  }
  list(model = fit, data = data)
}

# ---------------------------------------------------------------------------
# internal: marginal estimates via marginaleffects
# ---------------------------------------------------------------------------
.EstimateMarginalEffects <- function(fit_result, predictor_term,
                                     approach, estimand = "effects",
                                     from_var = NULL, to_var = NULL) {
  vcov_arg <- if (approach == "clustered") ~Dyad_ID else NULL

  if (predictor_term == "pair_type") {
    if (estimand == "effects") {
      me <- marginaleffects::avg_comparisons(
        fit_result$model,
        variables = "pair_type",
        vcov      = vcov_arg,
        newdata   = fit_result$data
      )
    } else {
      me <- marginaleffects::avg_predictions(
        fit_result$model,
        by       = "pair_type",
        vcov     = vcov_arg,
        newdata  = fit_result$data
      )
    }
  } else {
    me <- marginaleffects::avg_slopes(
      fit_result$model,
      variables = c(from_var, to_var),
      vcov      = vcov_arg,
      newdata   = fit_result$data
    )
  }
  as.data.frame(me)
}

# ---------------------------------------------------------------------------
# internal: build the unpaired per-clone connection summary
# ---------------------------------------------------------------------------
.BuildUnpairedSummary <- function(dyad_df, from_var, to_var, groupColumn) {
  clone_groups <- unique(dyad_df[, c("from_clone", from_var), drop = FALSE])
  names(clone_groups) <- c("clone", groupColumn)
  connected_dyads <- dyad_df[dyad_df$connected == 1L, , drop = FALSE]
  if (nrow(connected_dyads) > 0L) {
    from_partners <- connected_dyads[, c("from_clone", to_var), drop = FALSE]
    to_partners   <- connected_dyads[, c("to_clone",   from_var), drop = FALSE]
    names(from_partners) <- c("clone", "partner_group")
    names(to_partners)   <- c("clone", "partner_group")
    all_partners <- rbind(from_partners, to_partners)
    all_partners <- all_partners[!is.na(all_partners$partner_group), , drop = FALSE]
    partner_summary <- do.call(rbind, lapply(
      split(all_partners, all_partners$clone),
      function(x) data.frame(
        clone          = x$clone[[1L]],
        partner_groups = paste(sort(unique(x$partner_group)), collapse = ","),
        stringsAsFactors = FALSE
      )
    ))
    out <- merge(clone_groups, partner_summary, by = "clone", all.x = TRUE)
  } else {
    out <- clone_groups
    out$partner_groups <- NA_character_
  }
  out$is_connected <- !is.na(out$partner_groups)
  out <- out[, c("clone", groupColumn, "is_connected", "partner_groups"), drop = FALSE]
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# internal: compute per-group unpaired clone rates with binomial CIs
# ---------------------------------------------------------------------------
.ComputeUnpairedRates <- function(unpaired_summary, groupColumn) {
  groups <- split(unpaired_summary, unpaired_summary[[groupColumn]])
  do.call(rbind, lapply(names(groups), function(g) {
    df <- groups[[g]]
    n <- nrow(df)
    n_unpaired <- sum(!df$is_connected)
    bt <- stats::binom.test(n_unpaired, n)
    data.frame(
      group      = g,
      estimate   = as.numeric(bt$estimate),
      conf.low   = bt$conf.int[1L],
      conf.high  = bt$conf.int[2L],
      n_clones   = n,
      n_unpaired = n_unpaired,
      stringsAsFactors = FALSE,
      row.names  = NULL
    )
  }))
}

#' @title Model Pairwise TCR Dyad Relationships
#'
#' @description Fits two models over all unique clone pairs derived from a TCR
#'   distance matrix and returns marginal estimates for the pairwise
#'   relationship between levels of \code{groupColumn}.
#'
#'   **Connectivity model** (binomial): the response is whether a pair is
#'   connected. When \code{communityMethod = "DIANA"}, DIANA clustering is
#'   run internally on the distance matrix at \code{distanceThreshold}
#'   (used as the dendrogram cut height), and two clones are connected if
#'   they belong to the same cluster. When
#'   \code{communityMethod = "threshold"}, a clone pair is connected if its
#'   tcrdist distance is \eqn{\leq} \code{distanceThreshold}.
#'
#'   **Distance model** (Poisson): the response is the raw tcrdist distance
#'   among connected pairs.
#'
#'   When \code{estimand = "effects"} (default), marginal effects are
#'   estimated via \code{\link[marginaleffects]{avg_comparisons}}, giving
#'   pairwise contrasts between \code{pair_type} levels. When
#'   \code{estimand = "means"}, marginal means are estimated via
#'   \code{\link[marginaleffects]{avg_predictions}}, giving predicted values
#'   for each \code{pair_type} level.
#'
#'   For categorical \code{groupColumn}, an unordered dyad-type factor
#'   \code{pair_type} is constructed from the lexicographically sorted
#'   combination of the two endpoint group labels, e.g. \code{"G1:G1"},
#'   \code{"G1:G2"}. Because dyads are symmetric (\code{G1:G2 == G2:G1}),
#'   each unique pairing receives one level. For numeric \code{groupColumn},
#'   \code{from_<var>} and \code{to_<var>} enter as separate additive fixed
#'   effects.
#'
#'   When \code{approach = "mixed"} (default), crossed random intercepts
#'   for the subject owning each endpoint are included via
#'   \code{\link[lme4]{glmer}}. When \code{approach = "clustered"}, a plain
#'   GLM is fitted and cluster-robust standard errors (clustered on a sorted
#'   subject-pair identifier \code{Dyad_ID}) are computed via
#'   \code{marginaleffects}.
#'
#' @param seuratObj_TCR A Seurat object with TCR distance matrices stored in
#'   \code{@misc$TCR_Distances} (output of \code{CalculateTcrDistances}).
#' @param chains Character vector of chain(s), e.g. \code{"TRB"} or
#'   \code{c("TRA","TRB")}.
#' @param groupColumn Character. Seurat metadata column defining the groups
#'   to compare, e.g. a treatment arm (\code{"Treatment"}) or phenotype
#'   (\code{"Phenotype"}).
#' @param distanceThreshold Numeric. Interpretation depends on
#'   \code{communityMethod}. For \code{"DIANA"}: the dendrogram cut height
#'   passed to \code{.DianaClustering} (analogous to \code{dianaHeight}
#'   in \code{RunTcrClustering}). For \code{"threshold"}: the maximum
#'   pairwise distance at which two clones are considered connected.
#' @param communityDistanceThreshold Numeric or \code{NULL}. When
#'   \code{communityMethod = "threshold"}, the maximum pairwise distance used
#'   to build the Louvain community graph. When \code{NULL} (default),
#'   \code{distanceThreshold} is reused. Ignored when
#'   \code{communityMethod = "DIANA"}.
#' @param cdr3Only Logical. Use CDR3-only distances. Default \code{FALSE}.
#' @param estimand Character. \code{"effects"} (default) estimates marginal
#'   effects via \code{\link[marginaleffects]{avg_comparisons}}, returning
#'   pairwise contrasts between \code{pair_type} levels. \code{"means"}
#'   estimates marginal means via
#'   \code{\link[marginaleffects]{avg_predictions}}, returning predicted
#'   values for each level.
#' @param communityMethod Character. \code{"DIANA"} (default) runs DIANA
#'   hierarchical clustering internally at \code{distanceThreshold} and
#'   defines connectivity by cluster co-membership. \code{"threshold"}
#'   defines connectivity via a raw distance cutoff and derives Louvain
#'   communities for annotation.
#' @param clusterSizeThreshold Integer. Minimum number of clones in a DIANA
#'   cluster for it to be retained; smaller clusters are treated as
#'   unclustered. Only used when \code{communityMethod = "DIANA"}.
#'   Default \code{2}.
#' @param referenceLevel Character or \code{NULL}. For categorical
#'   \code{groupColumn}, the \code{pair_type} level to use as the baseline.
#'   When \code{NULL}, R's default factor ordering is used.
#' @param approach Character. \code{"mixed"} (default) fits GLMMs with
#'   crossed random intercepts for subject at each endpoint.
#'   \code{"clustered"} fits plain GLMs with cluster-robust standard errors
#'   via \code{marginaleffects}.
#' @param subjectIdCol Character. Metadata column containing subject IDs.
#'   Default \code{"SubjectId"}.
#' @param verbose Logical. Emit progress messages. Default \code{FALSE}.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{connectivity}}{A data frame of marginal estimates for the
#'       connectivity model. When \code{estimand = "effects"}, contains
#'       pairwise contrasts between \code{pair_type} levels. When
#'       \code{estimand = "means"}, contains predicted probabilities per
#'       level.}
#'     \item{\code{distance}}{A data frame of marginal estimates for the
#'       distance model.}
#'     \item{\code{models}}{A list with elements \code{connectivity} and
#'       \code{distance}, each containing \code{model} (the fitted model
#'       object) and \code{data} (the data frame used for fitting).}
#'     \item{\code{unpaired_summary}}{For categorical \code{groupColumn}, a
#'       data frame with one row per clone giving connectivity partner
#'       information. \code{NULL} for numeric predictors.}
#'     \item{\code{unpaired_rates}}{For categorical \code{groupColumn}, a
#'       data frame with one row per group level giving the proportion of
#'       clones that have no connections, with exact binomial 95\% confidence
#'       intervals. Columns: \code{group}, \code{estimate}, \code{conf.low},
#'       \code{conf.high}, \code{n_clones}, \code{n_unpaired}.
#'       \code{NULL} for numeric predictors.}
#'   }
#'
#' @seealso \code{\link{CalculateTcrDistances}}, \code{\link{TCRDistanceNetwork}}
#' @importFrom lme4 glmer
#' @importFrom stats as.formula setNames glm relevel binom.test
#' @importFrom igraph cluster_louvain membership V ecount vcount components
#' @importFrom tidygraph activate
#' @importFrom marginaleffects avg_predictions avg_comparisons avg_slopes
#' @export
ModelTCRDyads <- function(
    seuratObj_TCR,
    chains,
    groupColumn,
    distanceThreshold,
    communityDistanceThreshold  = NULL,
    cdr3Only                    = FALSE,
    estimand                    = "effects",
    approach                    = "mixed",
    communityMethod             = "DIANA",
    clusterSizeThreshold        = 2,
    referenceLevel              = NULL,
    subjectIdCol                = "SubjectId",
    verbose                     = FALSE
) {
  if (!estimand %in% c("effects", "means")) {
    stop("estimand must be 'effects' or 'means', got: ", estimand)
  }
  if (!approach %in% c("mixed", "clustered")) {
    stop("approach must be 'mixed' or 'clustered', got: ", approach)
  }
  if (!communityMethod %in% c("DIANA", "threshold")) {
    stop("communityMethod must be 'DIANA' or 'threshold', got: ", communityMethod)
  }

  # resolve column names
  chainPrefix <- if (length(chains) > 1L) .get_chain_field_prefix(chains) else chains
  cloneIdxCol <- paste0(chainPrefix, "_CloneIdx")
  dist_m      <- as.matrix(GetDistanceMatrix(seuratObj_TCR, chains, cdr3Only = cdr3Only))

  # build full dyad table (upper triangle)
  idx     <- which(upper.tri(dist_m), arr.ind = TRUE)
  dyad_df <- data.frame(
    from_clone = rownames(dist_m)[idx[, 1L]],
    to_clone   = colnames(dist_m)[idx[, 2L]],
    distance   = dist_m[idx],
    stringsAsFactors = FALSE
  )

  # join Seurat metadata for both endpoints
  if (cloneIdxCol %in% colnames(seuratObj_TCR@meta.data)) {
    meta_src         <- seuratObj_TCR@meta.data
    meta_src         <- meta_src[!duplicated(meta_src[[cloneIdxCol]]), , drop = FALSE]
    from_meta        <- meta_src
    to_meta          <- meta_src
    names(from_meta) <- paste0("from_", names(from_meta))
    names(to_meta)   <- paste0("to_",   names(to_meta))
    from_key         <- paste0("from_", cloneIdxCol)
    to_key           <- paste0("to_",   cloneIdxCol)
    dyad_df <- dplyr::left_join(
      dyad_df, from_meta,
      by = stats::setNames(from_key, "from_clone")
    )
    dyad_df <- dplyr::left_join(
      dyad_df, to_meta,
      by = stats::setNames(to_key, "to_clone")
    )
  }

  # drop dyads where either endpoint lacks metadata for model columns
  from_var  <- paste0("from_", groupColumn)
  to_var    <- paste0("to_",   groupColumn)
  from_subj <- paste0("from_", subjectIdCol)
  to_subj   <- paste0("to_",   subjectIdCol)
  required_cols <- c(from_var, to_var, from_subj, to_subj)
  present_cols  <- intersect(required_cols, names(dyad_df))
  if (length(present_cols) > 0L) {
    complete <- stats::complete.cases(dyad_df[, present_cols, drop = FALSE])
    n_drop   <- sum(!complete)
    if (n_drop > 0L) {
      if (verbose) {
        message("Dropping ", n_drop, " dyads (out of ", nrow(dyad_df),
                ") where endpoint metadata is NA for ", groupColumn,
                " or ", subjectIdCol)
      }
      dyad_df <- dyad_df[complete, , drop = FALSE]
    }
  }

  # define connectivity and annotate community membership
  if (communityMethod == "DIANA") {
    # run DIANA clustering internally using distanceThreshold as cut height
    diana_result <- .DianaClustering(dist_m, cutHeight = distanceThreshold, verbose = verbose)
    cluster_vec  <- .ThresholdClustersBySize(diana_result$clustering,
                                             minClusterSize = clusterSizeThreshold,
                                             verbose = verbose)
    if (all(cluster_vec == 0L)) {
      stop("No DIANA clusters formed at distanceThreshold = ", distanceThreshold,
           " with clusterSizeThreshold = ", clusterSizeThreshold,
           ". Try increasing distanceThreshold or lowering clusterSizeThreshold.")
    }
    # map cluster assignments onto dyad endpoints
    dyad_df$from_diana_community <- as.character(cluster_vec[dyad_df$from_clone])
    dyad_df$to_diana_community   <- as.character(cluster_vec[dyad_df$to_clone])
    # treat cluster 0 (unclustered) as NA
    dyad_df$from_diana_community[dyad_df$from_diana_community == "0"] <- NA_character_
    dyad_df$to_diana_community[dyad_df$to_diana_community == "0"]     <- NA_character_
    dyad_df$same_diana_community <- as.integer(
      !is.na(dyad_df$from_diana_community) & !is.na(dyad_df$to_diana_community) &
        dyad_df$from_diana_community == dyad_df$to_diana_community
    )
    # connectivity = same DIANA cluster
    dyad_df$connected <- dyad_df$same_diana_community
    dyad_df$connected[is.na(dyad_df$connected)] <- 0L
  } else {
    dyad_df$connected <- as.integer(dyad_df$distance <= distanceThreshold)
    community_dist <- if (is.null(communityDistanceThreshold)) distanceThreshold else communityDistanceThreshold
    tg <- BuildTCRDistanceGraph(
      seuratObj_TCR, chains = chains, cdr3Only = cdr3Only,
      distanceThreshold = community_dist, edgeType = "binary", verbose = FALSE
    )
    g_comm <- tidygraph::as.igraph(tg)
    if (igraph::ecount(g_comm) > 0L) {
      set.seed(42L)
      comms    <- igraph::cluster_louvain(g_comm)
      comm_vec <- stats::setNames(as.character(igraph::membership(comms)),
                                  igraph::V(g_comm)$name)
    } else {
      comm_vec <- stats::setNames(rep(NA_character_, igraph::vcount(g_comm)),
                                  igraph::V(g_comm)$name)
    }
    dyad_df$from_community <- comm_vec[dyad_df$from_clone]
    dyad_df$to_community   <- comm_vec[dyad_df$to_clone]
    dyad_df$same_community <- as.integer(
      !is.na(dyad_df$from_community) & !is.na(dyad_df$to_community) &
        dyad_df$from_community == dyad_df$to_community
    )
    if (sum(dyad_df$connected) == 0L) {
      stop("No connected pairs at distanceThreshold = ", distanceThreshold,
           ". No clones are within this distance. Try increasing distanceThreshold.")
    }
  }

  if (verbose) {
    message("Total dyads: ", nrow(dyad_df),
            " | connected: ", sum(dyad_df$connected))
  }

  # construct predictor term
  if (!from_var %in% names(dyad_df)) {
    stop("groupColumn '", groupColumn, "' not found in Seurat metadata")
  }

  if (is.numeric(dyad_df[[from_var]])) {
    predictor_term   <- paste(from_var, to_var, sep = " + ")
    unpaired_summary <- NULL
  } else {
    g_from <- as.character(dyad_df[[from_var]])
    g_to   <- as.character(dyad_df[[to_var]])
    pair_labels <- ifelse(
      is.na(g_from) | is.na(g_to),
      NA_character_,
      paste(pmin(g_from, g_to), pmax(g_from, g_to), sep = ":")
    )
    all_levels <- sort(unique(pair_labels[!is.na(pair_labels)]))
    dyad_df$pair_type <- factor(pair_labels, levels = all_levels)
    if (!is.null(referenceLevel)) {
      if (!referenceLevel %in% all_levels) {
        stop("referenceLevel '", referenceLevel, "' is not a valid pair_type level. ",
             "Available levels: ", paste(all_levels, collapse = ", "))
      }
      dyad_df$pair_type <- stats::relevel(dyad_df$pair_type, ref = referenceLevel)
    }
    predictor_term   <- "pair_type"
    unpaired_summary <- .BuildUnpairedSummary(dyad_df, from_var, to_var, groupColumn)
  }

  # compute per-group unpaired rates
  unpaired_rates <- if (!is.null(unpaired_summary)) {
    .ComputeUnpairedRates(unpaired_summary, groupColumn)
  } else {
    NULL
  }

  # validate subject column
  if (!from_subj %in% names(dyad_df)) {
    stop("subjectIdCol '", subjectIdCol, "' not found in Seurat metadata; ",
         "ensure the column is present before fitting models.")
  }

  # --- connectivity model (binomial, all dyad pairs) ---
  if (verbose) message("Fitting connectivity model...")
  conn_fit <- .FitDyadModel(
    data           = dyad_df,
    response       = "connected",
    family         = "binomial",
    predictor_term = predictor_term,
    approach       = approach,
    from_subj      = from_subj,
    to_subj        = to_subj,
    verbose        = verbose
  )

  # --- build edge subset for distance model ---
  # for DIANA, connected already means intra-cluster so no extra filtering needed
  edges <- dyad_df[dyad_df$connected == 1L, , drop = FALSE]

  if (nrow(edges) == 0L) {
    stop("No connected edges remain for the distance model. ",
         "Try increasing distanceThreshold, lowering clusterSizeThreshold, ",
         "or changing communityMethod.")
  }

  # propagate pair_type levels (including referenceLevel) to the edge subset
  if ("pair_type" %in% names(edges)) {
    edges$pair_type <- factor(edges$pair_type, levels = levels(dyad_df$pair_type))
  }

  # propagate Dyad_ID if already computed by the connectivity model
  if (approach == "clustered" && !"Dyad_ID" %in% names(edges)) {
    edges$Dyad_ID <- paste0(
      pmin(edges[[from_subj]], edges[[to_subj]]), "__",
      pmax(edges[[from_subj]], edges[[to_subj]])
    )
  }

  if (verbose) message("Fitting distance model on ", nrow(edges), " edges...")
  dist_fit <- .FitDyadModel(
    data           = edges,
    response       = "distance",
    family         = "poisson",
    predictor_term = predictor_term,
    approach       = approach,
    from_subj      = from_subj,
    to_subj        = to_subj,
    verbose        = verbose
  )

  # --- marginal estimates ---
  conn_me <- .EstimateMarginalEffects(
    conn_fit, predictor_term, approach,
    estimand = estimand, from_var = from_var, to_var = to_var
  )
  dist_me <- .EstimateMarginalEffects(
    dist_fit, predictor_term, approach,
    estimand = estimand, from_var = from_var, to_var = to_var
  )

  list(
    connectivity     = conn_me,
    distance         = dist_me,
    models           = list(
      connectivity = conn_fit,
      distance     = dist_fit
    ),
    unpaired_summary = unpaired_summary,
    unpaired_rates   = unpaired_rates
  )
}

#' @title Forest Plot of TCR Dyad Model Estimates
#'
#' @description Generates a forest plot from the output of
#'   \code{\link{ModelTCRDyads}}. The function detects whether the result
#'   contains marginal effects (\code{estimand = "effects"}, columns from
#'   \code{\link[marginaleffects]{avg_comparisons}}) or marginal means
#'   (\code{estimand = "means"}, columns from
#'   \code{\link[marginaleffects]{avg_predictions}}) and adjusts the plot
#'   accordingly.
#'
#'   For marginal effects, contrast labels appear on the y-axis with a dashed
#'   reference line at zero. For marginal means, \code{pair_type} levels
#'   appear on the y-axis with no reference line.
#'
#'   When \code{which = "connectivity"}, axis labels are tailored to
#'   emphasize that the estimates represent probabilities. When
#'   \code{showUnpaired = TRUE} (and \code{which = "connectivity"}), a
#'   second facet is added showing the empirical probability that a clone
#'   in each group has no connections, with exact binomial 95\% confidence
#'   intervals computed from the full dyad table prior to subsetting to
#'   connected pairs.
#'
#' @param result The list returned by \code{\link{ModelTCRDyads}}.
#' @param which Character. Which model to plot: \code{"connectivity"}
#'   (default) or \code{"distance"}.
#' @param showUnpaired Logical. When \code{TRUE} and
#'   \code{which = "connectivity"}, adds a facet showing the per-group
#'   probability that a clone is unpaired (has no connections). Requires
#'   a categorical \code{groupColumn} in \code{ModelTCRDyads()}.
#'   Default \code{FALSE}.
#' @param title Character or \code{NULL}. Plot title. When \code{NULL}
#'   (default), an informative title is generated automatically.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{ModelTCRDyads}}
#' @importFrom ggplot2 ggplot aes geom_pointrange geom_vline labs theme_bw ggtitle facet_wrap
#' @export
PlotDyadEstimates <- function(result, which = "connectivity",
                              showUnpaired = FALSE, title = NULL) {
  if (!which %in% c("connectivity", "distance")) {
    stop("'which' must be 'connectivity' or 'distance', got: ", which)
  }
  df <- result[[which]]
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L) {
    stop("No estimates found in result$", which,
         ". Ensure ModelTCRDyads() ran successfully.")
  }

  is_effects <- "contrast" %in% colnames(df)
  is_conn    <- which == "connectivity"

  if (is_effects) {
    y_var  <- "contrast"
    x_lab  <- if (is_conn) "Difference in Connectivity Probability" else "Difference Estimate"
    y_lab  <- "Contrast"
    auto_title <- paste0(
      tools::toTitleCase(which), ": Pairwise Contrasts"
    )
    auto_sub <- if (is_conn) {
      "Estimated difference in probability (95% CI)"
    } else {
      "Estimated difference (95% CI)"
    }
  } else {
    y_var  <- "pair_type"
    x_lab  <- if (is_conn) "Predicted Connectivity Probability" else "Predicted Value"
    y_lab  <- "Pair Type"
    auto_title <- paste0(
      tools::toTitleCase(which), ": Marginal Means"
    )
    auto_sub <- if (is_conn) {
      "Predicted probability (95% CI)"
    } else {
      "Predicted value (95% CI)"
    }
  }

  if (!y_var %in% colnames(df)) {
    stop("Expected column '", y_var, "' not found in result$", which)
  }

  plot_title <- if (!is.null(title)) title else auto_title

  # determine whether to include the unpaired facet
  include_unpaired <- showUnpaired && is_conn &&
    !is.null(result$unpaired_rates) && nrow(result$unpaired_rates) > 0L

  if (showUnpaired && !is_conn) {
    warning("showUnpaired is ignored when which = 'distance'")
  }
  if (showUnpaired && is_conn && is.null(result$unpaired_rates)) {
    warning("No unpaired_rates found in result; was groupColumn numeric? ",
            "showUnpaired requires a categorical groupColumn.")
  }

  if (include_unpaired) {
    conn_panel <- data.frame(
      label     = as.character(df[[y_var]]),
      estimate  = df[["estimate"]],
      conf.low  = df[["conf.low"]],
      conf.high = df[["conf.high"]],
      panel     = "Connectivity Probability",
      stringsAsFactors = FALSE
    )
    unpaired_panel <- data.frame(
      label     = result$unpaired_rates$group,
      estimate  = result$unpaired_rates$estimate,
      conf.low  = result$unpaired_rates$conf.low,
      conf.high = result$unpaired_rates$conf.high,
      panel     = "Probability Unpaired",
      stringsAsFactors = FALSE
    )
    combined <- rbind(conn_panel, unpaired_panel)
    combined$panel <- factor(combined$panel,
                             levels = c("Connectivity Probability",
                                        "Probability Unpaired"))

    p <- ggplot2::ggplot(combined, ggplot2::aes(
      x    = .data[["estimate"]],
      y    = stats::reorder(.data[["label"]], .data[["estimate"]]),
      xmin = .data[["conf.low"]],
      xmax = .data[["conf.high"]]
    )) +
      ggplot2::geom_pointrange(size = 0.5) +
      ggplot2::facet_wrap(~ panel, scales = "free") +
      ggplot2::labs(
        title    = plot_title,
        subtitle = auto_sub,
        x        = "Probability",
        y        = NULL
      ) +
      ggplot2::theme_bw()

    if (is_effects) {
      ref_df <- data.frame(
        panel      = factor("Connectivity Probability",
                            levels = levels(combined$panel)),
        xintercept = 0
      )
      p <- p + ggplot2::geom_vline(
        data = ref_df,
        ggplot2::aes(xintercept = .data[["xintercept"]]),
        linetype = "dashed", color = "firebrick"
      )
    }
  } else {
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x    = .data[["estimate"]],
      y    = stats::reorder(.data[[y_var]], .data[["estimate"]]),
      xmin = .data[["conf.low"]],
      xmax = .data[["conf.high"]]
    )) +
      ggplot2::geom_pointrange(size = 0.5) +
      ggplot2::labs(
        title    = plot_title,
        subtitle = auto_sub,
        x        = x_lab,
        y        = y_lab
      ) +
      ggplot2::theme_bw()

    if (is_effects) {
      p <- p + ggplot2::geom_vline(
        xintercept = 0, linetype = "dashed", color = "firebrick"
      )
    }
  }

  p
}
