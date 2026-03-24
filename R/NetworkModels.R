#' @include Visualization.R
NULL

#' @title Model TCR Clone Connection Probability
#'
#' @description Fits a regression over all unique clone pairs derived from the
#'   TCR distance matrix. The binary response indicates whether a pair is
#'   connected (distance \eqn{\leq} \code{distanceThreshold}). For categorical
#'   \code{metadataVar}, a \code{same_<var>} indicator (1 if both endpoints
#'   share the same value) is used as the fixed effect; for numeric variables,
#'   \code{from_<var>} and \code{to_<var>} enter as separate additive fixed
#'   effects.
#'
#'   When \code{approach = "mixed"} (default), crossed random intercepts are
#'   included for the subject owning each endpoint via
#'   \code{\link[lme4]{glmer}} —
#'   \code{(1|from_<subjectIdCol>) + (1|to_<subjectIdCol>)}. When
#'   \code{approach = "clustered"}, a plain binomial \code{\link[stats]{glm}}
#'   is fitted and cluster-robust standard errors (clustered on subject-pair
#'   \code{Dyad_ID}) are computed via \code{\link[sandwich]{vcovCL}}; this
#'   avoids convergence issues at the cost of assuming independence across
#'   subject pairs.
#'
#'   When \code{communityMethod = "DIANA"} (default), each dyad is additionally
#'   annotated with \code{diana_community} (the DIANA cluster label for
#'   \code{from_clone} and \code{to_clone} encoded as \code{same_diana_community})
#'   so models can condition on whether two clones belong to the same pre-computed
#'   cluster. When \code{communityMethod = "threshold"}, a Louvain community over
#'   the thresholded graph is used instead.
#'
#' @param seuratObj_TCR A Seurat object with TCR distance matrices stored in
#'   \code{@misc$TCR_Distances} (output of \code{CalculateTcrDistances}).
#' @param chains Character vector of chain(s), e.g. \code{"TRB"} or
#'   \code{c("TRA","TRB")}. Must match the chains used in
#'   \code{CalculateTcrDistances}.
#' @param metadataVar Character. Seurat metadata column to test as the fixed
#'   effect, e.g. \code{"outcome"} or \code{"age"}.
#' @param distanceThreshold Numeric. Distance cutoff that defines connectivity;
#'   should match the value used in \code{TCRDistanceNetwork}.
#' @param cdr3Only Logical. Use CDR3-only distances. Default \code{FALSE}.
#' @param communityMethod Character. \code{"DIANA"} (default) reads
#'   pre-computed DIANA cluster assignments from \code{RunTcrClustering} and
#'   annotates dyads with \code{from_diana_community}, \code{to_diana_community},
#'   and \code{same_diana_community}. \code{"threshold"} derives communities via
#'   Louvain detection on the thresholded graph instead.
#' @param approach Character. \code{"mixed"} (default) fits a binomial GLMM
#'   via \code{\link[lme4]{glmer}} with crossed random intercepts for subject
#'   at each dyad endpoint. \code{"clustered"} fits a plain binomial
#'   \code{\link[stats]{glm}} and applies cluster-robust standard errors
#'   (clustered on \code{Dyad_ID}, a sorted subject-pair identifier) via
#'   \code{\link[sandwich]{vcovCL}}; use this when the mixed model fails to
#'   converge.
#' @param subjectIdCol Character. Metadata column containing subject IDs.
#'   Used for random intercepts (\code{approach = "mixed"}) or as the basis
#'   for \code{Dyad_ID} clustering (\code{approach = "clustered"}).
#'   Default \code{"SubjectId"}.
#' @param verbose Logical. Emit progress messages. Default \code{FALSE}.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{model}}{The fitted \code{glmerMod} (mixed) or \code{glm}
#'       (clustered) object.}
#'     \item{\code{data}}{The full dyad data frame (one row per unique clone
#'       pair) used for fitting. Columns include \code{distance},
#'       \code{connected} (0/1), \code{from_<metadataVar>},
#'       \code{to_<metadataVar>}, \code{same_<metadataVar>} (categorical
#'       only), and all \code{from_}/\code{to_}-prefixed Seurat metadata
#'       columns. When \code{approach = "clustered"}, a \code{Dyad_ID} column
#'       is also added.}
#'     \item{\code{coeftest}}{For \code{approach = "clustered"}, the
#'       \code{\link[lmtest]{coeftest}} result with cluster-robust standard
#'       errors. \code{NULL} when \code{approach = "mixed"}.}
#'   }
#'
#' @note The dyad table contains \eqn{N(N-1)/2} rows where \eqn{N} is the
#'   number of unique clones; this can be memory-intensive for large datasets.
#'   Mixed models may not converge with few subjects or highly imbalanced data;
#'   use \code{approach = "clustered"} as a fallback.
#'
#' @seealso \code{\link{TCRDistanceNetwork}}, \code{\link{ModelTCREdgeDistances}}
#' @importFrom lme4 glmer
#' @importFrom stats as.formula setNames glm
#' @importFrom igraph cluster_louvain membership V ecount vcount
#' @importFrom sandwich vcovCL
#' @importFrom lmtest coeftest
#' @export
ModelTCRConnectivity <- function(
    seuratObj_TCR,
    chains,
    metadataVar,
    distanceThreshold,
    cdr3Only        = FALSE,
    approach        = "mixed",
    communityMethod = "DIANA",
    subjectIdCol    = "SubjectId",
    verbose         = FALSE
) {
  if (!approach %in% c("mixed", "clustered")) {
    stop("approach must be 'mixed' or 'clustered', got: ", approach)
  }
  if (!communityMethod %in% c("DIANA", "threshold")) {
    stop("communityMethod must be 'DIANA' or 'threshold', got: ", communityMethod)
  }
  chainPrefix <- if (length(chains) > 1L) .get_chain_field_prefix(chains) else chains
  assayName   <- paste0(chainPrefix, "_", ifelse(cdr3Only, "cdr3", "fl"))
  dist_m      <- as.matrix(GetDistanceMatrix(seuratObj_TCR, chains, cdr3Only = cdr3Only))
  cloneIdxCol <- paste0(chainPrefix, "_CloneIdx")

  if (communityMethod == "DIANA") {
    clusterIdxCol <- paste0(assayName, "_ClusterIdx")
    if (!clusterIdxCol %in% colnames(seuratObj_TCR@meta.data)) {
      stop("communityMethod = 'DIANA' requires pre-computed cluster assignments. ",
           "Column '", clusterIdxCol, "' not found. Run RunTcrClustering() first.")
    }
  }

  # all unique clone pairs (upper triangle)
  idx     <- which(upper.tri(dist_m), arr.ind = TRUE)
  dyad_df <- data.frame(
    from_clone = rownames(dist_m)[idx[, 1L]],
    to_clone   = colnames(dist_m)[idx[, 2L]],
    distance   = dist_m[idx],
    stringsAsFactors = FALSE
  )
  dyad_df$connected <- as.integer(dyad_df$distance <= distanceThreshold)

  if (verbose) {
    message("Total dyads: ", nrow(dyad_df),
            " | connected: ", sum(dyad_df$connected))
  }

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
    dyad_df          <- dplyr::left_join(
      dyad_df, from_meta,
      by = stats::setNames(from_key, "from_clone")
    )
    dyad_df          <- dplyr::left_join(
      dyad_df, to_meta,
      by = stats::setNames(to_key, "to_clone")
    )
  }

  # annotate community membership
  if (communityMethod == "DIANA") {
    from_clust <- paste0("from_", clusterIdxCol)
    to_clust   <- paste0("to_",   clusterIdxCol)
    # treat cluster 0 (noise) as NA
    dyad_df[[from_clust]][!is.na(dyad_df[[from_clust]]) & dyad_df[[from_clust]] == "0"] <- NA_character_
    dyad_df[[to_clust]][!is.na(dyad_df[[to_clust]])   & dyad_df[[to_clust]]   == "0"] <- NA_character_
    dyad_df$from_diana_community <- dyad_df[[from_clust]]
    dyad_df$to_diana_community   <- dyad_df[[to_clust]]
    dyad_df$same_diana_community <- as.integer(
      !is.na(dyad_df$from_diana_community) & !is.na(dyad_df$to_diana_community) &
        dyad_df$from_diana_community == dyad_df$to_diana_community
    )
  } else {
    # threshold community via Louvain on thresholded graph
    tg <- BuildTCRDistanceGraph(
      seuratObj_TCR, chains = chains, cdr3Only = cdr3Only,
      distanceThreshold = distanceThreshold, edgeType = "binary", verbose = FALSE
    )
    g_full <- tidygraph::as.igraph(tg)
    if (igraph::ecount(g_full) > 0L) {
      set.seed(42L)
      comms        <- igraph::cluster_louvain(g_full)
      comm_vec     <- stats::setNames(as.character(igraph::membership(comms)),
                                      igraph::V(g_full)$name)
    } else {
      comm_vec <- stats::setNames(rep(NA_character_, igraph::vcount(g_full)),
                                  igraph::V(g_full)$name)
    }
    dyad_df$from_community      <- comm_vec[dyad_df$from_clone]
    dyad_df$to_community        <- comm_vec[dyad_df$to_clone]
    dyad_df$same_community      <- as.integer(
      !is.na(dyad_df$from_community) & !is.na(dyad_df$to_community) &
        dyad_df$from_community == dyad_df$to_community
    )
  }

  from_var <- paste0("from_", metadataVar)
  to_var   <- paste0("to_",   metadataVar)
  if (!from_var %in% names(dyad_df)) {
    stop("metadataVar '", metadataVar, "' not found in Seurat metadata")
  }

  # fixed-effect term: same-group indicator for categorical, raw values for numeric
  if (is.numeric(dyad_df[[from_var]])) {
    predictor_term <- paste(from_var, to_var, sep = " + ")
  } else {
    same_col            <- paste0("same_", metadataVar)
    dyad_df[[same_col]] <- as.integer(
      !is.na(dyad_df[[from_var]]) & !is.na(dyad_df[[to_var]]) &
        dyad_df[[from_var]] == dyad_df[[to_var]]
    )
    predictor_term <- same_col
  }

  # subject columns (used for random intercepts or Dyad_ID clustering)
  from_subj <- paste0("from_", subjectIdCol)
  to_subj   <- paste0("to_",   subjectIdCol)
  if (!from_subj %in% names(dyad_df)) {
    stop("subjectIdCol '", subjectIdCol, "' not found in Seurat metadata; ",
         "ensure the column is present before fitting connectivity models.")
  }

  if (approach == "mixed") {
    rand_term <- paste0("(1|", from_subj, ") + (1|", to_subj, ")")
    fml_str   <- paste0("connected ~ ", predictor_term, " + ", rand_term)
    if (verbose) message("Fitting mixed model: ", fml_str)
    fit      <- lme4::glmer(stats::as.formula(fml_str), data = dyad_df, family = "binomial")
    coef_out <- NULL
  } else {
    # sorted subject-pair identifier — groups all clone dyads between the same two subjects
    dyad_df$Dyad_ID <- paste0(
      pmin(dyad_df[[from_subj]], dyad_df[[to_subj]]), "__",
      pmax(dyad_df[[from_subj]], dyad_df[[to_subj]])
    )
    fml_str <- paste0("connected ~ ", predictor_term)
    if (verbose) message("Fitting with clustered SEs: ", fml_str)
    fit      <- stats::glm(stats::as.formula(fml_str), data = dyad_df, family = "binomial")
    vcov_cl  <- sandwich::vcovCL(fit, cluster = dyad_df$Dyad_ID)
    coef_out <- lmtest::coeftest(fit, vcov = vcov_cl)
  }

  list(model = fit, data = dyad_df, coeftest = coef_out)
}


#' @title Model Distances Within TCR Network Components
#'
#' @description Fits a Poisson regression over connected clone pairs in the TCR
#'   distance network. The integer response is the raw tcrdist distance. For
#'   categorical \code{metadataVar}, a \code{same_<var>} indicator is derived;
#'   for numeric variables, \code{from_<var>} and \code{to_<var>} enter as
#'   separate additive fixed effects.
#'
#'   When \code{approach = "mixed"} (default), crossed random intercepts are
#'   included for the subject owning each endpoint via
#'   \code{\link[lme4]{glmer}} —
#'   \code{(1|from_<subjectIdCol>) + (1|to_<subjectIdCol>)}. When
#'   \code{approach = "clustered"}, a plain Poisson \code{\link[stats]{glm}}
#'   is fitted and cluster-robust standard errors (clustered on subject-pair
#'   \code{Dyad_ID}) are computed via \code{\link[sandwich]{vcovCL}}.
#'
#'   The primary path builds the edge table directly from \code{seuratObj_TCR},
#'   \code{chains}, and \code{distanceThreshold}. When
#'   \code{communityMethod = "DIANA"} (default), only edges where both endpoints
#'   share the same DIANA cluster are retained before fitting, ensuring the
#'   distance model is conditioned on membership in the same pre-computed cluster.
#'   \code{communityMethod = "threshold"} retains all thresholded edges regardless
#'   of community. Supplying a pre-computed \code{edges} data frame skips graph
#'   construction entirely as a time-saving backup.
#'
#' @param seuratObj_TCR A Seurat object with TCR distance matrices stored in
#'   \code{@misc$TCR_Distances} (output of \code{CalculateTcrDistances}).
#'   Required unless \code{edges} is supplied.
#' @param chains Character vector of chain(s), e.g. \code{"TRB"} or
#'   \code{c("TRA","TRB")}. Required unless \code{edges} is supplied.
#' @param metadataVar Character. Name of the metadata variable to test as a
#'   fixed effect, e.g. \code{"outcome"}.
#' @param distanceThreshold Numeric. Maximum pairwise distance for an edge.
#'   Required unless \code{edges} is supplied.
#' @param cdr3Only Logical. Use CDR3-only distances. Default \code{FALSE}.
#'   Ignored when \code{edges} is supplied.
#' @param communityMethod Character. \code{"DIANA"} (default) filters the edge
#'   table to retain only intra-cluster edges (both endpoints sharing the same
#'   DIANA cluster from \code{RunTcrClustering}). \code{"threshold"} retains all
#'   thresholded edges. Ignored when \code{edges} is supplied.
#' @param edges Data frame or \code{NULL}. When non-\code{NULL}, use this
#'   pre-computed edge table (e.g. \code{result$edges} from
#'   \code{TCRDistanceNetwork}) instead of building from scratch. All of
#'   \code{seuratObj_TCR}, \code{chains}, \code{distanceThreshold}, and
#'   \code{cdr3Only} are then ignored.
#' @param components Integer vector or \code{NULL}. When supplied, only edges
#'   whose \code{component} value is in this vector are retained before fitting.
#'   Default \code{NULL} retains all edges.
#' @param approach Character. \code{"mixed"} (default) fits a Poisson GLMM
#'   via \code{\link[lme4]{glmer}} with crossed random intercepts for subject
#'   at each edge endpoint. \code{"clustered"} fits a plain Poisson
#'   \code{\link[stats]{glm}} and applies cluster-robust standard errors
#'   (clustered on \code{Dyad_ID}, a sorted subject-pair identifier) via
#'   \code{\link[sandwich]{vcovCL}}; use this when the mixed model fails to
#'   converge.
#' @param subjectIdCol Character. Subject ID column name without the
#'   \code{from_}/\code{to_} prefix. Used for random intercepts
#'   (\code{approach = "mixed"}) or as the basis for \code{Dyad_ID} clustering
#'   (\code{approach = "clustered"}). Default \code{"SubjectId"}.
#' @param verbose Logical. Emit progress messages. Default \code{FALSE}.
#'
#' @return A named list:
#'   \describe{
#'     \item{\code{model}}{The fitted \code{glmerMod} (mixed) or \code{glm}
#'       (clustered) object.}
#'     \item{\code{data}}{The edge data frame used for fitting (after optional
#'       component filtering and predictor construction). When
#'       \code{approach = "clustered"}, a \code{Dyad_ID} column is also added.}
#'     \item{\code{coeftest}}{For \code{approach = "clustered"}, the
#'       \code{\link[lmtest]{coeftest}} result with cluster-robust standard
#'       errors. \code{NULL} when \code{approach = "mixed"}.}
#'   }
#'
#' @note Mixed models may not converge with very few edges or subjects. Use the
#'   \code{components} argument to constrain the analysis to specific connected
#'   subgraphs, or use \code{approach = "clustered"} to avoid convergence issues.
#'
#' @seealso \code{\link{TCRDistanceNetwork}}, \code{\link{ModelTCRConnectivity}}
#' @importFrom lme4 glmer
#' @importFrom stats as.formula glm
#' @importFrom igraph components
#' @importFrom tidygraph activate
#' @importFrom sandwich vcovCL
#' @importFrom lmtest coeftest
#' @export
ModelTCREdgeDistances <- function(
    seuratObj_TCR     = NULL,
    chains            = NULL,
    metadataVar,
    distanceThreshold = NULL,
    cdr3Only          = FALSE,
    approach          = "mixed",
    communityMethod   = "DIANA",
    edges             = NULL,
    components        = NULL,
    subjectIdCol      = "SubjectId",
    verbose           = FALSE
) {
  if (!approach %in% c("mixed", "clustered")) {
    stop("approach must be 'mixed' or 'clustered', got: ", approach)
  }
  if (!communityMethod %in% c("DIANA", "threshold")) {
    stop("communityMethod must be 'DIANA' or 'threshold', got: ", communityMethod)
  }
  has_primary <- !is.null(seuratObj_TCR) && !is.null(chains) && !is.null(distanceThreshold)
  has_backup  <- !is.null(edges)

  if (!has_primary && !has_backup) {
    stop("Provide 'seuratObj_TCR', 'chains', and 'distanceThreshold' (primary path) ",
         "or a pre-computed 'edges' data frame from TCRDistanceNetwork() (backup path).")
  }

  # primary path: build edge table from scratch via BuildTCRDistanceGraph
  if (is.null(edges)) {
    chainPrefix <- if (length(chains) > 1L) .get_chain_field_prefix(chains) else chains
    assayName   <- paste0(chainPrefix, "_", ifelse(cdr3Only, "cdr3", "fl"))

    if (communityMethod == "DIANA") {
      clusterIdxCol <- paste0(assayName, "_ClusterIdx")
      if (!clusterIdxCol %in% colnames(seuratObj_TCR@meta.data)) {
        stop("communityMethod = 'DIANA' requires pre-computed cluster assignments. ",
             "Column '", clusterIdxCol, "' not found. Run RunTcrClustering() first.")
      }
    }

    tcrGraph  <- BuildTCRDistanceGraph(
      seuratObj_TCR,
      chains            = chains,
      cdr3Only          = cdr3Only,
      distanceThreshold = distanceThreshold,
      edgeType          = "binary",
      verbose           = verbose
    )
    edge_data <- as.data.frame(tidygraph::activate(tcrGraph, edges))
    node_data <- as.data.frame(tidygraph::activate(tcrGraph, nodes))
    comps     <- igraph::components(tidygraph::as.igraph(tcrGraph))

    names(node_data)[names(node_data) == "name"] <- "clone"
    from_meta        <- node_data
    to_meta          <- node_data
    names(from_meta) <- paste0("from_", names(from_meta))
    names(to_meta)   <- paste0("to_",   names(to_meta))

    edge_cols <- setdiff(names(edge_data), c("from", "to"))
    if (nrow(edge_data) > 0L) {
      edges <- cbind(
        edge_data[edge_cols],
        component = as.integer(comps$membership)[edge_data$from],
        from_meta[edge_data$from, , drop = FALSE],
        to_meta  [edge_data$to,   , drop = FALSE]
      )
      rownames(edges) <- NULL
    } else {
      edges <- cbind(
        edge_data[edge_cols],
        component = integer(0L),
        from_meta[integer(0L), , drop = FALSE],
        to_meta  [integer(0L), , drop = FALSE]
      )
    }

    # DIANA: retain only intra-cluster edges
    if (communityMethod == "DIANA" && nrow(edges) > 0L) {
      from_clust <- paste0("from_", clusterIdxCol)
      to_clust   <- paste0("to_",   clusterIdxCol)
      edges[[from_clust]][!is.na(edges[[from_clust]]) & edges[[from_clust]] == "0"] <- NA_character_
      edges[[to_clust]][!is.na(edges[[to_clust]])   & edges[[to_clust]]   == "0"] <- NA_character_
      same_cluster <- !is.na(edges[[from_clust]]) & !is.na(edges[[to_clust]]) &
        edges[[from_clust]] == edges[[to_clust]]
      edges <- edges[same_cluster, , drop = FALSE]
      if (verbose) message("DIANA intra-cluster edges retained: ", nrow(edges))
    }
  }

  # optional component filtering (applies to both paths)
  if (!is.null(components)) {
    edges <- edges[!is.na(edges$component) & edges$component %in% components, , drop = FALSE]
  }

  if (nrow(edges) == 0L) {
    stop("No edges remain after component filtering.")
  }

  if (verbose) message("Edges used for modelling: ", nrow(edges))

  from_var <- paste0("from_", metadataVar)
  to_var   <- paste0("to_",   metadataVar)
  if (!from_var %in% names(edges)) {
    stop("metadataVar '", metadataVar, "' not found in edges; ",
         "ensure it was a metadata column when building the network.")
  }

  # fixed-effect term
  if (is.numeric(edges[[from_var]])) {
    predictor_term <- paste(from_var, to_var, sep = " + ")
  } else {
    same_col          <- paste0("same_", metadataVar)
    edges[[same_col]] <- as.integer(
      !is.na(edges[[from_var]]) & !is.na(edges[[to_var]]) &
        edges[[from_var]] == edges[[to_var]]
    )
    predictor_term <- same_col
  }

  # subject columns (used for random intercepts or Dyad_ID clustering)
  from_subj <- paste0("from_", subjectIdCol)
  to_subj   <- paste0("to_",   subjectIdCol)
  if (!from_subj %in% names(edges)) {
    stop("subjectIdCol '", subjectIdCol, "' not found in edges; ",
         "ensure it was a metadata column when building the network.")
  }

  if (approach == "mixed") {
    rand_term <- paste0("(1|", from_subj, ") + (1|", to_subj, ")")
    fml_str   <- paste0("distance ~ ", predictor_term, " + ", rand_term)
    if (verbose) message("Fitting mixed model: ", fml_str)
    fit      <- lme4::glmer(stats::as.formula(fml_str), data = edges, family = "poisson")
    coef_out <- NULL
  } else {
    # sorted subject-pair identifier — groups all clone dyads between the same two subjects
    edges$Dyad_ID <- paste0(
      pmin(edges[[from_subj]], edges[[to_subj]]), "__",
      pmax(edges[[from_subj]], edges[[to_subj]])
    )
    fml_str <- paste0("distance ~ ", predictor_term)
    if (verbose) message("Fitting with clustered SEs: ", fml_str)
    fit      <- stats::glm(stats::as.formula(fml_str), data = edges, family = "poisson")
    vcov_cl  <- sandwich::vcovCL(fit, cluster = edges$Dyad_ID)
    coef_out <- lmtest::coeftest(fit, vcov = vcov_cl)
  }

  list(model = fit, data = edges, coeftest = coef_out)
}
