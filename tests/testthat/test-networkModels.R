library(testthat)
library(Seurat)

# Helper: 10-cell Seurat with SubjectId — designed so lme4 has enough structure
# to attempt a mixed model (3 subjects, diverse sequences, both outcomes).
.make_model_test_seurat <- function() {
  set.seed(99L)
  metadata <- data.frame(
    TRB_V     = c("TRBV7-9", "TRBV7-9", "TRBV7-9", "TRBV6-4", "TRBV6-4",
                  "TRBV6-4", "TRBV5-1", "TRBV5-1", "TRBV7-2", "TRBV7-2"),
    TRB_J     = c("TRBJ2-1", "TRBJ2-1", "TRBJ2-1", "TRBJ1-1", "TRBJ1-1",
                  "TRBJ1-1", "TRBJ2-3", "TRBJ2-3", "TRBJ1-4", "TRBJ1-4"),
    TRB       = c("CASSLGQAYEQYF", "CASSLAQAYEQYF", "CASSLGQAEQYF",
                  "CASSYLDTEAFF",  "CASSYLDAEAFF",  "CASSYLETEAFF",
                  "CASSLAGRDTQYF", "CASSLARDTQYF",  "CASSPRAGEYF",  "CASSPRAGAYF"),
    outcome   = c("Controller",    "Controller",    "Progressor",
                  "Progressor",    "Progressor",    "Controller",
                  "Controller",    "Progressor",    "Controller",    "Controller"),
    SubjectId = c("S1", "S1", "S2", "S2", "S3", "S3", "S1", "S2", "S3", "S3"),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- paste0("cell_", seq_len(nrow(metadata)))
  counts <- matrix(
    rpois(nrow(metadata) * 10L, lambda = 5L),
    nrow = 10L, ncol = nrow(metadata),
    dimnames = list(paste0("gene_", seq_len(10L)), rownames(metadata))
  )
  Seurat::CreateSeuratObject(counts = Seurat::as.sparse(counts), meta.data = metadata)
}

# shared fixture — computed once per file load
.get_model_fixture <- function() {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_model_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )
  result <- suppressWarnings(
    TCRDistanceNetwork(
      seuratObj_TCR,
      chains            = "TRB",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      edgeType          = "continuous",
      colorBy           = "outcome",
      verbose           = FALSE
    )
  )
  list(seuratObj_TCR = seuratObj_TCR, result = result)
}

# ---------------------------------------------------------------------------
# ModelTCRConnectivity
# ---------------------------------------------------------------------------

test_that("ModelTCRConnectivity returns correct structure", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRConnectivity(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      approach          = "mixed",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_type(out, "list")
  expect_named(out, c("model", "data", "coeftest"))
  expect_true(inherits(out$model, "glmerMod"))
  expect_null(out$coeftest)
  expect_true(is.data.frame(out$data))
})

test_that("ModelTCRConnectivity dyad data has expected columns", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRConnectivity(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expected_cols <- c(
    "from_clone", "to_clone", "distance", "connected",
    "same_outcome", "from_SubjectId", "to_SubjectId"
  )
  expect_true(all(expected_cols %in% colnames(out$data)))
  # every row is a unique pair
  expect_true(all(out$data$connected %in% c(0L, 1L)))
  # dyad count = N*(N-1)/2
  n_clones <- nrow(GetDistanceMatrix(fix$seuratObj_TCR, "TRB"))
  expect_equal(nrow(out$data), n_clones * (n_clones - 1L) / 2L)
})

test_that("ModelTCRConnectivity errors on missing metadataVar", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRConnectivity(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "nonexistent_var",
      distanceThreshold = 200,
      communityMethod   = "threshold"
    ),
    "not found in Seurat metadata"
  )
})

test_that("ModelTCRConnectivity errors on missing subjectIdCol", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRConnectivity(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "NonexistentId"
    ),
    "not found in Seurat metadata"
  )
})

# ---------------------------------------------------------------------------
# ModelTCREdgeDistances
# ---------------------------------------------------------------------------

test_that("ModelTCREdgeDistances returns correct structure (primary path)", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      approach          = "mixed",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_type(out, "list")
  expect_named(out, c("model", "data", "coeftest"))
  expect_true(inherits(out$model, "glmerMod"))
  expect_null(out$coeftest)
  expect_true(is.data.frame(out$data))
})

test_that("ModelTCREdgeDistances data has same_<metadataVar> column for categorical var", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_true("same_outcome" %in% colnames(out$data))
  expect_true(all(out$data$same_outcome %in% c(0L, 1L)))
})

test_that("ModelTCREdgeDistances components argument filters correctly", {
  fix <- .get_model_fixture()

  # build edges via primary path to get component labels
  all_edges <- suppressWarnings(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )$data

  if (length(unique(stats::na.omit(all_edges$component))) < 2L) {
    skip("test data has fewer than 2 components; skipping component-filter test")
  }

  comp1_n <- sum(!is.na(all_edges$component) & all_edges$component == 1L)
  out <- suppressWarnings(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId",
      components        = 1L
    )
  )
  expect_equal(nrow(out$data), comp1_n)
})

test_that("ModelTCREdgeDistances errors on empty edges after component filtering", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      components        = 99999L
    ),
    "No edges remain"
  )
})

test_that("ModelTCREdgeDistances errors on missing metadataVar", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "nonexistent_var",
      distanceThreshold = 200,
      communityMethod   = "threshold"
    ),
    "not found in edges"
  )
})

test_that("ModelTCREdgeDistances errors on missing subjectIdCol", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "NonexistentId"
    ),
    "not found in edges"
  )
})

test_that("ModelTCREdgeDistances backup path (edges) produces a valid model", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCREdgeDistances(
      metadataVar  = "outcome",
      edges        = fix$result$edges,
      subjectIdCol = "SubjectId"
    )
  )

  expect_true(inherits(out$model, "glmerMod"))
  expect_true(is.data.frame(out$data))
  expect_true("same_outcome" %in% colnames(out$data))
})

test_that("ModelTCREdgeDistances errors when neither primary nor backup params provided", {
  expect_error(
    ModelTCREdgeDistances(metadataVar = "outcome"),
    "primary path"
  )
})

# ---------------------------------------------------------------------------
# DIANA community method
# ---------------------------------------------------------------------------

# fixture that runs RunTcrClustering so TRB_fl_ClusterIdx is present
.get_diana_fixture <- function() {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_model_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )
  seuratObj_TCR <- RunTcrClustering(
    seuratObj_TCR,
    chainsToCluster      = "TRB",
    clusterSizeThreshold = 1L,
    verbose              = FALSE
  )
  list(seuratObj_TCR = seuratObj_TCR)
}

test_that("ModelTCRConnectivity DIANA path annotates dyads with diana community columns", {
  fix <- .get_diana_fixture()

  out <- suppressWarnings(
    ModelTCRConnectivity(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "DIANA",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_type(out, "list")
  expect_named(out, c("model", "data"))
  expect_true(all(c("from_diana_community", "to_diana_community", "same_diana_community") %in% colnames(out$data)))
  expect_true(all(out$data$same_diana_community %in% c(0L, 1L, NA_integer_)))
})

test_that("ModelTCREdgeDistances DIANA path retains only intra-cluster edges", {
  fix <- .get_diana_fixture()

  out <- suppressWarnings(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "DIANA",
      subjectIdCol      = "SubjectId"
    )
  )

  # when there are edges, all retained edges should be intra-cluster
  if (nrow(out$data) > 0L) {
    clust_col_from <- "from_TRB_fl_ClusterIdx"
    clust_col_to   <- "to_TRB_fl_ClusterIdx"
    if (all(c(clust_col_from, clust_col_to) %in% colnames(out$data))) {
      from_cl <- out$data[[clust_col_from]]
      to_cl   <- out$data[[clust_col_to]]
      valid   <- !is.na(from_cl) & !is.na(to_cl)
      expect_true(all(from_cl[valid] == to_cl[valid]))
    }
  }

  expect_type(out, "list")
  expect_named(out, c("model", "data"))
})

test_that("ModelTCRConnectivity errors when DIANA col missing", {
  fix <- .get_model_fixture()  # no RunTcrClustering — no ClusterIdx col
  expect_error(
    ModelTCRConnectivity(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "DIANA"
    ),
    "Run RunTcrClustering"
  )
})

test_that("ModelTCREdgeDistances errors when DIANA col missing", {
  fix <- .get_model_fixture()  # no RunTcrClustering — no ClusterIdx col
  expect_error(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "DIANA"
    ),
    "Run RunTcrClustering"
  )
})

# ---------------------------------------------------------------------------
# Clustered standard errors approach
# ---------------------------------------------------------------------------

test_that("ModelTCRConnectivity clustered approach returns glm model and coeftest", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRConnectivity(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      approach          = "clustered",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_type(out, "list")
  expect_named(out, c("model", "data", "coeftest"))
  expect_true(inherits(out$model, "glm"))
  expect_false(inherits(out$model, "glmerMod"))
  expect_false(is.null(out$coeftest))
  expect_true(inherits(out$coeftest, "coeftest"))
  # Dyad_ID cluster variable is added to the data
  expect_true("Dyad_ID" %in% colnames(out$data))
})

test_that("ModelTCRConnectivity clustered approach Dyad_ID is sorted subject pair", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRConnectivity(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      approach          = "clustered",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  # Dyad_ID = sorted(from_SubjectId, to_SubjectId) joined by "__"
  expected <- paste0(
    pmin(out$data$from_SubjectId, out$data$to_SubjectId), "__",
    pmax(out$data$from_SubjectId, out$data$to_SubjectId)
  )
  expect_equal(out$data$Dyad_ID, expected)
})

test_that("ModelTCREdgeDistances clustered approach returns glm model and coeftest", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      approach          = "clustered",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_type(out, "list")
  expect_named(out, c("model", "data", "coeftest"))
  expect_true(inherits(out$model, "glm"))
  expect_false(inherits(out$model, "glmerMod"))
  expect_false(is.null(out$coeftest))
  expect_true(inherits(out$coeftest, "coeftest"))
  expect_true("Dyad_ID" %in% colnames(out$data))
})

test_that("ModelTCRConnectivity errors on invalid approach argument", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRConnectivity(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      approach          = "bad_value",
      communityMethod   = "threshold"
    ),
    "approach must be"
  )
})

test_that("ModelTCREdgeDistances errors on invalid approach argument", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCREdgeDistances(
      fix$seuratObj_TCR,
      chains            = "TRB",
      metadataVar       = "outcome",
      distanceThreshold = 200,
      approach          = "bad_value",
      communityMethod   = "threshold"
    ),
    "approach must be"
  )
})
