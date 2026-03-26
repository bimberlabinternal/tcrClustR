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

# shared fixture: only CalculateTcrDistances needed (DIANA runs inline)
.get_model_fixture <- function() {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_model_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )
  list(seuratObj_TCR = seuratObj_TCR)
}

# ---------------------------------------------------------------------------
# ModelTCRDyads — return structure
# ---------------------------------------------------------------------------

test_that("ModelTCRDyads returns correct structure (threshold, mixed)", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      approach          = "mixed",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_type(out, "list")
  expect_named(out, c("connectivity", "distance", "models",
                      "unpaired_summary", "unpaired_rates"))

  # marginal estimate data frames (default estimand = "effects" -> avg_comparisons)
  expect_true(is.data.frame(out$connectivity))
  expect_true(is.data.frame(out$distance))
  expect_true("estimate" %in% colnames(out$connectivity))
  expect_true("estimate" %in% colnames(out$distance))

  # models sub-list
  expect_named(out$models, c("connectivity", "distance"))
  expect_true(inherits(out$models$connectivity$model, "glmerMod"))
  expect_true(inherits(out$models$distance$model, "glmerMod"))
  expect_true(is.data.frame(out$models$connectivity$data))
  expect_true(is.data.frame(out$models$distance$data))
})

test_that("ModelTCRDyads connectivity data has expected columns", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  conn_data <- out$models$connectivity$data
  expected_cols <- c(
    "from_clone", "to_clone", "distance", "connected",
    "pair_type", "from_SubjectId", "to_SubjectId"
  )
  expect_true(all(expected_cols %in% colnames(conn_data)))
  expect_true(all(conn_data$connected %in% c(0L, 1L)))

  # dyad count = N*(N-1)/2
  n_clones <- nrow(GetDistanceMatrix(fix$seuratObj_TCR, "TRB"))
  expect_equal(nrow(conn_data), n_clones * (n_clones - 1L) / 2L)
})

test_that("ModelTCRDyads unpaired_summary is a data frame for categorical var", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_true(is.data.frame(out$unpaired_summary))
  expect_true(all(c("clone", "outcome", "is_connected", "partner_groups") %in%
                    colnames(out$unpaired_summary)))
})

# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------

test_that("ModelTCRDyads errors on missing groupColumn", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "nonexistent_var",
      distanceThreshold = 200,
      communityMethod   = "threshold"
    ),
    "not found in Seurat metadata"
  )
})

test_that("ModelTCRDyads errors on missing subjectIdCol", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "NonexistentId"
    ),
    "not found in Seurat metadata"
  )
})

test_that("ModelTCRDyads errors on invalid approach argument", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      approach          = "bad_value",
      communityMethod   = "threshold"
    ),
    "approach must be"
  )
})

test_that("ModelTCRDyads errors on invalid communityMethod argument", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "bad_value"
    ),
    "communityMethod must be"
  )
})

test_that("ModelTCRDyads errors on invalid estimand argument", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      estimand          = "bad_value"
    ),
    "estimand must be"
  )
})

# ---------------------------------------------------------------------------
# DIANA community method
# ---------------------------------------------------------------------------

test_that("ModelTCRDyads DIANA path annotates dyads with diana community columns", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains               = "TRB",
      groupColumn          = "outcome",
      distanceThreshold    = 200,
      communityMethod      = "DIANA",
      clusterSizeThreshold = 1L,
      subjectIdCol         = "SubjectId"
    )
  )

  conn_data <- out$models$connectivity$data
  expect_true(all(c("from_diana_community", "to_diana_community",
                     "same_diana_community") %in% colnames(conn_data)))
  expect_true(all(conn_data$same_diana_community %in% c(0L, 1L, NA_integer_)))

  # for DIANA, connected == same_diana_community (cluster co-membership)
  valid <- !is.na(conn_data$same_diana_community)
  expect_equal(conn_data$connected[valid], conn_data$same_diana_community[valid])
})

test_that("ModelTCRDyads DIANA connected edges are all intra-cluster", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains               = "TRB",
      groupColumn          = "outcome",
      distanceThreshold    = 200,
      communityMethod      = "DIANA",
      clusterSizeThreshold = 1,
      subjectIdCol         = "SubjectId"
    )
  )

  # distance model data == connected edges, which for DIANA == intra-cluster
  dist_data <- out$models$distance$data
  if (nrow(dist_data) > 0L) {
    expect_true(all(!is.na(dist_data$from_diana_community)))
    expect_true(all(!is.na(dist_data$to_diana_community)))
    expect_true(all(dist_data$from_diana_community == dist_data$to_diana_community))
  }
})

test_that("ModelTCRDyads DIANA does not require RunTcrClustering", {
  # base fixture has no RunTcrClustering - DIANA clustering should run
  fix <- .get_model_fixture()
  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains               = "TRB",
      groupColumn          = "outcome",
      distanceThreshold    = 200,
      communityMethod      = "DIANA",
      clusterSizeThreshold = 1L,
      subjectIdCol         = "SubjectId"
    )
  )
  expect_type(out, "list")
  expect_named(out, c("connectivity", "distance", "models",
                      "unpaired_summary", "unpaired_rates"))
})

test_that("ModelTCRDyads DIANA errors when no clusters form", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains               = "TRB",
      groupColumn          = "outcome",
      distanceThreshold    = 0.001,
      communityMethod      = "DIANA",
      clusterSizeThreshold = 100L,
      subjectIdCol         = "SubjectId"
    ),
    "No DIANA clusters formed"
  )
})

test_that("ModelTCRDyads threshold errors when no connections exist", {
  fix <- .get_model_fixture()
  expect_error(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 0.001,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    ),
    "No connected pairs"
  )
})

# ---------------------------------------------------------------------------
# Clustered standard errors approach
# ---------------------------------------------------------------------------

test_that("ModelTCRDyads clustered approach returns glm models", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      approach          = "clustered",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_true(inherits(out$models$connectivity$model, "glm"))
  expect_false(inherits(out$models$connectivity$model, "glmerMod"))
  expect_true(inherits(out$models$distance$model, "glm"))
  expect_false(inherits(out$models$distance$model, "glmerMod"))

  # Dyad_ID cluster variable is added to the data
  expect_true("Dyad_ID" %in% colnames(out$models$connectivity$data))
  expect_true("Dyad_ID" %in% colnames(out$models$distance$data))
})

test_that("ModelTCRDyads clustered approach Dyad_ID is sorted subject pair", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      approach          = "clustered",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  conn_data <- out$models$connectivity$data
  expected <- paste0(
    pmin(conn_data$from_SubjectId, conn_data$to_SubjectId), "__",
    pmax(conn_data$from_SubjectId, conn_data$to_SubjectId)
  )
  expect_equal(conn_data$Dyad_ID, expected)
})

test_that("ModelTCRDyads clustered marginal estimates have std.error", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      approach          = "clustered",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_true("std.error" %in% colnames(out$connectivity))
  expect_true("std.error" %in% colnames(out$distance))
})

# ---------------------------------------------------------------------------
# communityDistanceThreshold (separate edge vs community threshold)
# ---------------------------------------------------------------------------

test_that("ModelTCRDyads accepts communityDistanceThreshold without error", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains                     = "TRB",
      groupColumn                = "outcome",
      distanceThreshold          = 200,
      communityDistanceThreshold = 50,
      approach                   = "clustered",
      communityMethod            = "threshold",
      subjectIdCol               = "SubjectId"
    )
  )

  expect_type(out, "list")
  expect_named(out, c("connectivity", "distance", "models",
                      "unpaired_summary", "unpaired_rates"))
  expect_true(is.data.frame(out$connectivity))
  expect_true("same_community" %in% colnames(out$models$connectivity$data))
})

# ---------------------------------------------------------------------------
# estimand: effects vs means
# ---------------------------------------------------------------------------

test_that("ModelTCRDyads estimand='means' returns pair_type column", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      estimand          = "means",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_true(is.data.frame(out$connectivity))
  expect_true("pair_type" %in% colnames(out$connectivity))
  expect_true("estimate" %in% colnames(out$connectivity))
  expect_true("pair_type" %in% colnames(out$distance))
})

# ---------------------------------------------------------------------------
# PlotDyadEstimates
# ---------------------------------------------------------------------------

test_that("PlotDyadEstimates returns ggplot for effects estimand", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      estimand          = "effects",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  p <- PlotDyadEstimates(out, which = "connectivity")
  expect_true(inherits(p, "ggplot"))

  p2 <- PlotDyadEstimates(out, which = "distance")
  expect_true(inherits(p2, "ggplot"))
})

test_that("PlotDyadEstimates returns ggplot for means estimand", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      estimand          = "means",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  p <- PlotDyadEstimates(out, which = "connectivity")
  expect_true(inherits(p, "ggplot"))
})

test_that("PlotDyadEstimates accepts custom title", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  p <- PlotDyadEstimates(out, which = "connectivity", title = "Custom Title")
  expect_true(inherits(p, "ggplot"))
  expect_equal(p$labels$title, "Custom Title")
})

test_that("PlotDyadEstimates errors on invalid which", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_error(PlotDyadEstimates(out, which = "bad"), "'which' must be")
})

# ---------------------------------------------------------------------------
# unpaired_rates and showUnpaired
# ---------------------------------------------------------------------------

test_that("ModelTCRDyads returns unpaired_rates with correct structure", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_true(is.data.frame(out$unpaired_rates))
  expect_true(all(c("group", "estimate", "conf.low", "conf.high",
                    "n_clones", "n_unpaired") %in% colnames(out$unpaired_rates)))
  expect_true(all(out$unpaired_rates$estimate >= 0 & out$unpaired_rates$estimate <= 1))
  expect_true(all(out$unpaired_rates$conf.low <= out$unpaired_rates$estimate))
  expect_true(all(out$unpaired_rates$conf.high >= out$unpaired_rates$estimate))
})

test_that("PlotDyadEstimates showUnpaired returns faceted ggplot", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      estimand          = "means",
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  p <- PlotDyadEstimates(out, which = "connectivity", showUnpaired = TRUE)
  expect_true(inherits(p, "ggplot"))
  # facet_wrap produces a FacetWrap object
  expect_true(inherits(p$facet, "FacetWrap"))
})

test_that("PlotDyadEstimates showUnpaired warns for distance", {
  fix <- .get_model_fixture()

  out <- suppressWarnings(
    ModelTCRDyads(
      fix$seuratObj_TCR,
      chains            = "TRB",
      groupColumn       = "outcome",
      distanceThreshold = 200,
      communityMethod   = "threshold",
      subjectIdCol      = "SubjectId"
    )
  )

  expect_warning(
    PlotDyadEstimates(out, which = "distance", showUnpaired = TRUE),
    "showUnpaired is ignored"
  )
})
