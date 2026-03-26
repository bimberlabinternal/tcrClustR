library(testthat)
library(Seurat)

# Helper: minimal 6-cell Seurat object with TRB metadata for network tests.
# Clone sizes are deliberately set to 1 so minimumCloneSize = 1 retains all clones.
.make_network_test_seurat <- function() {
  set.seed(42L)
  metadata <- data.frame(
    TRB_V   = c("TRBV7-9",       "TRBV7-9",       "TRBV6-4",       "TRBV6-4",       "TRBV5-1",       "TRBV5-1"),
    TRB_J   = c("TRBJ2-1",       "TRBJ2-1",       "TRBJ1-1",       "TRBJ1-1",       "TRBJ2-3",       "TRBJ2-3"),
    TRB     = c("CASSLGQAYEQYF", "CASSLAQAYEQYF", "CASSYLDTEAFF",  "CASSYLDAEAFF",  "CASSLAGRDTQYF", "CASSLAGRDTQYF"),
    outcome = c("Controller",    "Controller",    "Progressor",    "Progressor",    "Controller",    "Controller"),
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

test_that("BuildTCRDistanceGraph returns a tbl_graph with correct node count", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  tg <- BuildTCRDistanceGraph(seuratObj_TCR, chains = "TRB", distanceThreshold = 200)

  expect_s3_class(tg, "tbl_graph")
  expect_true(inherits(tg, "igraph"))

  expected_clones <- nrow(GetDistanceMatrix(seuratObj_TCR, "TRB"))
  expect_equal(igraph::vcount(tg), expected_clones)
})

test_that("BuildTCRDistanceGraph joins Seurat metadata onto node table", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  tg      <- BuildTCRDistanceGraph(seuratObj_TCR, chains = "TRB", distanceThreshold = 200)
  node_df <- as.data.frame(tg)
  expect_true("outcome" %in% colnames(node_df))
})

test_that("BuildTCRDistanceGraph edgeType = 'continuous' produces weight and norm_weight columns", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  tg      <- BuildTCRDistanceGraph(seuratObj_TCR, chains = "TRB",
                                   distanceThreshold = 200, edgeType = "continuous")
  edge_df <- igraph::as_data_frame(tidygraph::as.igraph(tg), what = "edges")

  expect_true("weight"      %in% colnames(edge_df))
  expect_true("norm_weight" %in% colnames(edge_df))
  expect_true(all(edge_df$norm_weight >= 0 & edge_df$norm_weight <= 1))
  expect_equal(max(edge_df$norm_weight), 1)
})

test_that("BuildTCRDistanceGraph respects distanceThreshold", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  tg_open   <- BuildTCRDistanceGraph(seuratObj_TCR, chains = "TRB", distanceThreshold = 300)
  tg_strict <- BuildTCRDistanceGraph(seuratObj_TCR, chains = "TRB", distanceThreshold = 1)

  expect_gte(igraph::ecount(tidygraph::as.igraph(tg_open)),
             igraph::ecount(tidygraph::as.igraph(tg_strict)))
})

test_that("BuildTCRDistanceGraph validates arguments", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  expect_error(
    BuildTCRDistanceGraph(seuratObj_TCR, "TRB", edgeType = "step"),
    "edgeType must be"
  )
  expect_error(
    BuildTCRDistanceGraph(seuratObj_TCR, "TRB", distanceThreshold = -1),
    "single positive number"
  )
})

test_that("TCRDistanceNetwork returns a ggplot for threshold + binary edges", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  result <- TCRDistanceNetwork(
    seuratObj_TCR,
    chains            = "TRB",
    distanceThreshold = 200,
    communityMethod   = "threshold",
    edgeType          = "binary",
    colorBy           = "outcome"
  )
  expect_s3_class(result$plot, "ggplot")
  expect_true(is.data.frame(result$layout))
  expect_true(all(c("x", "y") %in% colnames(result$layout)))
  expect_true(is.data.frame(result$edges))
  expect_true(all(c("distance", "weight", "from_clone", "to_clone", "component") %in% colnames(result$edges)))
})

test_that("TCRDistanceNetwork returns a ggplot for continuous edges and colorByCommunity", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  result <- TCRDistanceNetwork(
    seuratObj_TCR,
    chains            = "TRB",
    distanceThreshold = 200,
    communityMethod   = "threshold",
    edgeType          = "continuous",
    colorByCommunity  = TRUE
  )
  expect_s3_class(result$plot, "ggplot")
  expect_true(is.data.frame(result$layout))
  expect_true(is.data.frame(result$edges))
  expect_true(all(c("distance", "weight", "norm_weight", "from_clone", "to_clone", "component") %in% colnames(result$edges)))
})

test_that("TCRDistanceNetwork works with communityMethod = 'DIANA'", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
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

  result <- TCRDistanceNetwork(
    seuratObj_TCR,
    chains            = "TRB",
    distanceThreshold = 200,
    communityMethod   = "DIANA",
    colorByCommunity  = TRUE
  )
  expect_s3_class(result$plot, "ggplot")
  expect_true(is.data.frame(result$layout))
  expect_true(is.data.frame(result$edges))
  expect_true(all(c("from_clone", "to_clone", "component") %in% colnames(result$edges)))
})

test_that("TCRDistanceNetwork works with showIsolated = TRUE", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  # very low threshold ensures some nodes will be isolated
  result <- TCRDistanceNetwork(
    seuratObj_TCR,
    chains            = "TRB",
    distanceThreshold = 1,
    communityMethod   = "threshold",
    edgeType          = "binary",
    showIsolated      = TRUE
  )
  expect_s3_class(result$plot, "ggplot")
  expect_true(is.data.frame(result$layout))
  expect_true(is.data.frame(result$edges))
  # all expected columns present even when there are zero edges
  expect_true(all(c("distance", "weight", "from_clone", "to_clone", "component") %in% colnames(result$edges)))
})

test_that("TCRDistanceNetwork communityDistanceThreshold uses a separate community threshold", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  # wide edge threshold, but restrict Louvain community graph to very tight distances
  result <- TCRDistanceNetwork(
    seuratObj_TCR,
    chains                     = "TRB",
    distanceThreshold          = 200,
    communityDistanceThreshold = 10,
    communityMethod            = "threshold",
    edgeType                   = "binary",
    colorByCommunity           = TRUE
  )
  expect_s3_class(result$plot, "ggplot")
  expect_true(is.data.frame(result$layout))
  # community column must be present on every node
  expect_true(".community" %in% colnames(result$layout))
})

test_that("TCRDistanceNetwork accepts a custom colorPalette passthrough", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  custom_pal <- c("Controller" = "#2166AC", "Progressor" = "#D6604D")
  result <- TCRDistanceNetwork(
    seuratObj_TCR,
    chains            = "TRB",
    distanceThreshold = 200,
    colorBy           = "outcome",
    colorPalette      = custom_pal
  )
  expect_s3_class(result$plot, "ggplot")
  expect_true(is.data.frame(result$layout))
  expect_true(is.data.frame(result$edges))
  expect_true(all(c("distance", "weight", "from_clone", "to_clone", "component") %in% colnames(result$edges)))
})

test_that("TCRDistanceNetwork validates arguments", {
  seuratObj_TCR <- CalculateTcrDistances(
    .make_network_test_seurat(),
    chains           = "TRB",
    minimumCloneSize = 1L,
    verbose          = FALSE
  )

  expect_error(
    TCRDistanceNetwork(seuratObj_TCR, "TRB", communityMethod = "kmeans"),
    "communityMethod must be"
  )
  expect_error(
    TCRDistanceNetwork(seuratObj_TCR, "TRB", edgeType = "weighted"),
    "edgeType must be"
  )
  expect_error(
    TCRDistanceNetwork(seuratObj_TCR, "TRB", communityMethod = "DIANA"),
    "Run RunTcrClustering"
  )
  expect_error(
    TCRDistanceNetwork(seuratObj_TCR, "TRB", colorBy = "nonexistent_col"),
    "colorBy"
  )
  expect_error(
    TCRDistanceNetwork(seuratObj_TCR, "TRB", labelBy = "nonexistent_col"),
    "labelBy"
  )
  expect_error(
    TCRDistanceNetwork(seuratObj_TCR, "TRB", nodeSizeBy = "nonexistent_col"),
    "nodeSizeBy"
  )
})
