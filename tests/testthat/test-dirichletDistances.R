library(testthat)
library(tcrClustR)

# helper to mock a seurat object with a distance matrix for tests
.make_dirichlet_seurat <- function(n_clones_per_group = 5, seed = 1) {
  set.seed(seed)

  clone_ids <- paste0("clone", seq_len(n_clones_per_group * 2))
  n <- length(clone_ids)

  # symmetric distance matrix with realistic-ish values (0-200 range)
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(clone_ids, clone_ids))
  upper_vals <- abs(rnorm(n * (n - 1) / 2, mean = 80, sd = 30))
  mat[upper.tri(mat)] <- upper_vals
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]

  dist_assay <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(mat))

  # one cell per clone, split evenly across two groups
  cell_ids <- paste0("cell", seq_len(n))
  meta <- data.frame(
    TRA_CloneIdx = clone_ids,
    group        = rep(c("GroupA", "GroupB"), each = n_clones_per_group),
    row.names    = cell_ids,
    stringsAsFactors = FALSE
  )

  dummy_counts <- matrix(1L, nrow = 2, ncol = n,
                         dimnames = list(c("gene1", "gene2"), cell_ids))
  seuratObj <- SeuratObject::CreateSeuratObject(
    counts    = Seurat::as.sparse(dummy_counts),
    meta.data = meta
  )
  seuratObj@misc$TCR_Distances <- list(TRA_fl = dist_assay)
  seuratObj
}


# ===========================================================================
# 1. .assay_to_chain_id()
# ===========================================================================
test_that(".assay_to_chain_id strips suffix correctly", {
  expect_equal(tcrClustR:::.assay_to_chain_id("TRA_fl"),       "TRA")
  expect_equal(tcrClustR:::.assay_to_chain_id("TRB_cdr3"),     "TRB")
  expect_equal(tcrClustR:::.assay_to_chain_id("TRA_TRB_fl"),   "TRA_TRB")
  expect_equal(tcrClustR:::.assay_to_chain_id("TRG_TRD_cdr3"), "TRG_TRD")
})


# ===========================================================================
# 2. .SampleDistanceVector()
# ===========================================================================
test_that(".SampleDistanceVector respects maxSamples", {
  out <- tcrClustR:::.SampleDistanceVector(seq_len(1000), maxSamples = 80, seed = 1)
  expect_lte(length(out), 80)
})

test_that(".SampleDistanceVector returns all values when input is small", {
  x <- 1:5
  out <- tcrClustR:::.SampleDistanceVector(x, maxSamples = 100, seed = 1)
  expect_equal(length(out), 5)
})

test_that(".SampleDistanceVector is reproducible with same seed", {
  x <- rnorm(500)
  out1 <- tcrClustR:::.SampleDistanceVector(x, maxSamples = 50, seed = 99)
  out2 <- tcrClustR:::.SampleDistanceVector(x, maxSamples = 50, seed = 99)
  expect_identical(out1, out2)
})


# ===========================================================================
# 3. ExtractGroupDistanceVectors() — validation
# ===========================================================================
test_that("ExtractGroupDistanceVectors errors without TCR_Distances", {
  seuratObj <- SeuratObject::CreateSeuratObject(
    counts = Seurat::as.sparse(matrix(1L, 2, 3,
      dimnames = list(c("g1", "g2"), paste0("c", 1:3))))
  )
  expect_error(ExtractGroupDistanceVectors(seuratObj, "TRA_fl", "group"),
               "TCR_Distances")
})

test_that("ExtractGroupDistanceVectors errors on missing assayName", {
  seuratObj <- .make_dirichlet_seurat()
  expect_error(ExtractGroupDistanceVectors(seuratObj, "TRB_fl", "group"), "TRB_fl")
})

test_that("ExtractGroupDistanceVectors errors on missing splitField", {
  seuratObj <- .make_dirichlet_seurat()
  expect_error(ExtractGroupDistanceVectors(seuratObj, "TRA_fl", "bad"), "bad")
})


# ===========================================================================
# 4. ExtractGroupDistanceVectors() — correct extraction
# ===========================================================================
test_that("ExtractGroupDistanceVectors returns one vector per group", {
  seuratObj <- .make_dirichlet_seurat(n_clones_per_group = 4)
  result <- ExtractGroupDistanceVectors(seuratObj, "TRA_fl", "group", verbose = FALSE)
  expect_type(result, "list")
  expect_setequal(names(result), c("GroupA", "GroupB"))
})

test_that("ExtractGroupDistanceVectors returns correct upper-triangle length", {
  n <- 4
  seuratObj <- .make_dirichlet_seurat(n_clones_per_group = n)
  result <- ExtractGroupDistanceVectors(seuratObj, "TRA_fl", "group", verbose = FALSE)
  expect_equal(length(result$GroupA), n * (n - 1) / 2)
})

test_that("ExtractGroupDistanceVectors drops groups below minClonesPerGroup", {
  seuratObj <- .make_dirichlet_seurat(n_clones_per_group = 1)
  result <- ExtractGroupDistanceVectors(seuratObj, "TRA_fl", "group",
                                        minClonesPerGroup = 2, verbose = FALSE)
  expect_equal(length(result), 0)
})


# ===========================================================================
# 5. .ExtractClusterParams()
# ===========================================================================
test_that(".ExtractClusterParams returns expected columns", {
  seuratObj <- .make_dirichlet_seurat(n_clones_per_group = 5)
  dists <- ExtractGroupDistanceVectors(seuratObj, "TRA_fl", "group", verbose = FALSE)
  m <- dirichletprocess::DirichletProcessGaussian(dists$GroupA)
  m <- dirichletprocess::Fit(m, 50)
  df <- tcrClustR:::.ExtractClusterParams(m, "GroupA")

  expect_true(is.data.frame(df))
  expect_named(df, c("Cluster", "Mu", "Sigma", "MixingProportion",
                      "PointsPerCluster", "Group"))
  expect_equal(nrow(df), m$numberClusters)
  expect_true(all(df$Mu > 0))
  expect_true(all(df$Sigma > 0))
  expect_equal(sum(df$MixingProportion), 1.0, tolerance = 1e-6)
})


# ===========================================================================
# 6. DirichletClusterAnalysis() — main entry point
# ===========================================================================
test_that("DirichletClusterAnalysis returns tcrDirichletResult with cluster_summary", {
  seuratObj <- .make_dirichlet_seurat(n_clones_per_group = 5)
  result <- DirichletClusterAnalysis(
    seuratObj   = seuratObj,
    assayName   = "TRA_fl",
    splitField  = "group",
    maxSamples  = 50,
    nIterations = 50,
    verbose     = FALSE
  )

  expect_s3_class(result, "tcrDirichletResult")
  expect_equal(result$assayName, "TRA_fl")
  expect_equal(result$splitField, "group")

  # cluster_summary is the primary deliverable
  cs <- result$cluster_summary
  expect_true(is.data.frame(cs))
  expect_named(cs, c("Cluster", "Mu", "Sigma", "MixingProportion",
                      "PointsPerCluster", "Group"))
  expect_setequal(unique(cs$Group), c("GroupA", "GroupB"))
})

test_that("DirichletClusterAnalysis accepts custom nBins and samplesPerBin", {
  seuratObj <- .make_dirichlet_seurat(n_clones_per_group = 5)
  result <- DirichletClusterAnalysis(
    seuratObj     = seuratObj,
    assayName     = "TRA_fl",
    splitField    = "group",
    maxSamples    = 50,
    nIterations   = 50,
    nBins         = 5,
    samplesPerBin = 10,
    verbose       = FALSE
  )
  expect_s3_class(result, "tcrDirichletResult")
  cs <- result$cluster_summary
  expect_true(nrow(cs) >= 2)  # at least one cluster per group
})

test_that("DirichletClusterAnalysis errors when no groups pass minClonesPerGroup", {
  seuratObj <- .make_dirichlet_seurat(n_clones_per_group = 1)
  expect_error(
    DirichletClusterAnalysis(seuratObj, "TRA_fl", "group",
                             minClonesPerGroup = 2, nIterations = 50, verbose = FALSE),
    "No groups"
  )
})


# ===========================================================================
# 7. PlotClusterMeans() and PlotMixingProportions()
# ===========================================================================
test_that("PlotClusterMeans returns a ggplot", {
  seuratObj <- .make_dirichlet_seurat(n_clones_per_group = 5)
  dp <- DirichletClusterAnalysis(seuratObj, "TRA_fl", "group",
                                  maxSamples = 50, nIterations = 50, verbose = FALSE)
  p <- PlotClusterMeans(dp)
  expect_s3_class(p, "ggplot")
  expect_true(grepl("TRA_fl", p$labels$title))
})

test_that("PlotMixingProportions returns a ggplot", {
  seuratObj <- .make_dirichlet_seurat(n_clones_per_group = 5)
  dp <- DirichletClusterAnalysis(seuratObj, "TRA_fl", "group",
                                  maxSamples = 50, nIterations = 50, verbose = FALSE)
  p <- PlotMixingProportions(dp)
  expect_s3_class(p, "ggplot")
  expect_true(grepl("TRA_fl", p$labels$title))
})

test_that("Plot functions error on wrong input type", {
  expect_error(PlotClusterMeans(list()),        "tcrDirichletResult")
  expect_error(PlotMixingProportions(list()),   "tcrDirichletResult")
})
