library(testthat)

#setup function to create a Seurat object with TCR data
setup_test_data <- function() {
  temp_dir <- tempdir()
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1)))

  rdsOutputPath <- file.path(temp_dir, "tcrdist3DistanceMatrices")

  seuratObj_TCR <- RunTcrdist3(inputData = seuratObj,
                               chains = c("TRA", "TRB"),
                               minimumClonesPerSubject = 2,
                               rdsOutputPath = rdsOutputPath)

  return(list(seuratObj = seuratObj, seuratObj_TCR = seuratObj_TCR))
}

#create test distance matrix
create_test_matrix <- function(n_samples = 20) {
  set.seed(1234)
  distance_matrix <- matrix(runif(n_samples^2, 0, 100), nrow = n_samples, ncol = n_samples)
  distance_matrix <- (distance_matrix + t(distance_matrix)) / 2
  diag(distance_matrix) <- 0
  return(distance_matrix)
}

test_that("ClusterTcrs core functionality works", {
  test_data <- setup_test_data()

  #test with standard PCA
  testthat::expect_no_error({
    result_pca <- ClusterTcrs(
      seuratObj_TCR = test_data$seuratObj_TCR,
      pcaComponents = 5,
      usePCA = TRUE,
      computeMultiChain = FALSE
    )
  })

  #test with kernel PCA
  testthat::expect_no_error({
    result_kpca <- ClusterTcrs(
      seuratObj_TCR = test_data$seuratObj_TCR,
      pcaComponents = 5,
      usePCA = FALSE,
      computeMultiChain = FALSE
    )
  })

  #validate results structure
  testthat::expect_true(is.list(result_pca))
  testthat::expect_true("singleChainSeuratObject" %in% names(result_pca))
  testthat::expect_true(is.list(result_kpca))
  testthat::expect_true("singleChainSeuratObject" %in% names(result_kpca))
})

test_that(".PcaAndClustering works with both PCA types", {
  distance_matrix <- create_test_matrix()

  #test standard PCA
  testthat::expect_no_error({
    result_pca <- .PcaAndClustering(
      distanceMatrix = distance_matrix,
      pcaComponents = 5,
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1
    )
  })

  #test kernel PCA
  testthat::expect_no_error({
    result_kpca <- .PcaAndClustering(
      distanceMatrix = distance_matrix,
      pcaComponents = 5,
      usePCA = FALSE,
      kpcaKernel = "rbfdot",
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1
    )
  })

  #validate results
  testthat::expect_true("graph" %in% names(result_pca))
  testthat::expect_true("pca_result" %in% names(result_pca))
  testthat::expect_true(inherits(result_pca$pca_result, "pca_result"))
  testthat::expect_true(inherits(result_kpca$pca_result, "kpca"))
  testthat::expect_true(igraph::is_igraph(result_pca$graph))
  testthat::expect_true(igraph::is_igraph(result_kpca$graph))
})

test_that(".TranslateGroupByVariablesToTcrdist3 works correctly", {
  #test basic translations
  test_variables <- c("TRA_V", "TRA_J", "TRA", "TRB_V", "TRB_J", "TRB")
  expected_output <- c("v_a_gene", "j_a_gene", "cdr3_a_aa", "v_b_gene", "j_b_gene", "cdr3_b_aa")

  result <- .TranslateGroupByVariablesToTcrdist3(test_variables)
  testthat::expect_equal(result, expected_output)

  #test gamma/delta chains
  test_variables_gd <- c("TRG_V", "TRG_J", "TRG", "TRD_V", "TRD_J", "TRD")
  expected_output_gd <- c("v_g_gene", "j_g_gene", "cdr3_g_aa", "v_d_gene", "j_d_gene", "cdr3_d_aa")

  result_gd <- .TranslateGroupByVariablesToTcrdist3(test_variables_gd)
  testthat::expect_equal(result_gd, expected_output_gd)

  #test empty input
  testthat::expect_equal(.TranslateGroupByVariablesToTcrdist3(character(0)), character(0))
})

test_that(".CreateTcrKeyLookup works correctly", {
  test_data <- setup_test_data()
  assays <- Seurat::Assays(test_data$seuratObj_TCR)

  for (assay in assays) {
    testthat::expect_no_error({
      lookup <- .CreateTcrKeyLookup(test_data$seuratObj_TCR, assay)
    })

    #validate lookup structure
    testthat::expect_true(is.data.frame(lookup))
    testthat::expect_true("key" %in% colnames(lookup))
    testthat::expect_true("matrix_rowname" %in% colnames(lookup))
    testthat::expect_gt(nrow(lookup), 0)
    testthat::expect_true(!all(is.na(lookup$key)))
  }
})

test_that("Parameter validation works correctly", {
  distance_matrix <- create_test_matrix(n_samples = 10)

  #test with too many components (should handle gracefully)
  testthat::expect_no_error({
    result <- .PcaAndClustering(
      distanceMatrix = distance_matrix,
      pcaComponents = 50,  # More than samples
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1
    )
  })

  #test edge cases for neighbor proportion
  testthat::expect_no_error({
    result <- .PcaAndClustering(
      distanceMatrix = distance_matrix,
      pcaComponents = 5,
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.9,  # High proportion
      jaccardIndexThreshold = 0.1
    )
  })
})

test_that("Multi-chain functionality works", {
  test_data <- setup_test_data()

  #test with multi-chain enabled
  testthat::expect_no_error({
    result <- .DistanceMatrixToClusteredGraphs(
      seuratObj_TCR = test_data$seuratObj_TCR,
      pcaComponents = 5,
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1,
      resolutions = c(0.1),
      seed = 1234,
      computeMultiChain = TRUE
    )
  })

  #validate structure
  testthat::expect_true(is.list(result))
  testthat::expect_true("singleChainSeuratObject" %in% names(result))
})

test_that("Spike-in functionality works", {
  #create spike-in dataframe
  spikeInDataframe <- data.frame(
    CloneNames = rep(1:3),
    TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
    TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
    TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
    TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
    TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
    TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY")
  )

  temp_dir <- tempdir()
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1)))

    rdsOutputPath <- file.path(temp_dir, "tcrdist3DistanceMatrices_spikeins")

  testthat::expect_no_error({
    seuratObj_TCR <- RunTcrdist3(inputData = seuratObj,
                                 chains = c("TRA", "TRB"),
                                 minimumClonesPerSubject = 2,
                                 rdsOutputPath = rdsOutputPath,
                                 spikeInDataframe = spikeInDataframe)
  })

  # TODO: verify spike-ins are present
  tcrdist3Input <- readr::read_csv(postFormattingMetadataCsvPath, show_col_types = FALSE)
  testthat::expect_true(sum(grepl("spikeIn", tcrdist3Input$subject)) == 3)
})
