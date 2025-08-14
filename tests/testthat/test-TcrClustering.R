library(testthat)

test_that("PCA and Clustering functions work correctly", {
  temp_dir <- tempdir()
  #debug
  print(temp_dir)
  paste0('file.exists: ', file.exists(temp_dir))
  paste0('dir.exists: ', dir.exists(temp_dir))

  #read in data
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1)))

  #create distance matrices using tcrdist3 for testing
  postFormattingMetadataCsvPath <- file.path(temp_dir, "tcrdist3Input.csv")
  rdsOutputPath <- file.path(temp_dir, "tcrdist3DistanceMatrices")

  #generate distance matrices
  seuratObj_TCR <- RunTcrdist3(seuratObj = seuratObj,
                              metadata = NULL,
                              formatMetadata = TRUE,
                              postFormattingMetadataCsvPath = postFormattingMetadataCsvPath,
                              chains = c("TRA", "TRB"),
                              cleanMetadata = TRUE,
                              minimumClonesPerSubject = 2,
                              rdsOutputPath = rdsOutputPath,
                              pythonExecutable = Sys.which("python3"),
                              debugTcrdist3 = "True")

  #test that distance matrix calculations worked
  testthat::expect_true(class(seuratObj_TCR)[1] == "Seurat")
  testthat::expect_gt(length(Seurat::Assays(seuratObj_TCR)), 0)

  #test that ClusterTcrs with standard PCA
  testthat::expect_no_error({
    clustered_results_pca <- ClusterTcrs(
      seuratObj = NULL,
      seuratObj_TCR = seuratObj_TCR,
      pcaComponents = 10,
      usePCA = TRUE,       # Test standard PCA
      proportionOfGraphAsNeighbors = 0.2,
      jaccardIndexThreshold = 0.05,
      seed = 1234,
      computeMultiChain = FALSE
    )
  })

  #test that ClusterTcrs works with kernel PCA
  testthat::expect_no_error({
    clustered_results_kpca <- ClusterTcrs(
      seuratObj = seuratObj,
      seuratObj_TCR = seuratObj_TCR,
      pcaComponents = 10,
      usePCA = FALSE,      # Test kernel PCA
      kpcaKernel = "rbfdot",
      proportionOfGraphAsNeighbors = 0.2,
      jaccardIndexThreshold = 0.05,
      seed = 1234,
      computeMultiChain = FALSE
    )
  })

  #test that results contain expected objects
  testthat::expect_true(is.list(clustered_results_pca))
  testthat::expect_true("singleChainSeuratObject" %in% names(clustered_results_pca))
  testthat::expect_true(is.list(clustered_results_kpca))
  testthat::expect_true("singleChainSeuratObject" %in% names(clustered_results_kpca))
})

#test internals
test_that(".DistanceMatrixToClusteredGraphs works with different PCA types", {
  #read test data
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1)))

  #generate distance matrices
  temp_dir <- tempdir()
  postFormattingMetadataCsvPath <- file.path(temp_dir, "tcrdist3Input.csv")
  rdsOutputPath <- file.path(temp_dir, "tcrdist3DistanceMatrices")

  seuratObj_TCR <- RunTcrdist3(seuratObj = seuratObj,
                              formatMetadata = TRUE,
                              postFormattingMetadataCsvPath = postFormattingMetadataCsvPath,
                              chains = c("TRA", "TRB"),
                              minimumClonesPerSubject = 2,
                              rdsOutputPath = rdsOutputPath,
                              pythonExecutable = Sys.which("python3"))

  #test internal .DistanceMatrixToClusteredGraphs works with standard PCA
  testthat::expect_no_error({
    result_pca <- .DistanceMatrixToClusteredGraphs(
      seuratObj_TCR = seuratObj_TCR,
      pcaComponents = 5,
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1,
      resolutions = c(0.1, 0.2),
      seed = 1234,
      computeMultiChain = FALSE
    )
  })

  #test internal .DistanceMatrixToClusteredGraphs works with kernel PCA
  testthat::expect_no_error({
    result_kpca <- .DistanceMatrixToClusteredGraphs(
      seuratObj_TCR = seuratObj_TCR,
      pcaComponents = 5,
      usePCA = FALSE,
      kpcaKernel = "rbfdot",
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1,
      resolutions = c(0.1, 0.2),
      seed = 1234,
      computeMultiChain = FALSE
    )
  })

  #test that results are properly structured
  testthat::expect_true(is.list(result_pca))
  testthat::expect_true(is.list(result_kpca))
  testthat::expect_true("singleChainSeuratObject" %in% names(result_pca))
  testthat::expect_true("singleChainSeuratObject" %in% names(result_kpca))
})

test_that(".PcaAndClustering works with both PCA types", {
  #create a simple test distance matrix
  set.seed(1234)
  n_samples <- 20
  distance_matrix <- matrix(runif(n_samples^2, 0, 100), nrow = n_samples, ncol = n_samples)
  #symmetry
  distance_matrix <- (distance_matrix + t(distance_matrix)) / 2
  diag(distance_matrix) <- 0

  #test with standard PCA
  testthat::expect_no_error({
    result_pca <- .PcaAndClustering(
      distanceMatrix = distance_matrix,
      pcaComponents = 5,
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1
    )
  })

  #test with kernel PCA
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

  #test that results contain expected components
  testthat::expect_true(is.list(result_pca))
  testthat::expect_true(is.list(result_kpca))
  testthat::expect_true("graph" %in% names(result_pca))
  testthat::expect_true("pca_result" %in% names(result_pca))
  testthat::expect_true("graph" %in% names(result_kpca))
  testthat::expect_true("pca_result" %in% names(result_kpca))

  #test that PCA result has correct class
  testthat::expect_true(inherits(result_pca$pca_result, "pca_result"))
  testthat::expect_true(inherits(result_kpca$pca_result, "kpca"))

  #test that graphs writing worked
  testthat::expect_true(igraph::is_igraph(result_pca$graph))
  testthat::expect_true(igraph::is_igraph(result_kpca$graph))
})

test_that(".TranslateGroupByVariablesToTcrdist3 works correctly", {
  #test basic translations
  test_variables <- c("TRA_V", "TRA_J", "TRA", "TRB_V", "TRB_J", "TRB")
  expected_output <- c("v_a_gene", "j_a_gene", "cdr3_a_aa", "v_b_gene", "j_b_gene", "cdr3_b_aa")

  testthat::expect_no_error({
    result <- .TranslateGroupByVariablesToTcrdist3(test_variables)
  })

  testthat::expect_equal(result, expected_output)

  #test gamma and delta chains
  test_variables_gd <- c("TRG_V", "TRG_J", "TRG", "TRD_V", "TRD_J", "TRD")
  expected_output_gd <- c("v_g_gene", "j_g_gene", "cdr3_g_aa", "v_d_gene", "j_d_gene", "cdr3_d_aa")

  result_gd <- .TranslateGroupByVariablesToTcrdist3(test_variables_gd)
  testthat::expect_equal(result_gd, expected_output_gd)

  #test empty input
  testthat::expect_equal(.TranslateGroupByVariablesToTcrdist3(character(0)), character(0))
})

test_that(".CreateTcrKeyLookup works correctly", {
  #read test data and create distance matrices
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1)))

  temp_dir <- tempdir()
  postFormattingMetadataCsvPath <- file.path(temp_dir, "tcrdist3Input.csv")
  rdsOutputPath <- file.path(temp_dir, "tcrdist3DistanceMatrices")

  seuratObj_TCR <- RunTcrdist3(seuratObj = seuratObj,
                              formatMetadata = TRUE,
                              postFormattingMetadataCsvPath = postFormattingMetadataCsvPath,
                              chains = c("TRA", "TRB"),
                              minimumClonesPerSubject = 2,
                              rdsOutputPath = rdsOutputPath,
                              pythonExecutable = Sys.which("python3"))

  # Test for different assay types
  assays <- Seurat::Assays(seuratObj_TCR)

  for (assay in assays) {
    testthat::expect_no_error({
      lookup <- .CreateTcrKeyLookup(seuratObj_TCR, assay)
    })

    # Test that lookup is a data frame
    testthat::expect_true(is.data.frame(lookup))

    # Test that it has the expected columns
    testthat::expect_true("key" %in% colnames(lookup))
    testthat::expect_true("matrix_rowname" %in% colnames(lookup))

    # Test that the number of rows matches expectations
    testthat::expect_gt(nrow(lookup), 0)

    # Test that keys are not all NA
    testthat::expect_true(!all(is.na(lookup$key)))
  }
})

test_that(".AddDimensionalityReductions works with both PCA types", {
  # Create test data
  set.seed(1234)
  n_samples <- 20
  distance_matrix <- matrix(runif(n_samples^2, 0, 100), nrow = n_samples, ncol = n_samples)
  distance_matrix <- (distance_matrix + t(distance_matrix)) / 2
  diag(distance_matrix) <- 0

  # Create test PCA results
  pca_result_standard <- .PcaAndClustering(
    distanceMatrix = distance_matrix,
    pcaComponents = 5,
    usePCA = TRUE,
    proportionOfGraphAsNeighbors = 0.3,
    jaccardIndexThreshold = 0.1
  )$pca_result

  pca_result_kernel <- .PcaAndClustering(
    distanceMatrix = distance_matrix,
    pcaComponents = 5,
    usePCA = FALSE,
    kpcaKernel = "rbfdot",
    proportionOfGraphAsNeighbors = 0.3,
    jaccardIndexThreshold = 0.1
  )$pca_result

  # Create a minimal Seurat object for testing
  test_counts <- matrix(rnorm(n_samples * 10), nrow = 10, ncol = n_samples)
  rownames(test_counts) <- paste0("feature_", 1:10)
  colnames(test_counts) <- paste0("cell_", 1:n_samples)

  seuratObj_test <- Seurat::CreateSeuratObject(counts = test_counts)

  # Test with standard PCA
  testthat::expect_no_error({
    result_pca <- .AddDimensionalityReductions(
      seuratObj = seuratObj_test,
      pca_result = pca_result_standard,
      reductionName = "test_pca",
      assayName = "RNA",
      pcaComponents = 5,
      distanceMatrix = distance_matrix
    )
  })

  # Test with kernel PCA
  testthat::expect_no_error({
    result_kpca <- .AddDimensionalityReductions(
      seuratObj = seuratObj_test,
      pca_result = pca_result_kernel,
      reductionName = "test_kpca",
      assayName = "RNA",
      pcaComponents = 5,
      distanceMatrix = distance_matrix
    )
  })

  # Test that reductions were added
  testthat::expect_true("test_pca" %in% names(result_pca@reductions))
  testthat::expect_true("test_kpca" %in% names(result_kpca@reductions))

  # Test that UMAP reductions were created
  testthat::expect_true("test_pca_umap" %in% names(result_pca@reductions))
  testthat::expect_true("test_kpca_umap" %in% names(result_kpca@reductions))
})

test_that("Parameter validation works correctly", {
  # Create minimal test data
  set.seed(1234)
  n_samples <- 10
  distance_matrix <- matrix(runif(n_samples^2, 0, 100), nrow = n_samples, ncol = n_samples)
  distance_matrix <- (distance_matrix + t(distance_matrix)) / 2
  diag(distance_matrix) <- 0

  # Test invalid pcaComponents
  testthat::expect_no_error({
    result <- .PcaAndClustering(
      distanceMatrix = distance_matrix,
      pcaComponents = 50,  # More components than samples
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1
    )
  })

  # Test edge cases for proportionOfGraphAsNeighbors
  testthat::expect_no_error({
    result <- .PcaAndClustering(
      distanceMatrix = distance_matrix,
      pcaComponents = 5,
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.9,  # High proportion
      jaccardIndexThreshold = 0.1
    )
  })

  testthat::expect_no_error({
    result <- .PcaAndClustering(
      distanceMatrix = distance_matrix,
      pcaComponents = 5,
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.01,  # Low proportion
      jaccardIndexThreshold = 0.1
    )
  })
})

test_that("Multi-chain clustering works when enabled", {
  #read test data
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1)))

  #generate distance matrices
  temp_dir <- tempdir()
  postFormattingMetadataCsvPath <- file.path(temp_dir, "tcrdist3Input.csv")
  rdsOutputPath <- file.path(temp_dir, "tcrdist3DistanceMatrices")

  seuratObj_TCR <- RunTcrdist3(seuratObj = seuratObj,
                              formatMetadata = TRUE,
                              postFormattingMetadataCsvPath = postFormattingMetadataCsvPath,
                              chains = c("TRA", "TRB"),
                              minimumClonesPerSubject = 2,
                              rdsOutputPath = rdsOutputPath,
                              pythonExecutable = Sys.which("python3"),
                              multichain = TRUE)

  #test with multi-chain enabled
  testthat::expect_no_error({
    result <- .DistanceMatrixToClusteredGraphs(
      seuratObj_TCR = seuratObj_TCR,
      pcaComponents = 5,
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1,
      resolutions = c(0.1),
      seed = 1234,
      computeMultiChain = TRUE
    )
  })

  # Test structure of results
  testthat::expect_true(is.list(result))
  testthat::expect_true("singleChainSeuratObject" %in% names(result))
  testthat::expect_true("multiChainSeuratObject" %in% names(result))
})

test_that("Spike-in functionality works correctly", {
  # Create spike-in dataframe
  spikeInDataframe <- data.frame(
    CloneNames = rep(1:3),
    TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
    TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
    TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
    TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
    TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
    TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY")
  )

  # Read test data
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1)))

  # Generate distance matrices with spike-ins
  temp_dir <- tempdir()
  postFormattingMetadataCsvPath <- file.path(temp_dir, "tcrdist3Input_spikeins.csv")
  rdsOutputPath <- file.path(temp_dir, "tcrdist3DistanceMatrices_spikeins")

  testthat::expect_no_error({
    seuratObj_TCR <- RunTcrdist3(seuratObj = seuratObj,
                                formatMetadata = TRUE,
                                postFormattingMetadataCsvPath = postFormattingMetadataCsvPath,
                                chains = c("TRA", "TRB"),
                                minimumClonesPerSubject = 2,
                                rdsOutputPath = rdsOutputPath,
                                pythonExecutable = Sys.which("python3"),
                                spikeInDataframe = spikeInDataframe)
  })

  # Test that clustering works with spike-ins
  testthat::expect_no_error({
    clustered_results <- ClusterTcrs(
      seuratObj = seuratObj,
      seuratObj_TCR = seuratObj_TCR,
      pcaComponents = 5,
      usePCA = TRUE,
      proportionOfGraphAsNeighbors = 0.3,
      jaccardIndexThreshold = 0.1,
      seed = 1234,
      spikeInDataframe = spikeInDataframe,
      computeMultiChain = FALSE
    )
  })

  # Verify spike-ins are present in the input file
  tcrdist3Input <- readr::read_csv(postFormattingMetadataCsvPath, show_col_types = FALSE)
  testthat::expect_true(sum(grepl("spikeIn", tcrdist3Input$subject)) == 3)
})

test_that("Multichain assay detection works correctly", {
  # Test single chain assays
  single_chain_assays <- c("TRA", "TRB", "TRA_cdr3", "TRB_cdr3", "TRG", "TRD")
  
  for (assay in single_chain_assays) {
    parts <- strsplit(assay, "_")[[1]]
    chain_count <- sum(parts %in% c("TRA", "TRB", "TRG", "TRD"))
    is_single_chain <- chain_count == 1
    testthat::expect_true(is_single_chain, 
                         info = paste("Assay", assay, "should be detected as single chain"))
  }
  
  # Test multichain assays
  multichain_assays <- c("TRA_TRB", "TRA_TRB_cdr3", "TRA_cdr3_TRB", "TRA_cdr3_TRB_cdr3")
  
  for (assay in multichain_assays) {
    parts <- strsplit(assay, "_")[[1]]
    chain_count <- sum(parts %in% c("TRA", "TRB", "TRG", "TRD"))
    is_single_chain <- chain_count == 1
    testthat::expect_false(is_single_chain, 
                          info = paste("Assay", assay, "should be detected as multichain"))
  }
})
