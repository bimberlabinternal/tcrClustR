library(testthat)
library(tcrClustR)

test_that("JoinClusteringResults handles single-chain parquet files", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("Seurat")
  
  # Create mock Seurat object with TCR metadata
  mock_metadata <- data.frame(
    row.names = paste0("cell_", 1:10),
    TRA_V = c("TRAV1-2", "TRAV8-1", "TRAV1-2", "TRAV8-1", "TRAV1-2",
              "TRAV8-1", "TRAV1-2", "TRAV8-1", "TRAV1-2", "TRAV8-1"),
    TRA_J = c("TRAJ33", "TRAJ21", "TRAJ33", "TRAJ21", "TRAJ33",
              "TRAJ21", "TRAJ33", "TRAJ21", "TRAJ33", "TRAJ21"),
    TRA = c("CAVRD", "CAVSL", "CAVRD", "CAVSL", "CAVRD",
            "CAVSL", "CAVRD", "CAVSL", "CAVRD", "CAVSL"),
    TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4", "TRBV6-4", "TRBV6-4",
              "TRBV6-4", "TRBV6-4", "TRBV6-4", "TRBV6-4", "TRBV6-4"),
    TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ1-1", "TRBJ2-1", "TRBJ1-1",
              "TRBJ2-1", "TRBJ1-1", "TRBJ2-1", "TRBJ1-1", "TRBJ2-1"),
    TRB = c("CASS1", "CASS2", "CASS1", "CASS2", "CASS1",
            "CASS2", "CASS1", "CASS2", "CASS1", "CASS2"),
    stringsAsFactors = FALSE
  )
  
  mock_seurat <- SeuratObject::CreateSeuratObject(counts = matrix(1:30, nrow = 3, ncol = 10,
                                                                   dimnames = list(paste0("gene", 1:3),
                                                                                 rownames(mock_metadata))),
                                                   meta.data = mock_metadata)
  
  # Create mock parquet file (single-chain TRA)
  temp_dir <- tempdir()
  parquet_file_tra <- file.path(temp_dir, "test_TRA_single.parquet")
  
  cluster_data_tra <- data.frame(
    chain = rep("TRA", 2),
    v_gene = c("TRAV1-2", "TRAV8-1"),
    j_gene = c("TRAJ33", "TRAJ21"),
    CDR3 = c("CAVRD", "CAVSL"),
    Cluster = c("1", "2"),
    Cluster_Size_Threshold = c(2, 2),
    Clustering_Method = c("DIANA", "DIANA"),
    stringsAsFactors = FALSE
  )
  
  arrow::write_parquet(cluster_data_tra, parquet_file_tra)
  
  # Test joining
  result_seurat <- JoinClusteringResults(mock_seurat,
                                          parquetFiles = parquet_file_tra,
                                          verbose = FALSE,
                                          stripAlleles = FALSE)
  
  # Check that new column was added
  expect_true("TcrFamily_TRA" %in% colnames(result_seurat@meta.data))
  
  # Check that result is a factor (JoinClusteringResults converts to ordered factors)
  expect_true(is.factor(result_seurat@meta.data$TcrFamily_TRA))
  
  # Check cluster assignments (cells with CAVRD should be cluster 1, CAVSL should be cluster 2)
  # Use as.character() to compare factor values to strings
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA[1]), "1")
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA[2]), "2")
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA[3]), "1")
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA[4]), "2")
  
  # Clean up
  unlink(parquet_file_tra)
})


test_that("JoinClusteringResults handles paired-chain parquet files (AND join)", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("Seurat")
  
  # Create mock Seurat object
  mock_metadata <- data.frame(
    row.names = paste0("cell_", 1:6),
    TRA_V = c("TRAV1-2", "TRAV8-1", "TRAV12-1", NA, "TRAV1-2", "TRAV8-1"),
    TRA_J = c("TRAJ33", "TRAJ21", "TRAJ45", NA, "TRAJ33", "TRAJ21"),
    TRA = c("CAVRD", "CAVSL", "CAVXX", NA, "CAVRD", "CAVSL"),
    TRB_V = c("TRBV6-4", "TRBV6-4", NA, "TRBV6-4", "TRBV6-4", "TRBV6-4"),
    TRB_J = c("TRBJ1-1", "TRBJ2-1", NA, "TRBJ2-1", "TRBJ1-1", "TRBJ2-1"),
    TRB = c("CASS1", "CASS2", NA, "CASS2", "CASS1", "CASS2"),
    stringsAsFactors = FALSE
  )
  
  mock_seurat <- SeuratObject::CreateSeuratObject(counts = matrix(1:18, nrow = 3, ncol = 6,
                                                                   dimnames = list(paste0("gene", 1:3),
                                                                                 rownames(mock_metadata))),
                                                   meta.data = mock_metadata)
  
  # Create mock parquet file (paired TRA_TRB)
  temp_dir <- tempdir()
  parquet_file_paired <- file.path(temp_dir, "test_TRA_TRB_paired.parquet")
  
  cluster_data_paired <- data.frame(
    chain_1 = c("TRA", "TRA"),
    v_gene_1 = c("TRAV1-2", "TRAV8-1"),
    j_gene_1 = c("TRAJ33", "TRAJ21"),
    CDR3_1 = c("CAVRD", "CAVSL"),
    chain_2 = c("TRB", "TRB"),
    v_gene_2 = c("TRBV6-4", "TRBV6-4"),
    j_gene_2 = c("TRBJ1-1", "TRBJ2-1"),
    CDR3_2 = c("CASS1", "CASS2"),
    Cluster = c("1", "2"),
    Cluster_Size_Threshold = c(2, 2),
    Clustering_Method = c("DIANA", "DIANA"),
    stringsAsFactors = FALSE
  )
  
  arrow::write_parquet(cluster_data_paired, parquet_file_paired)
  
  # Test joining
  result_seurat <- JoinClusteringResults(mock_seurat,
                                          parquetFiles = parquet_file_paired,
                                          verbose = FALSE,
                                          stripAlleles = FALSE)
  
  # Check that new column was added
  expect_true("TcrFamily_TRA_TRB" %in% colnames(result_seurat@meta.data))
  
  # Check that result is a factor (JoinClusteringResults converts to ordered factors)
  expect_true(is.factor(result_seurat@meta.data$TcrFamily_TRA_TRB))
  
  # Check cluster assignments under AND logic:
  # Cell 1: TRA=CAVRD + TRB=CASS1 -> cluster 1 (complete match)
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA_TRB[1]), "1")

  # Cell 2: TRA=CAVSL + TRB=CASS2 -> cluster 2 (complete match)
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA_TRB[2]), "2")

  # Cell 3: Only TRA present (TRB missing) -> "No_TCR_Data" (no paired TCR data for TRB)
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA_TRB[3]), "No_TCR_Data")

  # Cell 4: Only TRB present (TRA missing) -> "No_TCR_Data" (no paired TCR data for TRA)
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA_TRB[4]), "No_TCR_Data")

  # Cell 5: TRA=CAVRD + TRB=CASS1 -> cluster 1 (complete match)
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA_TRB[5]), "1")

  # Cell 6: TRA=CAVSL + TRB=CASS2 -> cluster 2 (complete match)
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA_TRB[6]), "2")
  
  # Clean up
  unlink(parquet_file_paired)
})


test_that("JoinClusteringResults strips alleles when requested", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("Seurat")
  
  # Create mock Seurat object with allele suffixes
  mock_metadata <- data.frame(
    row.names = paste0("cell_", 1:4),
    TRA_V = c("TRAV1-2*01", "TRAV8-1*02", "TRAV1-2*01", "TRAV8-1*02"),
    TRA_J = c("TRAJ33*01", "TRAJ21*01", "TRAJ33*01", "TRAJ21*01"),
    TRA = c("CAVRD", "CAVSL", "CAVRD", "CAVSL"),
    stringsAsFactors = FALSE
  )
  
  mock_seurat <- SeuratObject::CreateSeuratObject(counts = matrix(1:12, nrow = 3, ncol = 4,
                                                                   dimnames = list(paste0("gene", 1:3),
                                                                                 rownames(mock_metadata))),
                                                   meta.data = mock_metadata)
  
  # Create mock parquet file (alleles already stripped by RunTcrClustering)
  temp_dir <- tempdir()
  parquet_file <- file.path(temp_dir, "test_alleles.parquet")
  
  cluster_data <- data.frame(
    chain = rep("TRA", 2),
    v_gene = c("TRAV1-2", "TRAV8-1"),  # No alleles
    j_gene = c("TRAJ33", "TRAJ21"),    # No alleles
    CDR3 = c("CAVRD", "CAVSL"),
    Cluster = c("1", "2"),
    Cluster_Size_Threshold = c(2, 2),
    Clustering_Method = c("DIANA", "DIANA"),
    stringsAsFactors = FALSE
  )
  
  arrow::write_parquet(cluster_data, parquet_file)
  
  # Test joining WITH allele stripping (should match)
  result_seurat <- JoinClusteringResults(mock_seurat,
                                          parquetFiles = parquet_file,
                                          verbose = FALSE,
                                          stripAlleles = TRUE)
  
  # Check that result is a factor
  expect_true(is.factor(result_seurat@meta.data$TcrFamily_TRA))
  
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA[1]), "1")
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA[2]), "2")
  
  # Test joining WITHOUT allele stripping (should NOT match because metadata has *01, parquet doesn't)
  mock_seurat2 <- SeuratObject::CreateSeuratObject(counts = matrix(1:12, nrow = 3, ncol = 4,
                                                                    dimnames = list(paste0("gene", 1:3),
                                                                                  rownames(mock_metadata))),
                                                    meta.data = mock_metadata)
  
  result_seurat2 <- JoinClusteringResults(mock_seurat2,
                                           parquetFiles = parquet_file,
                                           verbose = FALSE,
                                           stripAlleles = FALSE)
  
  # All should be "LowFrequency" because alleles don't match (cells have valid TCR data but no cluster match)
  # The function sets unmatched cells with valid TCR data to "0", which gets remapped to "LowFrequency"
  expect_true(all(as.character(result_seurat2@meta.data$TcrFamily_TRA) == "LowFrequency"))
  
  # Clean up
  unlink(parquet_file)
})


test_that("JoinClusteringResults handles missing columns gracefully", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("Seurat")
  
  # Create mock Seurat object WITHOUT TRA columns
  mock_metadata <- data.frame(
    row.names = paste0("cell_", 1:4),
    TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4", "TRBV6-4"),
    TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ1-1", "TRBJ2-1"),
    TRB = c("CASS1", "CASS2", "CASS1", "CASS2"),
    stringsAsFactors = FALSE
  )
  
  mock_seurat <- SeuratObject::CreateSeuratObject(counts = matrix(1:12, nrow = 3, ncol = 4,
                                                                   dimnames = list(paste0("gene", 1:3),
                                                                                 rownames(mock_metadata))),
                                                   meta.data = mock_metadata)
  
  # Create parquet file for TRA (which doesn't exist in metadata)
  temp_dir <- tempdir()
  parquet_file <- file.path(temp_dir, "test_missing_cols.parquet")
  
  cluster_data <- data.frame(
    chain = rep("TRA", 2),
    v_gene = c("TRAV1-2", "TRAV8-1"),
    j_gene = c("TRAJ33", "TRAJ21"),
    CDR3 = c("CAVRD", "CAVSL"),
    Cluster = c("1", "2"),
    Cluster_Size_Threshold = c(2, 2),
    Clustering_Method = c("DIANA", "DIANA"),
    stringsAsFactors = FALSE
  )
  
  arrow::write_parquet(cluster_data, parquet_file)
  
  # Should handle gracefully with warning and return NAs -> remapped to "No_TCR_Data"
  expect_warning({
    result_seurat <- JoinClusteringResults(mock_seurat,
                                            parquetFiles = parquet_file,
                                            verbose = FALSE)
  }, "Missing metadata columns")
  
  # Column should be added but all "No_TCR_Data" (since original NAs get remapped)
  expect_true("TcrFamily_TRA" %in% colnames(result_seurat@meta.data))
  expect_true(is.factor(result_seurat@meta.data$TcrFamily_TRA))
  expect_true(all(as.character(result_seurat@meta.data$TcrFamily_TRA) == "No_TCR_Data"))
  
  # Clean up
  unlink(parquet_file)
})


test_that("JoinClusteringResults auto-detects parquet files", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("Seurat")
  
  # Create mock Seurat object
  mock_metadata <- data.frame(
    row.names = paste0("cell_", 1:4),
    TRA_V = c("TRAV1-2", "TRAV8-1", "TRAV1-2", "TRAV8-1"),
    TRA_J = c("TRAJ33", "TRAJ21", "TRAJ33", "TRAJ21"),
    TRA = c("CAVRD", "CAVSL", "CAVRD", "CAVSL"),
    stringsAsFactors = FALSE
  )
  
  mock_seurat <- SeuratObject::CreateSeuratObject(counts = matrix(1:12, nrow = 3, ncol = 4,
                                                                   dimnames = list(paste0("gene", 1:3),
                                                                                 rownames(mock_metadata))),
                                                   meta.data = mock_metadata)
  
  # Create temporary directory with parquet files
  temp_dir <- file.path(tempdir(), "test_parquet_dir")
  dir.create(temp_dir, showWarnings = FALSE)
  
  parquet_file1 <- file.path(temp_dir, "clustering_TRA.parquet")
  cluster_data <- data.frame(
    chain = rep("TRA", 2),
    v_gene = c("TRAV1-2", "TRAV8-1"),
    j_gene = c("TRAJ33", "TRAJ21"),
    CDR3 = c("CAVRD", "CAVSL"),
    Cluster = c("1", "2"),
    Cluster_Size_Threshold = c(2, 2),
    Clustering_Method = c("DIANA", "DIANA"),
    stringsAsFactors = FALSE
  )
  arrow::write_parquet(cluster_data, parquet_file1)
  
  # Test auto-detection
  result_seurat <- JoinClusteringResults(mock_seurat,
                                          parquetDir = temp_dir,
                                          verbose = FALSE,
                                          stripAlleles = FALSE)
  
  expect_true("TcrFamily_TRA" %in% colnames(result_seurat@meta.data))
  expect_true(is.factor(result_seurat@meta.data$TcrFamily_TRA))
  expect_equal(as.character(result_seurat@meta.data$TcrFamily_TRA[1]), "1")
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
})


test_that("JoinClusteringResults respects overwriteExisting parameter", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("Seurat")
  
  # Create mock Seurat object with existing cluster column
  mock_metadata <- data.frame(
    row.names = paste0("cell_", 1:4),
    TRA_V = c("TRAV1-2", "TRAV8-1", "TRAV1-2", "TRAV8-1"),
    TRA_J = c("TRAJ33", "TRAJ21", "TRAJ33", "TRAJ21"),
    TRA = c("CAVRD", "CAVSL", "CAVRD", "CAVSL"),
    TcrFamily_TRA = c("old1", "old2", "old1", "old2"),  # Existing cluster column
    stringsAsFactors = FALSE
  )
  
  mock_seurat <- SeuratObject::CreateSeuratObject(counts = matrix(1:12, nrow = 3, ncol = 4,
                                                                   dimnames = list(paste0("gene", 1:3),
                                                                                 rownames(mock_metadata))),
                                                   meta.data = mock_metadata)
  
  # Create parquet file with new cluster assignments (use numeric cluster IDs like real output)
  temp_dir <- tempdir()
  parquet_file <- file.path(temp_dir, "test_overwrite.parquet")
  
  cluster_data <- data.frame(
    chain = rep("TRA", 2),
    v_gene = c("TRAV1-2", "TRAV8-1"),
    j_gene = c("TRAJ33", "TRAJ21"),
    CDR3 = c("CAVRD", "CAVSL"),
    Cluster = c("1", "2"),  # Use numeric-style cluster IDs
    Cluster_Size_Threshold = c(2, 2),
    Clustering_Method = c("DIANA", "DIANA"),
    stringsAsFactors = FALSE
  )
  
  arrow::write_parquet(cluster_data, parquet_file)
  
  # Test with overwriteExisting = FALSE (should keep old values)
  result_seurat1 <- JoinClusteringResults(mock_seurat,
                                           parquetFiles = parquet_file,
                                           verbose = FALSE,
                                           overwriteExisting = FALSE,
                                           stripAlleles = FALSE)
  
  # Column should be unchanged (still a character, not factor) because skip happened
  expect_equal(as.character(result_seurat1@meta.data$TcrFamily_TRA[1]), "old1")  # Unchanged
  
  # Test with overwriteExisting = TRUE (should replace with new values)
  result_seurat2 <- JoinClusteringResults(mock_seurat,
                                           parquetFiles = parquet_file,
                                           verbose = FALSE,
                                           overwriteExisting = TRUE,
                                           stripAlleles = FALSE)
  
  # New values should be applied (and converted to factor)
  expect_true(is.factor(result_seurat2@meta.data$TcrFamily_TRA))
  expect_equal(as.character(result_seurat2@meta.data$TcrFamily_TRA[1]), "1")  # Updated
  
  # Clean up
  unlink(parquet_file)
})
