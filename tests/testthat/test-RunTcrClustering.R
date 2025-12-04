library(testthat)
library(Seurat)
library(arrow)

test_that("RunTcrClustering filters invalid clones correctly", {
  skip_if_not(nzchar(Sys.getenv("RETICULATE_PYTHON")) || nzchar(Sys.which("python3")) || nzchar(Sys.which("python")), "Python not available")
  
  set.seed(123)
  
  # Create test data with valid and invalid clones
  metadata <- rbind(
    data.frame(TRB_V = "TRBV7-9", TRB_J = "TRBJ2-1", TRB = "CASSLGQAYEQYF", SubjectId = "S1"),
    data.frame(TRB_V = "TRBV7-9", TRB_J = "TRBJ2-1", TRB = "CASSLAQAYEQYF", SubjectId = "S2"),
    data.frame(TRB_V = "TRBV6-4", TRB_J = "TRBJ1-1", TRB = NA_character_, SubjectId = "S3"),  # NA CDR3
    data.frame(TRB_V = "TRBV5-1", TRB_J = "TRBJ2-3", TRB = "CASSA,CASSB", SubjectId = "S4"), # Concatenated
    data.frame(TRB_V = "TRBV1-1;TRBV2-1", TRB_J = "TRBJ1-1", TRB = "CASSC", SubjectId = "S5"), # Concatenated V
    data.frame(TRB_V = "TRBV6-4", TRB_J = "TRBJ1-1", TRB = "CASSYLDTEAFF", SubjectId = "S6"),
    data.frame(TRB_V = "TRBV6-4", TRB_J = "TRBJ1-1", TRB = "CASSYLDAEAFF", SubjectId = "S7")
  )
  rownames(metadata) <- paste0("clone_", seq_len(nrow(metadata)))
  
  dummy_counts <- matrix(rnorm(nrow(metadata) * 10), nrow = 10, ncol = nrow(metadata),
                         dimnames = list(paste0("gene_", 1:10), rownames(metadata)))
  seurat_obj <- CreateSeuratObject(counts = dummy_counts, meta.data = metadata)
  
  # Run tcrdist3
  tcr_obj <- RunTcrdist3(seurat_obj, chains = "TRB",
                         calculateChainPairs = FALSE,
                         minimumCloneSize = 1,
                         verbose = FALSE,
                         )
  
  # Test with filtering
  output_dir <- file.path(tempdir(), "test_filtered")
  unlink(output_dir, recursive = TRUE)
  dir.create(output_dir, recursive = TRUE)
  
  results <- RunTcrClustering(
    seuratObj = tcr_obj,
    assays = "TRB",
    clusteringMethod = "DIANA",
    dianaHeight = 20,
    clusterSizeThreshold = 2,
    outputDir = output_dir,
    outputPrefix = "with_filter",
    filterInvalidClones = TRUE,
    verbose = FALSE
  )
  
  df <- read_parquet(results$parquet_files[1])
  
  # Should have 4 valid clones (S1, S2, S6, S7)
  expect_equal(nrow(df), 4)
  expect_false(any(is.na(df$CDR3)))
  expect_false(any(grepl(",|;", df$v_gene)))
  expect_false(any(grepl(",|;", df$j_gene)))
  expect_false(any(grepl(",|;", df$CDR3)))
})

test_that("RunTcrClustering strips alleles correctly", {
  skip_if_not(nzchar(Sys.getenv("RETICULATE_PYTHON")) || nzchar(Sys.which("python3")) || nzchar(Sys.which("python")), "Python not available")
  
  set.seed(456)
  
  metadata <- rbind(
    data.frame(TRB_V = "TRBV7-9", TRB_J = "TRBJ2-1", TRB = "CASSA", SubjectId = "S1"),
    data.frame(TRB_V = "TRBV6-4", TRB_J = "TRBJ1-1", TRB = "CASSB", SubjectId = "S2"),
    data.frame(TRB_V = "TRBV5-1", TRB_J = "TRBJ2-3", TRB = "CASSC", SubjectId = "S3")
  )
  rownames(metadata) <- paste0("clone_", seq_len(nrow(metadata)))
  
  dummy_counts <- matrix(rnorm(nrow(metadata) * 10), nrow = 10, ncol = nrow(metadata),
                         dimnames = list(paste0("gene_", 1:10), rownames(metadata)))
  seurat_obj <- CreateSeuratObject(counts = dummy_counts, meta.data = metadata)
  
  tcr_obj <- RunTcrdist3(seurat_obj, chains = "TRB",
                         calculateChainPairs = FALSE,
                         minimumCloneSize = 1,
                         verbose = FALSE,
                         )
  
  # With allele stripping (default)
  output_dir_stripped <- file.path(tempdir(), "test_no_alleles")
  unlink(output_dir_stripped, recursive = TRUE)
  dir.create(output_dir_stripped, recursive = TRUE)
  
  results_stripped <- RunTcrClustering(
    seuratObj = tcr_obj,
    assays = "TRB",
    clusteringMethod = "DIANA",
    dianaHeight = 20,
    clusterSizeThreshold = 1,
    outputDir = output_dir_stripped,
    outputPrefix = "no_alleles",
    stripAlleles = TRUE,
    verbose = FALSE
  )
  
  df_stripped <- read_parquet(results_stripped$parquet_files[1])
  expect_false(any(grepl("\\*", df_stripped$v_gene)))
  expect_false(any(grepl("\\*", df_stripped$j_gene)))
  
  # Without allele stripping
  output_dir_with <- file.path(tempdir(), "test_with_alleles")
  unlink(output_dir_with, recursive = TRUE)
  dir.create(output_dir_with, recursive = TRUE)
  
  results_with <- RunTcrClustering(
    seuratObj = tcr_obj,
    assays = "TRB",
    clusteringMethod = "DIANA",
    dianaHeight = 20,
    clusterSizeThreshold = 1,
    outputDir = output_dir_with,
    outputPrefix = "with_alleles",
    stripAlleles = FALSE,
    verbose = FALSE
  )
  
  df_with <- read_parquet(results_with$parquet_files[1])
  expect_true(all(grepl("\\*", df_with$v_gene)))
  expect_true(all(grepl("\\*", df_with$j_gene)))
})

test_that(".filter_clones_for_assay handles single vs multi-chain correctly", {
  # Test filtering helper directly
  test_metadata <- data.frame(
    TRA_V = c("TRAV1-1", NA, "TRAV1-1", NA),
    TRA_J = c("TRAJ1-1", NA, "TRAJ1-1", NA),
    TRA = c("CAAA", NA, "CBBB", NA),
    TRB_V = c("TRBV1-1", "TRBV1-1", NA, NA),
    TRB_J = c("TRBJ1-1", "TRBJ1-1", NA, NA),
    TRB = c("CASSX", "CASSY", NA, NA),
    row.names = paste0("clone_", 1:4)
  )
  
  # Single-chain TRA: should keep clones 1, 3 (have complete TRA)
  keep_tra <- .filter_clones_for_assay(test_metadata, c("TRA"), verbose = FALSE)
  expect_equal(sum(keep_tra), 2)
  expect_equal(unname(which(keep_tra)), c(1, 3))
  
  # Single-chain TRB: should keep clones 1, 2 (have complete TRB)
  keep_trb <- .filter_clones_for_assay(test_metadata, c("TRB"), verbose = FALSE)
  expect_equal(sum(keep_trb), 2)
  expect_equal(unname(which(keep_trb)), c(1, 2))
  
  # Multi-chain TRA+TRB: should keep clones 1,2,3 (at least one chain complete)
  keep_multi <- .filter_clones_for_assay(test_metadata, c("TRA", "TRB"), verbose = FALSE)
  expect_equal(sum(keep_multi), 3)
  expect_equal(unname(which(keep_multi)), c(1, 2, 3))
})

test_that("RunTcrClustering respects outputPrefix parameter", {
  skip_if_not(nzchar(Sys.getenv("RETICULATE_PYTHON")) || nzchar(Sys.which("python3")) || nzchar(Sys.which("python")), "Python not available")
  
  set.seed(789)
  
  metadata <- data.frame(
    TRB_V = c("TRBV7-9", "TRBV6-4"),
    TRB_J = c("TRBJ2-1", "TRBJ1-1"),
    TRB = c("CASSA", "CASSB"),
    SubjectId = c("S1", "S2")
  )
  rownames(metadata) <- paste0("clone_", seq_len(nrow(metadata)))
  
  dummy_counts <- matrix(rnorm(nrow(metadata) * 10), nrow = 10, ncol = nrow(metadata),
                         dimnames = list(paste0("gene_", 1:10), rownames(metadata)))
  seurat_obj <- CreateSeuratObject(counts = dummy_counts, meta.data = metadata)
  
  tcr_obj <- RunTcrdist3(seurat_obj, chains = "TRB",
                         calculateChainPairs = FALSE,
                         minimumCloneSize = 1,
                         verbose = FALSE,
                         )
  
  output_dir <- file.path(tempdir(), "test_prefix")
  unlink(output_dir, recursive = TRUE)
  dir.create(output_dir, recursive = TRUE)
  
  results <- RunTcrClustering(
    seuratObj = tcr_obj,
    assays = "TRB",
    clusteringMethod = "DIANA",
    dianaHeight = 20,
    clusterSizeThreshold = 1,
    outputDir = output_dir,
    outputPrefix = "custom_prefix",
    verbose = FALSE
  )
  
  filename <- basename(results$parquet_files[1])
  expect_true(grepl("^custom_prefix_", filename))
})

test_that(".is_concatenated detects concatenation correctly", {
  expect_true(.is_concatenated("TRAV1-1,TRAV2-1"))
  expect_true(.is_concatenated("CASS;CAST"))
  expect_true(.is_concatenated("A,B,C"))
  expect_false(.is_concatenated("TRAV1-1"))
  expect_false(.is_concatenated("CASSLGQAYEQYF"))
  expect_false(.is_concatenated(NA_character_))
})

test_that(".strip_allele_suffix strips alleles correctly", {
  expect_equal(.strip_allele_suffix("TRBV7-9*01"), "TRBV7-9")
  expect_equal(.strip_allele_suffix("TRAJ1-1*02"), "TRAJ1-1")
  expect_equal(.strip_allele_suffix("TRBV7-9"), "TRBV7-9")  # No change if no allele
  expect_true(is.na(.strip_allele_suffix(NA_character_)))
})
