library(testthat)
library(Seurat)
library(arrow)

test_that("RunTcrClustering filters invalid clones correctly", {
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
                         dimnames = list(paste0("gene-", 1:10), rownames(metadata)))
  seuratObj <- CreateSeuratObject(counts = Seurat::as.sparse(dummy_counts), meta.data = metadata)
  
  seuratObj <- CalculateTcrDistances(seuratObj, 
                              chains = "TRB",
                              calculateChainPairs = FALSE,
                              minimumCloneSize = 1,
                              verbose = FALSE
  )
  
  seuratObj_NoFilter <- RunTcrClustering(
    seuratObj = seuratObj,
    chainsToCluster = 'TRB',
    clusterSizeThreshold = 1,
    verbose = FALSE
  )
  expect_equal(length(unique(seuratObj_NoFilter$TRB_fl_ClusterIdx)), 3)
  expect_equal(length(unique(seuratObj_NoFilter$TRB_cdr3_ClusterIdx)), 3)
  expect_equal(max(seuratObj_NoFilter$TRB_CloneSize), 1)
  
  seuratObj_Filter <- RunTcrClustering(
    seuratObj = seuratObj,
    chainsToCluster = 'TRB',
    clusterSizeThreshold = 5,
    verbose = FALSE
  )
  expect_equal(length(unique(seuratObj_Filter$TRB_fl_ClusterIdx)), 2)
  expect_equal(length(unique(seuratObj_Filter$TRB_cdr3_ClusterIdx)), 2)
})