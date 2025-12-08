library(testthat)

test_that("CalculateTcrDistances works", {
  #read in a small Seurat object with TCR data
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1 )))

  #test that the function runs without errors
  seuratObj_TCR <- CalculateTcrDistances(
           inputData = seuratObj,
           chains = c("TRA", "TRB"),
           minimumCloneSize = 2,
           debugMode = TRUE
  )

  testthat::expect_equal(max(seuratObj_TCR$TRA_CloneSize), 1045)
  testthat::expect_equal(max(seuratObj_TCR$TRA_CloneSize), 1045)
  testthat::expect_equal(length(names(seuratObj_TCR@misc$TCR_Distances)), 4)
  testthat::expect_equal(length(names(seuratObj_TCR@assays)), 1)
  testthat::expect_equal(ncol(seuratObj_TCR@assays$RNA), ncol(seuratObj@assays$RNA))
  testthat::expect_true('umap' %in% names(seuratObj_TCR@reductions))
  
  spikeInDataframe <- data.frame(CloneNames = rep(1:3),
                                 TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
                                 TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
                                 TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
                                 TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
                                 TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
                                 TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY")
  )

  # Test that spiking in TCRs works:
  seuratObj_TCR <- CalculateTcrDistances(
           inputData = seuratObj,
           chains = c("TRA", "TRB"),
           minimumCloneSize = 2,
           debugMode = TRUE,
           spikeInDataframe = spikeInDataframe
  )

  testthat::expect_true(!'RNA' %in% names(seuratObj_TCR@assays))
  testthat::expect_equal(length(names(seuratObj_TCR@assays)), 1)
  testthat::expect_equal(sum(!is.na(seuratObj_TCR$CloneNames)), 3)
  testthat::expect_equal(sum(seuratObj_TCR$IsSpikeInClone), 3)
  testthat::expect_equal(length(names(seuratObj_TCR@misc$TCR_Distances)), 4)

  seuratObj_TCR_multichain <- CalculateTcrDistances(
           inputData = seuratObj,
           chains = c("TRA", "TRB"),
           minimumCloneSize = 2,
           debugMode = TRUE,
           calculateChainPairs = TRUE,
           verbose = TRUE
  )
  
  testthat::expect_equal(SeuratObject::Assays(seuratObj_TCR_multichain), c("RNA"))
  testthat::expect_equal(names(seuratObj_TCR_multichain@misc$TCR_Distances), c("TRA_fl", "TRA_cdr3", "TRB_fl", "TRB_cdr3", "TRA_TRB_fl", "TRA_TRB_cdr3"))
  
  joint_matrix <- GetDistanceMatrix(seuratObj_TCR_multichain, chains = c("TRA", "TRB"))
  testthat::expect_equal(nrow(joint_matrix), 21)
  testthat::expect_equal(ncol(joint_matrix), 21)
  
  tra_matrix <- GetDistanceMatrix(seuratObj_TCR_multichain, chains = c("TRA"))
  testthat::expect_equal(nrow(tra_matrix), 28)
  testthat::expect_equal(ncol(tra_matrix), 28)
  
  trb_matrix <- GetDistanceMatrix(seuratObj_TCR_multichain, chains = c("TRB"))
  testthat::expect_equal(nrow(tra_matrix), 28)
  testthat::expect_equal(ncol(tra_matrix), 28)
})

# test_that("CalculateTcrDistances works for rhesus", {
#   #read in a small Seurat object with TCR data
#   seuratObj <- readRDS("../testdata/small_RIRA.rds")
#   seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1 )))
#
#   #test that the function runs without errors
#   seuratObj_TCR <- CalculateTcrDistances(
#     inputData = seuratObj,
#     chains = c("TRA", "TRB"),
#     minimumCloneSize = 2,
#     organism = 'rhesus',
#     debugMode = TRUE
#   )
#
#   testthat::expect_equal(length(names(seuratObj_TCR@misc$TCR_Distances)), 4)
#   testthat::expect_equal(length(names(seuratObj_TCR@assays)), 1)
#   testthat::expect_equal(ncol(seuratObj_TCR@assays$RNA), ncol(seuratObj@assays$RNA))
#   testthat::expect_true('umap' %in% names(seuratObj_TCR@reductions))
# })