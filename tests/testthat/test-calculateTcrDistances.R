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
           debugTcrdist3 = TRUE
  )


  testthat::expect_equal(length(names(seuratObj_TCR@misc$TCR_Distances)), 4)
  testthat::expect_equal(length(names(seuratObj_TCR@assays)), 4)
  
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
           debugTcrdist3 = TRUE,
           spikeInDataframe = spikeInDataframe
  )

  testthat::expect_equal(sum(grepl("SpikeIn", seuratObj_TCR$SubjectId)), 3)

  seuratObj_TCR_multichain <- CalculateTcrDistances(
           inputData = seuratObj,
           chains = c("TRA", "TRB"),
           minimumCloneSize = 2,
           debugTcrdist3 = TRUE,
           calculateChainPairs = TRUE,
           verbose = TRUE
  )
  
  testthat::expect_equal(SeuratObject::Assays(seuratObj_TCR_multichain), c("TRA_fl", "TRA_cdr3", "TRB_fl", "TRB_cdr3", "TRA.TRB_fl", "TRA.TRB_cdr3"))
  
  #test that joint distance matrices exist and have reasonable properties
  joint_matrix <- SeuratObject::GetAssayData(seuratObj_TCR_multichain, assay = "TRA.TRB_fl", layer = "counts")
  tra_matrix <- SeuratObject::GetAssayData(seuratObj_TCR_multichain, assay = "TRA_fl", layer = "counts")
  trb_matrix <- SeuratObject::GetAssayData(seuratObj_TCR_multichain, assay = "TRB_fl", layer = "counts")
  
  #joint matrix should exist and be non-empty
  testthat::expect_gt(nrow(joint_matrix), 0)
  testthat::expect_gt(ncol(joint_matrix), 0)
  
  #TODO: expand test cases
})
