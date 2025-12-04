library(testthat)

test_that("tcrdist3 works", {
  #read in a small Seurat object with TCR data
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1 )))

  #test that the function runs without errors
  tcrOutputs <- RunTcrdist3(
           inputData = seuratObj,
           chains = c("TRA", "TRB"),
           minimumCloneSize = 2,
           debugTcrdist3 = TRUE
  )

  print(list.files(tempdir()))
  testthat::expect_equal(length(names(tcrOutputs)), 5)

  spikeInDataframe <- data.frame(CloneNames = rep(1:3),
                                 TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
                                 TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
                                 TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
                                 TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
                                 TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
                                 TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY")
  )

  # Test that spiking in TCRs works:
  seuratObj_TCR <- RunTcrdist3(
           inputData = seuratObj,
           chains = c("TRA", "TRB"),
           minimumCloneSize = 2,
           debugTcrdist3 = TRUE,
           spikeInDataframe = spikeInDataframe
  )

  testthat::expect_true(sum(grepl("spikeIn", seuratObj_TCR$SubjectId)) == 3)

  #test calculateChainPairs functionality
  seuratObj_TCR_multichain <- NULL
  testthat::expect_no_error(
    seuratObj_TCR_multichain <- RunTcrdist3(
             inputData = seuratObj,
             chains = c("TRA", "TRB"),
             minimumCloneSize = 2,
             debugTcrdist3 = TRUE,
             calculateChainPairs = TRUE,
             verbose = TRUE)
  )
  
  #test that calculateChainPairs assays were created
  assay_names <- SeuratObject::Assays(seuratObj_TCR_multichain)
  
  #we should have multiple joint combinations:
  # TRA_TRB (full TRA + full TRB)
  # TRA_TRB_cdr3 (full TRA + CDR3 TRB) 
  # TRA_cdr3_TRB (CDR3 TRA + full TRB)
  # TRA_cdr3_TRB_cdr3 (CDR3 TRA + CDR3 TRB)
  testthat::expect_true("TRA_TRB" %in% assay_names)
  testthat::expect_true("TRA_TRB_cdr3" %in% assay_names)
  testthat::expect_true("TRA_cdr3_TRB" %in% assay_names) 
  testthat::expect_true("TRA_cdr3_TRB_cdr3" %in% assay_names)
  
  #test that joint distance matrices exist and have reasonable properties
  joint_matrix <- SeuratObject::GetAssayData(seuratObj_TCR_multichain, assay = "TRA_TRB", layer = "counts")
  tra_matrix <- SeuratObject::GetAssayData(seuratObj_TCR_multichain, assay = "TRA", layer = "counts")
  trb_matrix <- SeuratObject::GetAssayData(seuratObj_TCR_multichain, assay = "TRB", layer = "counts")
  
  #joint matrix should exist and be non-empty
  testthat::expect_gt(nrow(joint_matrix), 0)
  testthat::expect_gt(ncol(joint_matrix), 0)
  
  #joint matrix should be square (distance matrix property)
  testthat::expect_equal(nrow(joint_matrix), ncol(joint_matrix))
  
  #joint matrix dimensions should be <= individual matrix dimensions 
  # (since it only includes observations with both chains)
  testthat::expect_lte(nrow(joint_matrix), nrow(tra_matrix))
  testthat::expect_lte(nrow(joint_matrix), nrow(trb_matrix))
  
  #all joint distances should be non-negative
  testthat::expect_true(all(joint_matrix >= 0))
  
  #test mixed CDR3/full combinations exist
  mixed_matrix1 <- SeuratObject::GetAssayData(seuratObj_TCR_multichain, assay = "TRA_TRB_cdr3", layer = "counts")
  mixed_matrix2 <- SeuratObject::GetAssayData(seuratObj_TCR_multichain, assay = "TRA_cdr3_TRB", layer = "counts")
  
  testthat::expect_gt(nrow(mixed_matrix1), 0)
  testthat::expect_gt(nrow(mixed_matrix2), 0)

})
