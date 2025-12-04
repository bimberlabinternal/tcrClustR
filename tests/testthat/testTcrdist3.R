library(testthat)

test_that("tcrdist3 works", {
  #define paths using a temporary directory
  temp_dir <- tempdir()

  print(temp_dir)
  #debug statements
  paste0('file.exists: ', file.exists(temp_dir))
  paste0('dir.exists: ', dir.exists(temp_dir))
  outFile <- tempfile(tmpdir =  temp_dir)
  print(paste0('outfile: ', outFile))
  paste0('file.exists: ', file.exists(outFile))
  print('creating')
  print(file.create(outFile))
  print(paste0('file.exists after create: ', file.exists(outFile)))
  print(paste0('python executable: ', Sys.which("python3")))

  rdsOutputPath <- file.path(temp_dir, "tcrdist3DistanceMatrices")
  filteredGeneSegmentsPath <- file.path(temp_dir, "filtered_TRB_gene_segments.csv")

  #read in a small Seurat object with TCR data
  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1 )))

  #test that the function runs without errors
  testthat::expect_no_error(
    seuratObj_TCR <- RunTcrdist3(
             inputData = seuratObj,
             chains = c("TRA", "TRB"),
             minimumCloneSize = 2,
             rdsOutputPath = rdsOutputPath,
             debugTcrdist3 = TRUE)
  )

  print(list.files(temp_dir))

  #test that the "missing TCRs file" was created and properly stores the TCRs missing from the db
  testthat::expect_true(file.exists(filteredGeneSegmentsPath))
  testthat::expect_gt(file.size(filteredGeneSegmentsPath), 16)
  
  #validate files
  filtered_content <- readr::read_csv(filteredGeneSegmentsPath, show_col_types = FALSE)
  testthat::expect_gt(nrow(filtered_content), 0)  # Should have some filtered gene segments
  testthat::expect_true("TRB_V" %in% colnames(filtered_content))  # Should have expected columns
  testthat::expect_true("TRB_J" %in% colnames(filtered_content))

  #test that the RDS distance matrices were created
  testthat::expect_true(file.exists(file.path(rdsOutputPath, "pw_cdr3_a_aa.rds")))


  spikeInDataframe <- data.frame(CloneNames = rep(1:3),
                                 TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
                                 TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
                                 TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
                                 TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
                                 TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
                                 TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY")
  )

  #test that spiking in TCRs works:
  seuratObj_TCR <- NULL
  testthat::expect_no_error(
    seuratObj_TCR <- RunTcrdist3(
             inputData = seuratObj,
             chains = c("TRA", "TRB"),
             minimumCloneSize = 2,
             rdsOutputPath = rdsOutputPath,
             debugTcrdist3 = TRUE,
             spikeInDataframe = spikeInDataframe)
  )

  testthat::expect_true(sum(grepl("spikeIn", seuratObj_TCR$SubjectId)) == 3)
  testthat::expect_true(file.exists(file.path(rdsOutputPath, "pw_cdr3_a_aa.rds")))

  #test calculateChainPairs functionality
  seuratObj_TCR_multichain <- NULL
  testthat::expect_no_error(
    seuratObj_TCR_multichain <- RunTcrdist3(
             inputData = seuratObj,
             chains = c("TRA", "TRB"),
             minimumCloneSize = 2,
             rdsOutputPath = rdsOutputPath,
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
