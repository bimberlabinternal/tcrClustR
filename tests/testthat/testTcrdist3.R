library(testthat)

test_that("tcrdist3 works", {
  #define paths using a temporary directory
  temp_dir <- tempdir()
  print(temp_dir)

  paste0('file.exists: ', file.exists(temp_dir))
  paste0('dir.exists: ', dir.exists(temp_dir))


  outFile <- tempfile(tmpdir =  temp_dir)
  print(paste0('outfile: ', outFile))
  paste0('file.exists: ', file.exists(outFile))
  print('creating')
  print(file.create(outFile))
  print(paste0('file.exists after create: ', file.exists(outFile)))

  postFormattingMetadataCsvPath <- file.path(temp_dir, "tcrdist3Input.csv")
  rdsOutputPath <- file.path(temp_dir, "tcrdist3DistanceMatrices")
  filteredGeneSegmentsPath <- file.path(temp_dir, "filtered_TRB_gene_segments.csv")

  seuratObj <- readRDS("../testdata/small_RIRA.rds")
  seuratObj <- subset(seuratObj, cells = SeuratObject::WhichCells(seuratObj, which(as.numeric(seuratObj$cDNA_ID) > 1 )))

  #test that the function runs without errors
  testthat::expect_no_error(
    seuratObj_TCR <- RunTcrdist3(seuratObj = seuratObj,
             metadata = NULL,
             formatMetadata = T,
             postFormattingMetadataCsvPath = postFormattingMetadataCsvPath,
             chains = c("TRA", "TRB"),
             cleanMetadata = T,
             minimumClonesPerSubject = 2,
             rdsOutputPath = rdsOutputPath,
             pythonExecutable = Sys.which("python3"),
             debugTcrdist3 = "True")
  )
  print(postFormattingMetadataCsvPath)
  print(list.files(temp_dir))
  #test that FormatMetadataForTcrDist3 worked
  testthat::expect_true(file.exists(postFormattingMetadataCsvPath))

  #test that the "missing TCRs file" was created and properly stores the TCRs missing from the db
  testthat::expect_true(file.exists(filteredGeneSegmentsPath))
  testthat::expect_gt(file.size(filteredGeneSegmentsPath), 16)
  testthat::expect_equal(file.size(filteredGeneSegmentsPath), 298)

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
    seuratObj_TCR <- RunTcrdist3(seuratObj = seuratObj,
             metadata = NULL,
             formatMetadata = T,
             postFormattingMetadataCsvPath = postFormattingMetadataCsvPath,
             chains = c("TRA", "TRB"),
             cleanMetadata = T,
             minimumClonesPerSubject = 2,
             rdsOutputPath = rdsOutputPath,
             pythonExecutable = system("which python3"),
             debugTcrdist3 = "True",
             spikeInDataframe = spikeInDataframe)
  )
  #read the resulting tcrdist3Input and ensure the spike-ins are present
  tcrdist3Input <- readr::read_csv(postFormattingMetadataCsvPath)
  testthat::expect_true(sum(grepl("spikeIn", tcrdist3Input$subject)) == 3)
  #test that the RDS distance matrices were created
  testthat::expect_true(file.exists(file.path(rdsOutputPath, "pw_cdr3_a_aa.rds")))

})
