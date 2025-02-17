context("conga functionality")

test_that("conga works", {
  seuratObj <- readRDS("../testdata/small_RIRA.rds")

  #test that the function runs without errors
  testthat::expect_no_error(
    RunConga(seuratObj = seuratObj,
             metadata = NULL,
             formatMetadata = T,
             postFormattingMetadataCsvPath = './congaInput.csv',
             chains = c("TRA", "TRB"),
             cleanMetadata = T,
             minimumClonesPerSubject = 2,
             rdsOutputPath = "./congaDistanceMatrices/",
             pythonExecutable = reticulate::py_exe(),
    )
  )

  #test that FormatMetadataForTcrDist3 worked
  testthat::expect_true(file.exists("./congaInput.csv"))

  #test that the "missing TCRs file" was created and properly stores the TCRs missing from the db
  testthat::expect_true(file.exists("./filtered_TRB_gene_segments.csv"))
  testthat::expect_gt(file.size("./filtered_TRB_gene_segments.csv"), 16)
  testthat::expect_equal(file.size("./filtered_TRB_gene_segments.csv"), 220)

  #test that the RDS distance matrices were created
  testthat::expect_true(file.exists("./congaDistanceMatrices/congaTcrDistances_TRA.rds"))


  spikeInDataframe <- data.frame(CloneNames = rep(1:3),
                                 TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
                                 TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
                                 TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
                                 TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
                                 TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
                                 TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY")
  )

  #test that spiking in TCRs works:
  testthat::expect_no_error(
    RunConga(seuratObj = seuratObj,
             metadata = NULL,
             formatMetadata = T,
             postFormattingMetadataCsvPath = './congaInput.csv',
             chains = c("TRA", "TRB"),
             cleanMetadata = T,
             minimumClonesPerSubject = 2,
             rdsOutputPath = "./congaDistanceMatrices/",
             pythonExecutable = reticulate::py_exe(),
             spikeInDataframe = spikeInDataframe
    )
  )
  #read the resulting congaInput and ensure the spike-ins are present
  #the conga input is a bit shambolic due to needing to be parsed literally by python,
  #so we just check that the spike-in sequences are present
  congaInput <-  readr::read_file("./congaInput.csv")
  testthat::expect_true(grepl("CASSAAAAAAAAFF", congaInput) &
                           grepl("CASSVVVVVVVVQF", congaInput) &
                           grepl("CASSWWWWWWWWQY", congaInput))
  #test that the RDS distance matrices were created
  testthat::expect_true(file.exists("./congaDistanceMatrices/congaTcrDistances_TRA.rds"))
})
