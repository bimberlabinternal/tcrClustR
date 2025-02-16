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

})
