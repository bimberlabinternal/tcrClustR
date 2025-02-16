context("tcrdist3 functionality")

test_that("tcrdist3 works", {

  seuratObj <- readRDS("../testdata/small_RIRA.rds")

  #test that the function runs without errors
  testthat::expect_no_error(
    RunTcrdist3(seuratObj = seuratObj,
             metadata = NULL,
             formatMetadata = T,
             postFormattingMetadataCsvPath = './tcrdist3Input.csv',
             chains = c("TRA", "TRB"),
             cleanMetadata = T,
             minimumClonesPerSubject = 2,
             rdsOutputPath = "./tcrdist3DistanceMatrices/",
             pythonExecutable = reticulate::py_exe(),
             debugTcrdist3 = "True")
  )

  #test that FormatMetadataForTcrDist3 worked
  testthat::expect_true(file.exists("./tcrdist3Input.csv"))

  #test that the "missing TCRs file" was created and properly stores the TCRs missing from the db
  testthat::expect_true(file.exists("./filtered_TRB_gene_segments.csv"))
  testthat::expect_gt(file.size("./filtered_TRB_gene_segments.csv"), 16)
  testthat::expect_equal(file.size("./filtered_TRB_gene_segments.csv"), 298)

  #test that the RDS distance matrices were created
  testthat::expect_true(file.exists("./tcrdist3DistanceMatrices/pw_cdr3_a_aa.rds"))

})
