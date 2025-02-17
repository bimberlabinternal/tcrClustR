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
    RunTcrdist3(seuratObj = seuratObj,
             metadata = NULL,
             formatMetadata = T,
             postFormattingMetadataCsvPath = './tcrdist3Input.csv',
             chains = c("TRA", "TRB"),
             cleanMetadata = T,
             minimumClonesPerSubject = 2,
             rdsOutputPath = "./tcrdist3DistanceMatrices/",
             pythonExecutable = reticulate::py_exe(),
             debugTcrdist3 = "True",
             spikeInDataframe = spikeInDataframe)
  )
  #read the resulting tcrdist3Input and ensure the spike-ins are present
  tcrdist3Input <- readr::read_csv("./tcrdist3Input.csv")
  testthat::expect_true(sum(grepl("spikeIn", tcrdist3Input$subject)) == 3)
  #test that the RDS distance matrices were created
  testthat::expect_true(file.exists("./tcrdist3DistanceMatrices/pw_cdr3_a_aa.rds"))

})
