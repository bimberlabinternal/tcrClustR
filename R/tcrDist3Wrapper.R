#' @title RunTcrdist3
#' @description This function runs the tcrdist3 pipeline on a set of CDR3 sequences in a Seurat Object.
#'
#' @param seuratObj Seurat Object containing TCR information. If NULL, metadata must be provided.
#' @param metadata Data frame containing metadata. If NULL, seuratObj must be provided.
#' @param organism Organism to use for tcrdist3. Default is 'human'.
#' @param formatMetadata Boolean controlling whether to format the metadata for tcrdist3 using the internal FormatMetadataForTcrDist3 function. Default is TRUE.
#' @param postFormattingMetadataCsvPath Path to the output CSV file from FormatMetadataForTcrDist3. Default is './tcrDist3Input.csv'.
#' @param chains Vector of TCR chains to include in the analysis. Default is c("TRA", "TRB").
#' @param cleanMetadata Pass-through boolean controlling whether to clean the metadata by removing rows with NA values or commas in the specified chains. Default is TRUE.
#' @param summarizeClones Pass-through boolean controlling whether to summarize clones into a frequency (by SubjectId, TRA, TRB, TRA_V, and TRB_J). Default is TRUE.
#' @param imputeCloneNames Pass-through boolean controlling whether to impute clone names if they are missing. Existing clone names will be inherited. Default is TRUE.
#' @param minimumClonesPerSubject Minimum number of clones per subject to include in the analysis. Default is 2.
#' @param rdsOutputPath Path to the output directory for the RDS files containing the distance matrices. Default is "./tcrdist3DistanceMatrices/".
#' @param pythonExecutable Path to the python executable. Default is reticulate::py_exe().
#' @param debugTcrdist3 String (to be passed to python and converted to boolean) controlling whether to run tcrdist3 in debug mode. Default is "True".
#'
#'@examples
#'\dontrun{
#'   RunTcrdist3(seuratObj = seuratObj,
#'               metadata = NULL,
#'               formatMetadata = T,
#'               postFormattingMetadataCsvPath = './tcrDist3Input.csv',
#'               chainsString = "alpha,beta",
#'               cleanMetadata = T,
#'               summarizeClones = T,
#'               imputeCloneNames = T,
#'               minimumClonesPerSubject = 2,
#'               rdsOutputPath = "./tcrdist3DistanceMatrices/",
#'               pythonExecutable = reticulate::py_exe(),
#'               debugTcrdist3 = "True")
#' }
#'
#'
#'


RunTcrdist3 <- function(seuratObj = NULL,
                        metadata = NULL,
                        organism = 'human',
                        formatMetadata = T,
                        postFormattingMetadataCsvPath = './tcrDist3Input.csv',
                        chains = c("TRA", "TRB"),
                        cleanMetadata = T,
                        summarizeClones = T,
                        imputeCloneNames = T,
                        minimumClonesPerSubject = 2,
                        rdsOutputPath = "./tcrdist3DistanceMatrices/",
                        pythonExecutable = reticulate::py_exe(),
                        debugTcrdist3 = "True") {
  #identify the metadata dataframe
  if (is.null(seuratObj) & is.null(metadata)) {
    stop("Please provide either a Seurat Object or the Seurat Object's metadata as input.")
  }
  #TODO: add a check for the metadata dataframe, or, if it's a csv file, read it.
  if (!is.null(seuratObj)) {
    metadata <- seuratObj@meta.data
  }

  #format metadata if necessary (a user may have done this already, so add optional flag)
  if (formatMetadata) {
    metadata <- FormatMetadataForTcrDist3(metadata = metadata,
                                          outputCsv = postFormattingMetadataCsvPath,
                                          chains = chains,
                                          cleanMetadata = cleanMetadata,
                                          summarizeClones = summarizeClones,
                                          imputeCloneNames = imputeCloneNames,
                                          minimumClonesPerSubject = minimumClonesPerSubject)

  }

  #convert chains into a string for parsing in python
  chainsString <- ""
  for (chain in chains) {
    if (chain == "TRA") {
      chainsString <- paste0(chainsString, "alpha")
      } else if (chain == "TRB") {
        chainsString <- paste0(chainsString, "beta")
      } else if (chain == "TRG") {
        chainsString <- paste0(chainsString, "gamma")
      } else if (chain == "TRD") {
        chainString <- paste0(chainsString, "delta")
      } else {
        warning(paste0("Chain Type ", chain, " not recognized. Skipping."))
      }
    }

  #run tcrdist3 in python and return individual RDS files with the results (they may get very large)

  #read the python script template from tcrClustR and write them to a tempfile
  template <- readr::read_file(system.file("scripts/tcrdist3TcrDistances.py", package = "tcrClustR"))
  script <- tempfile()
  readr::write_file(template, script)

  #convert paths to absolute paths
  postFormattingMetadataCsvPath <- R.utils::getAbsolutePath(postFormattingMetadataCsvPath)
  rdsOutputPath <- R.utils::getAbsolutePath(rdsOutputPath)
  pythonExecutable <- R.utils::getAbsolutePath(pythonExecutable)

  #create distance matrix output directory if it doesn't exist
  if (!dir.exists(rdsOutputPath)) {
    dir.create(rdsOutputPath)
  }

  #format and write the python function to the end of the script
  command <- paste0("writeTcrDistances(csv_path = '", postFormattingMetadataCsvPath,
                   "', organism = '", organism,
                   "', chainsString = '", chainsString,
                   "', db_file = 'alphabeta_gammadelta_db.tsv", #I don't think we can change the db file, but perhaps we'd want to someday?
                   "', rds_output_path = '", rdsOutputPath,
                   "', debug ='", debugTcrdist3,
                   "')")
  readr::write_file(command, script, append = TRUE)
  #execute
  system2(pythonExecutable, script)
}






