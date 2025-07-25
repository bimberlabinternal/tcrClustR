
#' @title RunConga
#' @description Run Conga's implementation of tcrdist on a Seurat object/metadata
#'
#' @param seuratObj Seurat object containing TCR information. If NULL, metadata must be provided.
#' @param metadata Data frame containing metadata. If NULL, seuratObj must be provided.
#' @param organism Organism to use for Conga. Default is 'human'.
#' @param formatMetadata Boolean controlling whether to format the metadata for Conga using the internal FormatMetadataForConga function. Default is TRUE.
#' @param postFormattingMetadataCsvPath Path to the output CSV file from FormatMetadataForConga. Default is './congaInput.csv'.
#' @param chains Vector of TCR chains to include in the analysis. Default is c("TRA", "TRB").
#' @param cleanMetadata Pass-through boolean controlling whether to clean the metadata by removing rows with NA values or commas in the specified chains. Default is TRUE.
#' @param minimumClonesPerSubject Minimum number of clones per subject to include in the analysis. Default is 2.
#' @param rdsOutputPath Path to the output directory for the RDS files containing the distance matrices. Default is "./tcrdist3DistanceMatrices/".
#' @param pythonExecutable Path to the python executable. Default is NULL, but imputes to reticulate::py_exe().
#' @param spikeInDataframe Data frame containing known CDR3s and gene segments to be included in the clustering. Default is NULL.
#'
#' @examples
#' \dontrun{
#'     RunConga(seuratObj = seuratObj,
#'              metadata = NULL,
#'              formatMetadata = T,
#'              postFormattingMetadataCsvPath = './congaInput.csv',
#'              chains = c("TRA", "TRB"),
#'              cleanMetadata = T,
#'              minimumClonesPerSubject = 2,
#'              rdsOutputPath = "./tcrdist3DistanceMatrices/",
#'              pythonExecutable = reticulate::py_exe()
#'              )
#'     spikeInDataframe <- data.frame(CloneNames = rep(1:3),
#'                                  TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
#'                                  TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
#'                                  TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
#'                                  TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
#'                                  TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
#'                                  TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY")
#' }
#'
#' @export

#TODO: flesh out examples demonstrating formatting requirements for spikeInDataframe

RunConga <- function(seuratObj = NULL,
                     metadata = NULL,
                     organism = 'human',
                     formatMetadata = T,
                     postFormattingMetadataCsvPath = './congaInput.csv',
                     chains = c("TRA", "TRB"),
                     cleanMetadata = T,
                     spikeInDataframe = NULL,
                     minimumClonesPerSubject = 2,
                     rdsOutputPath = "./congaDistanceMatrices/",
                     pythonExecutable = NULL
                     ) {
  #identify the metadata dataframe
  if (is.null(seuratObj) & is.null(metadata)) {
    stop("Please provide either a Seurat Object or the Seurat Object's metadata as input.")
  }
  #TODO: add a check for the metadata dataframe, or, if it's a csv file, read it.
  if (!is.null(seuratObj)) {
    metadata <- seuratObj@meta.data
  }
  if (is.null(pythonExecutable)) {
    pythonExecutable <- reticulate::py_exe()
  }

  for (chain in chains) {
    #call the formatting function to write the metadata to a csv file
    #see FormatMetadataForConga in congaUtils.R for details
    #TODO: potentially make this write to a tempfile?
    FormatMetadataForConga(metadata,
                           outputCsv = postFormattingMetadataCsvPath,
                           organism = organism,
                           chains = chain,
                           cleanMetadata = cleanMetadata,
                           minimumClonesPerSubject = minimumClonesPerSubject,
                           spikeInDataframe = spikeInDataframe)

    #run conga in python and return individual RDS files with the results (they may get very large)

    #read the python script template from tcrClustR and write them to a tempfile
    template <- readr::read_file(system.file("scripts/congaTcrDistances.py", package = "tcrClustR"))
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
    command <- paste0("getTcrDistances(tcrFile = '", postFormattingMetadataCsvPath,
                      "', organism = '", organism,
                      "', chain = '", chain,
                      "', rds_output_path = '", rdsOutputPath,
                      "')")
    readr::write_file(command, script, append = TRUE)
    #execute
    system2(pythonExecutable, script)
  }
  #return a seurat object, with the distance matrices implemented as assays
  for (chain in chains) {
    #read the RDS file
    rdsFile <- paste0(rdsOutputPath, "/congaTcrDistances_", chain, ".rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("RDS file not found: ", rdsFile))
    }
    distanceMatrix <- readRDS(rdsFile)
    colnames(distanceMatrix) <- paste0(chain, "_", seq_along(1:ncol(distanceMatrix)))
    #add the distance matrix to the seurat object
    if (is.null(seuratObj_TCR)){
      seuratObj_TCR <- SeuratObject::CreateSeuratObject(counts = distanceMatrix,
                                                        assay =  chain)
    } else {
      seuratObj_TCR_subsequentChain <- SeuratObject::CreateSeuratObject(counts = distanceMatrix, assay = chain)
      seuratObj_TCR <- merge(seuratObj_TCR, seuratObj_TCR_subsequentChain)
    }
  }
  return(seuratObj_TCR)
}
