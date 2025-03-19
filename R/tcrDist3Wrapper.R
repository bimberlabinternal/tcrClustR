
utils::globalVariables(
  names = c('SubjectId', 'TRA_V', 'TRA_J', 'TRB_V', 'TRB_J', 'CloneNames', 'count'),
  package = 'tcrClustR',
  add = TRUE
)

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
#' @param spikeInDataframe Data frame containing known CDR3s and gene segments to be included in the clustering. Default is NULL.
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
#'
#'     spikeInDataframe <- data.frame(CloneNames = rep(1:3),
#'                                  TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
#'                                  TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
#'                                  TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
#'                                  TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
#'                                  TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
#'                                  TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY"))
#' }
#'
#' @export

#TODO: flesh out examples demonstrating formatting requirements for spikeInDataframe


RunTcrdist3 <- function(seuratObj = NULL,
                        metadata = NULL,
                        organism = 'human',
                        formatMetadata = T,
                        postFormattingMetadataCsvPath = './tcrDist3Input.csv',
                        chains = c("TRA", "TRB"),
                        cleanMetadata = T,
                        spikeInDataframe = NULL,
                        summarizeClones = T,
                        imputeCloneNames = T,
                        minimumClonesPerSubject = 2,
                        rdsOutputPath = "./tcrdist3DistanceMatrices/",
                        pythonExecutable = NULL,
                        debugTcrdist3 = "True") {
  #TODO: allow for more direct control of all of the files that will be written
  #FormatMetadata can write several different kinds of files (filtered gene segments, database)
  #a "metadata directory" is probably the cleanest way to organize those files.
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

  #format metadata if necessary (a user may have done this already, so add optional flag)
  if (formatMetadata) {
    formatted_metadata <- FormatMetadataForTcrDist3(metadata = metadata,
                                          outputCsv = postFormattingMetadataCsvPath,
                                          chains = chains,
                                          cleanMetadata = cleanMetadata,
                                          summarizeClones = summarizeClones,
                                          imputeCloneNames = imputeCloneNames,
                                          minimumClonesPerSubject = minimumClonesPerSubject,
                                          spikeInDataframe = spikeInDataframe)

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
  print(paste0("Creating tcrdist3 distance matrices in the following directory: ", rdsOutputPath))
  #format and write the python function to the end of the script
  command <- paste0("writeTcrDistances(csv_path = '", postFormattingMetadataCsvPath,
                   "', organism = '", organism,
                   "', chainsString = '", chainsString,
                   "', db_file = 'alphabeta_gammadelta_db.tsv", #TODO: I don't think we can change the db file, but perhaps we'd want to someday?
                   "', rds_output_path = '", rdsOutputPath,
                   "', debug ='", debugTcrdist3,
                   "')")
  readr::write_file(command, script, append = TRUE)
  #execute
  system2(pythonExecutable, script)

  #return a seurat object, with the distance matrices implemented as assays
  for (chain in chains) {
    #read the RDS file
    if (chain == "TRA") {
      chain_tcrdist3 <- "alpha"
    } else if (chain == "TRB") {
      chain_tcrdist3 <- "beta"
    } else if (chain == "TRG") {
      chain_tcrdist3 <- "gamma"
    } else if (chain == "TRD") {
      chain_tcrdist3 <- "delta"
    } else {
      warning(paste0("Chain Type ", chain, " not recognized. Skipping."))
    }
    #process the full length TCR distance matrix
    rdsFile <- paste0(rdsOutputPath, "/pw_", chain_tcrdist3, ".rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise 'full' tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_full_length <- readRDS(rdsFile)
    colnames(distanceMatrix_full_length) <- paste0(chain, "_", seq_along(1:ncol(distanceMatrix_full_length)))
    rownames(distanceMatrix_full_length) <- paste0(chain, "_", seq_along(1:nrow(distanceMatrix_full_length)))

    #process the CDR3 only TCR distance matrix
    #grab the first letter of the chain (distance matrices for the cdr3 are stored as "pw_cdr3_b_aa.rds" for a beta chain)
    chain_cdr3_id <- tolower(strsplit(chain_tcrdist3, split = "")[[1]][1])
    rdsFile <- paste0(rdsOutputPath, "/pw_cdr3_", chain_cdr3_id, "_aa.rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise CDR3 tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_CDR3 <- readRDS(rdsFile)
    colnames(distanceMatrix_CDR3) <- paste0(chain, "_", seq_along(1:ncol(distanceMatrix_CDR3)), "_cdr3")
    rownames(distanceMatrix_CDR3) <- paste0(chain, "_", seq_along(1:nrow(distanceMatrix_CDR3)), "_cdr3")

    #add the distance matrices to the Seurat object
    if (is.null(seuratObj_TCR)){
      seuratObj_TCR <- SeuratObject::CreateSeuratObject(counts = distanceMatrix_full_length,
                                                        assay =  chain)
      #TODO: linked TODO with L231 in tcrdistUtils.R, this currently only works for joint TRA+TRBs.
      seuratObj_TCR <- Seurat::AddMetaData(seuratObj_TCR, metadata = formatted_metadata)
      seuratObj_TCR_CDR3 <- SeuratObject::CreateSeuratObject(counts = distanceMatrix_CDR3, assay = paste0(chain, "_cdr3"))
      seuratObj_TCR_CDR3 <- Seurat::AddMetaData(seuratObj_TCR_CDR3, metadata = formatted_metadata)
      seuratObj_TCR <- merge(seuratObj_TCR, seuratObj_TCR_CDR3)
    } else {
      seuratObj_TCR_subsequentChain <- SeuratObject::CreateSeuratObject(counts = distanceMatrix_full_length,
                                                                        assay =  chain)
      seuratObj_TCR_subsequentChain <- Seurat::AddMetaData(seuratObj_TCR_subsequentChain, metadata = formatted_metadata)
      seuratObj_TCR_CDR3_subsequentChain <- SeuratObject::CreateSeuratObject(counts = distanceMatrix_CDR3, assay = paste0(chain, "_cdr3"))
      seuratObj_TCR_CDR3_subsequentChain <- Seurat::AddMetaData(seuratObj_TCR_CDR3_subsequentChain, metadata = formatted_metadata)
      seuratObj_TCR_subsequentChain <- merge(seuratObj_TCR_subsequentChain, seuratObj_TCR_CDR3_subsequentChain)
      seuratObj_TCR <- merge(seuratObj_TCR, seuratObj_TCR_subsequentChain)
    }
  }
  return(seuratObj_TCR)
}






