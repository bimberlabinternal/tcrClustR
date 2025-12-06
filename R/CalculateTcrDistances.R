#' @include Utils.R
#' @include RunTcrdist3.R
#' @include FormatMetadata.R
#'
#' @title CalculateTcrDistances
#' @description This function calculates pairwise TCR distances and stores them in a Seurat Object.
#'
#' @param inputData Either a seurat Object containing TCR information or a dataframe containing metadata
#' @param organism Organism to use for reference segments. Default is 'human'.
#' @param chains Vector of TCR chains to include in the analysis. Default is c("TRA", "TRB").
#' @param spikeInDataframe An optional data frame containing additional CDR3s and gene segments to be included in the clustering, such as published clonotypes.
#' @param minimumCloneSize Minimum number of clones per subject to include in the analysis. Default is 2.
#' @param calculateChainPairs Boolean controlling whether to compute joint multi-chain distance matrices for observed chain combinations.
#' @param debugMode Boolean that will result in more debugging information to be printed
#' @param verbose Boolean controlling whether to print additional information about procesing
#' @param method The method for clustering. Currently only tcrdist3 is supported
#' @import Matrix
#' @importFrom methods as
#' @importFrom methods is
#' @examples
#' \dontrun{
#'   resultList <- CalculateTcrDistances(inputData = seuratObj@meta.data,
#'               chains = c("TRA", "TRB"),
#'               minimumCloneSize = 2,
#'               calculateChainPairs = FALSE,
#'               verbose = FALSE)
#'
#'   # Example with calculateChainPairs analysis
#'   resultList <- CalculateTcrDistances(inputData = seuratObj,
#'               chains = c("TRA", "TRB", "TRG", "TRD"),
#'               calculateChainPairs = TRUE,
#'               verbose = TRUE)
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
CalculateTcrDistances <- function(inputData = NULL,
  organism = 'human',
  chains = c("TRA", "TRB"),
  spikeInDataframe = NULL,
  minimumCloneSize = 2,
  calculateChainPairs = FALSE,
  debugMode = FALSE,
  verbose = FALSE,
  method = 'tcrdist3'
) {
  if (is.null(inputData)) {
    stop("Please provide either a Seurat Object or the Seurat Object's metadata as input.")
  }

  if (typeof(inputData) == 'S4' && class(inputData)[1] == 'Seurat' ){
    metadata <- inputData@meta.data
  }

  if (!is.data.frame(metadata)) {
    stop('Expected inputData to be either a Seurat object or a data.frame')
  }

  if (method == 'tcrdist3') {
    return(RunTcrdist3(
      inputData = inputData,
      organism = organism,
      chains = chains,
      spikeInDataframe = spikeInDataframe,
      minimumCloneSize = minimumCloneSize,
      calculateChainPairs = calculateChainPairs,
      debugTcrdist3 = debugMode,
      verbose = verbose
    ))
  } else {
    stop(paste0('Unknown method: ', method))
  }
}