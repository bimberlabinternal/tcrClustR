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
#'   # When providing a dataframe, a new seurat object will be returned:
#'   seuratObj_TCR <- CalculateTcrDistances(inputData = seuratObj@meta.data,
#'               chains = c("TRA", "TRB"),
#'               minimumCloneSize = 2,
#'               calculateChainPairs = FALSE,
#'               verbose = FALSE)
#'
#'   # When providing a seuratObj, the original object will be returned, with the clustering results added as assays:
#'   seuratObj <- CalculateTcrDistances(inputData = seuratObj,
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

  isSeuratObj <- FALSE
  if (typeof(inputData) == 'S4' && class(inputData)[1] == 'Seurat' ){
    metadata <- inputData@meta.data
    isSeuratObj <- TRUE
  } else {
    metadata <- inputData
  }

  if (!is.data.frame(metadata)) {
    stop('Expected inputData to be either a Seurat object or a data.frame')
  }

  formatted_metadata <- .FormatMetadata(metadata = metadata,
                                        chains = chains,
                                        spikeInDataframe = spikeInDataframe,
                                        calculateChainPairs = calculateChainPairs,
                                        verbose = verbose
  )

  if (nrow(formatted_metadata[!formatted_metadata$IsSpikeInClone,]) != nrow(metadata)) {
    stop('Incorrect row count after .FormatMetadata')
  }

  if (any(rownames(formatted_metadata[!formatted_metadata$IsSpikeInClone,]) != rownames(metadata))) {
    stop('Rownames did not match after .FormatMetadata')
  }

  seuratObj <- NULL
  if (isSeuratObj && is.null(spikeInDataframe)) {
    if (verbose) {
      message('The results will be written to the original seurat object')
    }
    seuratObj <- inputData

    cols <- sort(c(
      grep(names(formatted_metadata), pattern = '_CloneIdx$', value = TRUE),
      grep(names(formatted_metadata), pattern = '_CloneSize$', value = TRUE),
      grep(names(formatted_metadata), pattern = '_ValidForClustering$', value = TRUE)
    ))
    toAdd <- formatted_metadata[cols]
    if (any(rownames(toAdd) != rownames(seuratObj@meta.data))) {
      stop('Row names did not match when applying metadata')
    }

    seuratObj <- Seurat::AddMetaData(seuratObj, toAdd)
  } else {
    if (isSeuratObj && !is.null(spikeInDataframe) && verbose) {
        message('The results will be written to a new object because spike-in data was used')
    }

    # NOTE: we need to create a dummy GEX assay, which we can drop later:
    dummy_mat <- matrix(rep(0, nrow(formatted_metadata)*2), nrow = 2, dimnames = list(c(LETTERS[1:2]), rownames(formatted_metadata)))
    seuratObj <- Seurat::CreateSeuratObject(counts = Seurat::as.sparse(dummy_mat), assay = '~PLACEHOLDER~', meta.data = formatted_metadata)
  }

  if ('TCR_Distances' %in% names(seuratObj@misc)) {
    seuratObj@misc$TCR_Distances <- list()
  }

  if (method == 'tcrdist3') {
    seuratObj <- RunTcrdist3(
      seuratObj = seuratObj,
      organism = organism,
      chains = chains,
      minimumCloneSize = minimumCloneSize,
      calculateChainPairs = calculateChainPairs,
      debugTcrdist3 = debugMode,
      verbose = verbose
    )
  } else {
    stop(paste0('Unknown method: ', method))
  }

  if ('~PLACEHOLDER~' %in% names(seuratObj@assays)) {
    if (length(names(seuratObj@assays)) > 1) {
      seuratObj[['~PLACEHOLDER~']] <- NULL
    } else {
      message('Seurat objects require one assay, so will not delete the PLACEHOLDER counts')
    }
  }

  return(seuratObj)
}