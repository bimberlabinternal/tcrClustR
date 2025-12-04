
utils::globalVariables(
  names = c('SubjectId', 'TRA_V', 'TRA_J', 'TRB_V', 'TRB_J', 'CloneNames', 'count'),
  package = 'tcrClustR',
  add = TRUE
)

#' @title RunTcrdist3
#' @description This function runs the tcrdist3 pipeline on a set of CDR3 sequences in a Seurat Object.
#'
#' @param inputData Either a seurat Object containing TCR information or a dataframe containing metadata
#' @param organism Organism to use for tcrdist3. Default is 'human'.
#' @param chains Vector of TCR chains to include in the analysis. Default is c("TRA", "TRB").
#' @param spikeInDataframe Data frame containing known CDR3s and gene segments to be included in the clustering. Default is NULL.
#' @param minimumCloneSize Minimum number of clones per subject to include in the analysis. Default is 2.
#' @param calculateChainPairs Boolean controlling whether to compute joint multi-chain distance matrices for observed chain combinations. Default is FALSE.
#' @param rdsOutputPath Path to the output directory for the RDS files containing the distance matrices. Default is "./tcrdist3DistanceMatrices/".
#' @param pythonExecutable Path to the python executable. Default is reticulate::py_exe().
#' @param debugTcrdist3 Boolean controlling whether to run tcrdist3 in debug mode. Default is TRUE.
#' @param verbose Boolean controlling whether to display processing steps. Default is FALSE.
#' @import Matrix
#' @importFrom methods as
#' @importFrom SeuratObject CreateSeuratObject Assays RenameAssays GetAssayData JoinLayers
#' @importFrom Seurat AddMetaData
#' @importFrom methods is
#' @examples
#' \dontrun{
#'   seuratObj_TCR<- RunTcrdist3(inputData = seuratObj@meta.data,
#'               chains = c("TRA", "TRB"),
#'               minimumCloneSize = 2,
#'               calculateChainPairs = FALSE,
#'               rdsOutputPath = "./tcrdist3DistanceMatrices/",
#'               pythonExecutable = reticulate::py_exe(),
#'               verbose = FALSE)
#'
#'   # Example with calculateChainPairs analysis
#'   seuratObj_TCR<- RunTcrdist3(inputData = seuratObj,
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

#TODO: flesh out examples demonstrating formatting requirements for spikeInDataframe


RunTcrdist3 <- function(inputData = NULL,
                        organism = 'human',
                        chains = c("TRA", "TRB"),
                        spikeInDataframe = NULL,
                        minimumCloneSize = 2,
                        calculateChainPairs = FALSE,
                        rdsOutputPath = "./tcrdist3DistanceMatrices/",
                        pythonExecutable = NULL,
                        debugTcrdist3 = TRUE,
                        verbose = FALSE) {

  if (is.null(inputData)) {
    stop("Please provide either a Seurat Object or the Seurat Object's metadata as input.")
  }

  if (typeof(inputData) == 'S4' && class(inputData)[1] == 'Seurat' ){
    metadata <- inputData@meta.data
  }

  if (!is.data.frame(metadata)) {
    stop('Expected inputData to be either a Seurat object or a data.frame')
  }

  if (is.null(pythonExecutable)) {
    pythonExecutable <- reticulate::py_exe()

    #if reticulate doesn't provide a valid executable, try system default
    if (is.null(pythonExecutable) || pythonExecutable == "" || !file.exists(pythonExecutable)) {
      pythonExecutable <- Sys.which("python3")
      if (pythonExecutable == "") {
        pythonExecutable <- Sys.which("python")
      }
    }

    if (pythonExecutable == "" || !file.exists(pythonExecutable)) {
      stop("No valid Python executable found. Please specify pythonExecutable parameter or ensure Python is in PATH.")
    }
  }

  #validate Python environment and required packages
  if (verbose) {
    message(paste("Using Python executable:", pythonExecutable))
  }

  # NOTE: The package is installed as 'tcrdist3' but imported as 'tcrdist'
  #test if required packages are available
  testCmd <- paste0(pythonExecutable, " -c \"import tcrdist; import rpy2; print('Python environment OK')\"")
  testResult <- tryCatch({
    system(testCmd, intern = TRUE)
  }, error = function(e) {
    warning(paste("Python environment test failed:", e$message))
  })

  if (all(is.null(testResult))) {
    stop("Python environment may not have required packages (tcrdist3, rpy2). This may cause failures.")
  }

  rdsOutputPath <- R.utils::getAbsolutePath(rdsOutputPath)
  rdsOutputPath <- gsub(rdsOutputPath, pattern = '\\\\', replacement = '/')

  pythonExecutable <- R.utils::getAbsolutePath(pythonExecutable)

  formatted_metadata <- FormatMetadataForTcrDist3(metadata = metadata,
                                                  chains = chains,
                                                  spikeInDataframe = spikeInDataframe,
                                                  calculateChainPairs = calculateChainPairs,
                                                  verbose = verbose
  )

  if (nrow(formatted_metadata[!formatted_metadata$IsSpikeInClone,]) != nrow(metadata)) {
    stop('Incorrect row count after FormatMetadataForTcrDist3')
  }

  if (any(rownames(formatted_metadata[!formatted_metadata$IsSpikeInClone,]) != rownames(metadata))) {
    print(head(rownames(formatted_metadata[!formatted_metadata$IsSpikeInClone,])))
    print(head(rownames(metadata)))
    stop('Rownames did not match after FormatMetadataForTcrDist3')
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
      chainsString <- paste0(chainsString, "delta")
    } else {
      warning(paste0("Chain Type ", chain, " not recognized. Skipping."))
    }
  }

  # Read the python script template from tcrClustR and write them to a tempfile
  template <- readr::read_file(system.file("scripts/tcrdist3TcrDistances.py", package = "tcrClustR"))
  script <- tempfile(fileext = ".py")
  readr::write_file(template, script)

  # Create distance matrix output directory if it doesn't exist
  if (!dir.exists(rdsOutputPath)) {
    dir.create(rdsOutputPath)
  }

  if (verbose) {
    message("Creating tcrdist3 distance matrices in the following directory: ", rdsOutputPath)
  }

  output_file <- tempfile(fileext = ".csv")
  output_file <- gsub(output_file, pattern = '\\\\', replacement = '/')

  dat <- .filter_and_group_for_tcrdist3(
    metadata = formatted_metadata,
    chains = chain,
    minimumCloneSize = minimumCloneSize
  )

  write.table(dat, file = output_file, sep = ',', row.names = FALSE, quote = TRUE)

  #format and write the python function to the end of the script
  command <- paste0("writeTcrDistances(csv_path = '", output_file,
                    "', organism = '", organism,
                    "', chainsString = '", chainsString,
                    # TODO: review what DB this is using, and probably conditionalize based on organisms
                    "', db_file = 'alphabeta_gammadelta_db.tsv",
                    "', rds_output_path = '", rdsOutputPath,
                    "', debug ='", ifelse(!is.null(debugTcrdist3) && debugTcrdist3, yes = 'True', no = 'False'),
                    "')")
  readr::write_file(command, script, append = TRUE)
  system2(pythonExecutable, script)

  unlink(output_file)

  #return a seurat object, with the distance matrices implemented as assays
  return_values <- list(
    metadata = formatted_metadata
  )

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
      next
    }
    #process the full length TCR distance matrix
    rdsFile <- paste0(rdsOutputPath, "/pw_", chain_tcrdist3, ".rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise 'full' tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_full_length <- readRDS(rdsFile)

    # TODO: GW, can we count on the row/column order being identical to our input??
    colnames(distanceMatrix_full_length) <- as.character(dat$CloneId)
    rownames(distanceMatrix_full_length) <- as.character(dat$CloneId)

    #process the CDR3 only TCR distance matrix
    #grab the first letter of the chain (distance matrices for the cdr3 are stored as "pw_cdr3_b_aa.rds" for a beta chain)
    chain_cdr3_id <- tolower(strsplit(chain_tcrdist3, split = "")[[1]][1])
    rdsFile <- paste0(rdsOutputPath, "/pw_cdr3_", chain_cdr3_id, "_aa.rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise CDR3 tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_CDR3 <- readRDS(rdsFile)
    colnames(distanceMatrix_CDR3) <- as.character(dat$CloneId)
    rownames(distanceMatrix_CDR3) <- as.character(dat$CloneId)

    return_values[[paste0(chain, '_fl')]] <- distanceMatrix_full_length
    rm(distanceMatrix_full_length)

    return_values[[paste0(chain, '_cdr3')]] <- distanceMatrix_CDR3
    rm(distanceMatrix_CDR3)
  }

  for (chain in c(chains)) {
    matName <- paste0(chain, '_cdr3')
    fieldName <- paste0(chain, '-CloneIdx')
    return_values[[paste0(chain, '_cdr3_seurat')]] <- .restore_matrix_to_seurat_obj(formatted_metadata, return_values[[matName]], fieldName)

    matName <- paste0(chain, '_fl')
    return_values[[paste0(chain, '_fl_seurat')]] <- restore_matrix_to_seurat_obj(formatted_metadata, return_values[[matName]], fieldName)
  }

  # TODO: GW, please check my logic on this
  if (calculateChainPairs) {
    message('Calculating joint-chain distances')

    if ('TRA' %in% chains && 'TRB' %in% chains) {
      return_values[[paste0('TRA_TRB_fl_fullSize')]] <-  Seurat::GetAssayData(return_values[[paste0('TRA_fl_fullSize')]], layer = 'counts') + Seurat::GetAssayData(return_values[[paste0('TRB_fl_fullSize')]], layer = 'counts')
      return_values[[paste0('TRA_TRB_cdr3_fl_fullSize')]] <-  Seurat::GetAssayData(return_values[[paste0('TRA_cdr3_fl_fullSize')]], layer = 'counts') + Seurat::GetAssayData(return_values[[paste0('TRB_cdr3_fl_fullSize')]], layer = 'counts')
    }

    if ('TRG' %in% chains && 'TRD' %in% chains) {
      return_values[[paste0('TRG_TRD_fl_fullSize')]] <-  Seurat::GetAssayData(return_values[[paste0('TRG_fl_fullSize')]], layer = 'counts') + Seurat::GetAssayData(return_values[[paste0('TRD_fl_fullSize')]], layer = 'counts')
      return_values[[paste0('TRG_TRD_cdr3_fl_fullSize')]] <-  Seurat::GetAssayData(return_values[[paste0('TRG_cdr3_fl_fullSize')]], layer = 'counts') + Seurat::GetAssayData(return_values[[paste0('TRD_cdr3_fl_fullSize')]], layer = 'counts')
    }
  }

  return(return_values)
}

# NOTE: this will probbaly not work correctly with spike-in values. We could consider subsetting metadata[!metadata$IsSpikeInClone,].
.restore_matrix_to_seurat_obj <- function(metadata, mat, cloneFieldName) {
  if (!cloneFieldName %in% names(metadata)) {
    stop(paste0('Missing field: ', cloneFieldName, ', fields present: ', paste0(sort(names(metadata)), collapse = ',')))
  }

  # NOTE: these can be duplicated:
  cloneNames <- metadata[[cloneFieldName]]

  if (!all(colnames(mat) %in% cloneNames)) {
    stop('The metadata cloneNames and matrix colnames did not match')
  }

  if (!all(rownames(mat) %in% cloneNames)) {
    stop('The metadata cloneNames and matrix rownames did not match')
  }

  # The idea here is to duplicate very clone for each time it appeared in one of the input rows (cells). We can leave the rows/features alone
  # TODO: GW: what's the right way to represent cells that were not compared (such as those lacking a given chain)? NAs might get converted to 0s in a sparse matrix. Zero implies 0-distance, wouldnt it?
  new_mat <- NULL
  for (colName in cloneNames) {
    if (is.na(colName)) {
      new_mat <- cbind(new_mat, tidyr::setNames(rep(nrow(mat), x = NA), colnames(mat)))
    }
    else if (colName %in% colnames(mat)) {
      new_mat <- cbind(new_mat, mat[,colName])
    } else {
      new_mat <- cbind(new_mat, tidyr::setNames(rep(nrow(mat), x = NA), colnames(mat)))
    }
  }
  # Treat these like cells
  colnames(new_mat) <- rownames(metadata)

  return(Seurat::CreateSeuratObject(counts = Seurat::as.Sparse(new_mat), assay = 'TCR', meta.data = metadata))
}

# This is an internal method that expects the dataframe produced by FormatMetadataForTcrDist3. This dataframe should contain columns to uniquely identify each
.filter_and_group_for_tcrdist3 <- function(metadata, chains, minimumCloneSize = 1) {
  initialRows <- nrow(metadata)

  chainId <- paste0(chains, collapse = '_')
  cloneIdxCol <- paste0(chainId, '-CloneIdx')
  if (! cloneIdxCol %in% names(metadata)) {
    stop(paste0('Metadata missing column: ', cloneIdxCol))
  }

  validCol <- paste0(chainId, '_ValidForClustering')
  if (! validCol %in% names(metadata)) {
    stop(paste0('Metadata missing column: ', validCol))
  }

  metadata <- metadata[!is.na(metadata[[cloneIdxCol]]) & metadata[[validCol]],]
  message(paste0('Initial rows: ', initialRows, ', after dropping invalid clones: ', nrow(metadata)))

  # Filter out unique/rare clones
  if (minimumCloneSize > 1) {
    if ('IsSpikeInClone' %in% names(metadata)) {
      metadata <- metadata |>
        dplyr::group_by(dplyr::across(dplyr::all_of(c(cloneIdxCol, 'IsSpikeInClone')))) |>
        dplyr::mutate(count = dplyr::n()) |>
        dplyr::filter(!IsSpikeInClone & count >= minimumCloneSize) |>
        dplyr::ungroup()
    } else {
      metadata <- metadata |>
        dplyr::group_by(dplyr::across(dplyr::all_of(cloneIdxCol))) |>
        dplyr::mutate(count = dplyr::n()) |>
        dplyr::filter(count >= minimumCloneSize) |>
        dplyr::ungroup()
    }

    message(paste0('Rows remaining after filtering clones with cloneSize less than ', minimumCloneSize, ': ', nrow(metadata)))
  }

  if ('IsSpikeInClone' %in% names(metadata)) {
    metadata <- metadata |>
      dplyr::select(-IsSpikeInClone)
  }

  # Construct the grouping columns based on the supplied chains. This deliberately drops TRB if only TRA is selected:
  grouping_columns <- c(cloneIdxCol, 'count')
  for (chain in chains) {
    grouping_columns <- c(grouping_columns, c(chain, paste0(chain, "_V"), paste0(chain, "_J")))
  }

  if (any(! grouping_columns %in% names(metadata))) {
    missing <- grouping_columns[! grouping_columns %in% names(metadata)]
    stop('Missing columns after grouping: ', paste0(missing, collapse = ','))
  }

  metadata <- metadata |>
    dplyr::group_by(dplyr::across(dplyr::all_of(grouping_columns))) |>
    unique.data.frame()

  message(paste0('Unique metadata rows after grouping: ', nrow(metadata)))

  # Reformat data for tcrDist3, iterate over chains specified:
  formatted_data <- data.frame(CloneId = metadata[[cloneIdxCol]], count = metadata$count)

  for (chain in chains){
    charName <- substring(chain, 3)
    if (! paste0('TR', charName) %in% names(metadata)) {
      stop(paste0('missing column: ', paste0('TR', charName)))
    }
    formatted_data[[paste0('v_', tolower(charName), '_gene')]] <- .add_gene_suffix(metadata, paste0('TR', charName, '_V'))
    formatted_data[[paste0('j_', tolower(charName), '_gene')]] <- .add_gene_suffix(metadata, paste0('TR', charName, '_J'))
    formatted_data[[paste0('cdr3_', tolower(charName), '_aa')]] <- metadata[[paste0('TR', charName)]]
    formatted_data[[paste0('cdr3_', tolower(charName), '_nucseq')]] <- sapply(metadata[[paste0('TR', charName)]], .reverse_translate_cdr3)
  }

  return(formatted_data)
}