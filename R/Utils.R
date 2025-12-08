.pythonConfig <- new.env(parent=emptyenv());

ASSAY_MULTI_CHAIN_DELIM <- '_'

#' @title SetPythonExecutable
#' @description Can be used to set the python executable used by tcrClustR
#'
#' @param pythonExecutable The path to the python executable
#' @export
SetPythonExecutable <- function(pythonExecutable = ""){
  if (pythonExecutable != "") {
    if (!!file.exists(pythonExecutable)) {
      stop(paste0('File does not exist: ', pythonExecutable))
    }

    pythonExecutable <- R.utils::getAbsolutePath(pythonExecutable)
    .pythonConfig$pythonExecutable <- pythonExecutable
  }
}


.get_python_executable <- function(verbose = FALSE) {
  if (!is.null(.pythonConfig$pythonExecutable) && .pythonConfig$pythonExecutable != '') {
    if (verbose) {
      message(paste0('Using user-supplied python executable: ', .pythonConfig$pythonExecutable))
    }

    return(.pythonConfig$pythonExecutable)
  }

  pythonExecutable <- reticulate::py_exe()
  if (is.null(pythonExecutable) || pythonExecutable == '') {
    pythonExecutable <- Sys.which("python3")
  }

  if (is.null(pythonExecutable) || pythonExecutable == "") {
    stop("No valid Python executable found. Please specify pythonExecutable parameter or ensure Python is in PATH.")
  }

  pythonExecutable <- R.utils::getAbsolutePath(pythonExecutable)
  if (!file.exists(pythonExecutable)) {
    stop(paste0('Python executable does not exist: ', pythonExecutable))
  }

  #validate Python environment and required packages
  if (verbose) {
    message(paste("Using Python executable:", pythonExecutable))
  }

  return(pythonExecutable)
}

#' @title GetDistanceMatrix
#' @description Can be used to set the python executable used by tcrClustR
#'
#' @param seuratObj_TCR The seurat object, created by CalculateTcrDistances
#' @param chains The chain or chains to return (e.g., TRA or c('TRG', 'TRD'))
#' @param cdr3Only If true, it will return the distances calculated on CDR3 alone. The default uses distances based on CDR3+V/J
#'
#' @export
GetDistanceMatrix <- function(seuratObj_TCR, chains, cdr3Only = FALSE) {
  if (length(chains) > 1) {
    chains <- .get_chain_field_prefix(chains)
  }

  key <- paste0(chains, '_', ifelse(cdr3Only, yes = 'cdr3', no = 'fl'))
  if (!'TCR_Distances' %in% names(seuratObj_TCR@misc)) {
    stop('This seuratObj does not contain TCR_Distances in the misc slot. Only seurat objects created by ')
  }

  if (!key %in% names(seuratObj_TCR@misc$TCR_Distances)) {
    stop(paste0('Matrix not found, expected: ', key))
  }

  return(Seurat::GetAssayData(seuratObj_TCR@misc$TCR_Distances[[key]], layer = 'counts'))
}

# The primary purpose of this function is to ensure consistent order for TRA/TRB and TRG/TRD on joint assays:
.get_chain_field_prefix <- function(chains, delim = '_') {
  if (length(chains) == 1){
    return(chains)
  }

  if ('TRA' %in% chains && 'TRB' %in% chains) {
    return(paste0('TRA', delim, 'TRB'))
  } else if ('TRG' %in% chains && 'TRD' %in% chains) {
    return(paste0('TRG', delim, 'TRD'))
  } else {
    stop(paste0('Unexpected chain combination: ', paste0(chains, collapse = ',')))
  }
}

#' @title ExpandDistancesToMatchSeuratObj
#' @description This takes the distance matrix, which is one row/col per clone, and joined this with the seuratObj in order to make a larger matrix with one column/cell
#'
#' @param seuratObj The seurat object
#' @param chains A single chain or vector of the joint chains to use
#' @export
ExpandDistancesToMatchSeuratObj <- function(seuratObj, chains) {
  if (length(chains) > 1) {
    chains <- .get_chain_field_prefix(chains)
  }
  cloneFieldName <- paste0(chains, '-ClonexIdx')
  if (!cloneFieldName %in% names(seuratObj@meta.data)) {
    stop(paste0('Missing field: ', cloneFieldName, ', fields present: ', paste0(sort(names(seuratObj@meta.data)), collapse = ',')))
  }

  mat <- GetDistanceMatrix(seuratObj, chains)

  # NOTE: these can be duplicated:
  cloneNames <- seuratObj@meta.data[[cloneFieldName]]

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
      new_mat <- cbind(new_mat, stats::setNames(rep(nrow(mat), x = NA), rownames(mat)))
    }
    else if (colName %in% colnames(mat)) {
      new_mat <- cbind(new_mat, mat[,colName])
    } else {
      new_mat <- cbind(new_mat, stats::setNames(rep(nrow(mat), x = NA), rownames(mat)))
    }
  }

  # Treat these like cells
  colnames(new_mat) <- rownames(seuratObj@meta.data)

  # These are features/genes, and must be unique:
  if (any(duplicated(rownames(new_mat)))) {
    dups <- unique(rownames(new_mat)[duplicated(rownames(new_mat))])
    stop(paste0('There were duplicated rownames on the final matrix: ', paste0(dups, collapse = ',')))
  }

  if (any(duplicated(colnames(new_mat)))) {
    dups <- unique(colnames(new_mat)[duplicated(colnames(new_mat))])
    stop(paste0('There were duplicated colnames on the final matrix: ', paste0(dups, collapse = ',')))
  }

  return(Seurat::as.sparse(new_mat))
}