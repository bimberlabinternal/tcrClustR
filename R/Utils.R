.pythonConfig <- new.env(parent=emptyenv());

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
    chains <- paste0(chains, collapse = '-')
  }

  key <- paste0(chains, '-', ifelse(cdr3Only, yes = 'cdr3', no = 'fl'))
  if (!'TCR_Distances' %in% names(seuratObj_TCR@misc)) {
    stop('This seuratObj does not contain TCR_Distances in the misc slot. Only seurat objects created by ')
  }

  if (!key %in% names(seuratObj_TCR@misc$TCR_Distances)) {
    stop(paste0('Matrix not found, expected: ', key))
  }

  return(seuratObj_TCR@misc$TCR_Distances[[key]])
}

.strip_allele_suffix <- function(){

}