
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
#' @param spikeInDataframe An optional data frame containing additional CDR3s and gene segments to be included in the clustering, such as published clonotypes.
#' @param minimumCloneSize Minimum number of clones per subject to include in the analysis. Default is 2.
#' @param calculateChainPairs Boolean controlling whether to compute joint multi-chain distance matrices for observed chain combinations.
#' @param pythonExecutable Path to the python executable. Default is reticulate::py_exe().
#' @param debugTcrdist3 Boolean controlling whether to run tcrdist3 in debug mode.
#' @param verbose Boolean controlling whether to display processing steps.
#' @import Matrix
#' @importFrom methods as
#' @importFrom methods is
#' @examples
#' \dontrun{
#'   resultList <- RunTcrdist3(inputData = seuratObj@meta.data,
#'               chains = c("TRA", "TRB"),
#'               minimumCloneSize = 2,
#'               calculateChainPairs = FALSE,
#'               verbose = FALSE)
#'
#'   # Example with calculateChainPairs analysis
#'   resultList <- RunTcrdist3(inputData = seuratObj,
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
                        pythonExecutable = NULL,
                        debugTcrdist3 = FALSE,
                        verbose = FALSE
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

  rdsOutputPath <- tempdir()
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

  seuratObj <- NULL

  for (chain in chains) {
    message(paste0('Preparing chain: ', chain))

    script <- tempfile(fileext = ".py")
    template <- readr::read_file(system.file("scripts/tcrdist3TcrDistances.py", package = "tcrClustR"))
    readr::write_file(template, script)

    input_data_file <- tempfile(fileext = ".csv")
    input_data_file <- gsub(input_data_file, pattern = '\\\\', replacement = '/')
    dat <- .filter_and_group_for_tcrdist3(
      metadata = formatted_metadata,
      chains = chain,
      minimumCloneSize = minimumCloneSize
    )

    write.table(dat, file = input_data_file, sep = ',', row.names = FALSE, quote = TRUE)

    #format and write the python function to the end of the script
    command <- paste0("writeTcrDistances(csv_path = '", input_data_file,
                      "', organism = '", organism,
                      "', chainsString = '", .convert_chain_for_python(chain),
                      # TODO: GW, we should review what DB this is using, and probably conditionalize based on organisms
                      "', db_file = 'alphabeta_gammadelta_db.tsv",
                      "', rds_output_path = '", rdsOutputPath,
                      "', debug ='", ifelse(!is.null(debugTcrdist3) && debugTcrdist3, yes = 'True', no = 'False'),
                      "')")
    readr::write_file(command, script, append = TRUE)
    system2(pythonExecutable, script)

    unlink(input_data_file)
    unlink(script)

    #read the RDS files
    chain_tcrdist3 <- .convert_chain_for_python(chain)
    rdsFile <- paste0(rdsOutputPath, "/pw_", chain_tcrdist3, ".rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise 'full' tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_full_length <- readRDS(rdsFile)

    # TODO: GW, can we count on the row/column order being identical to our input??
    if (nrow(distanceMatrix_full_length) != nrow(dat) || ncol(distanceMatrix_full_length) != nrow(dat)) {
      stop('Expected tcrdist3 fl matrix to have the same dimensions as input')
    }
    colnames(distanceMatrix_full_length) <- as.character(dat$CloneId)
    rownames(distanceMatrix_full_length) <- as.character(dat$CloneId)

    unlink(rdsFile)

    #process the CDR3 only TCR distance matrix
    #grab the first letter of the chain (distance matrices for the cdr3 are stored as "pw_cdr3_b_aa.rds" for a beta chain)
    chain_cdr3_id <- tolower(strsplit(chain_tcrdist3, split = "")[[1]][1])
    rdsFile <- paste0(rdsOutputPath, "/pw_cdr3_", chain_cdr3_id, "_aa.rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise CDR3 tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_CDR3 <- readRDS(rdsFile)
    if (nrow(distanceMatrix_CDR3) != nrow(dat) || ncol(distanceMatrix_CDR3) != nrow(dat)) {
      stop('Expected tcrdist3 fl matrix to have the same dimensions as input')
    }

    colnames(distanceMatrix_CDR3) <- as.character(dat$CloneId)
    rownames(distanceMatrix_CDR3) <- as.character(dat$CloneId)
    unlink(rdsFile)

    fieldName <- paste0(chain, '-CloneIdx')
    seuratObj <- .restore_matrix_to_seurat_obj(formatted_metadata, distanceMatrix_full_length, fieldName, seuratObj = seuratObj, assayName = paste0(chain, '_fl'), assayMeta = dat)
    seuratObj <- .restore_matrix_to_seurat_obj(formatted_metadata, distanceMatrix_CDR3, fieldName, seuratObj = seuratObj, assayName = paste0(chain, '_cdr3'), assayMeta = dat)

    rm(distanceMatrix_CDR3)
    rm(distanceMatrix_full_length)
  }

  # TODO: GW, please check my logic on this
  if (calculateChainPairs) {
    message('Calculating joint-chain distances')

    # The idea below is to subset to clones that are valid in both chains:
    if ('TRA' %in% chains && 'TRB' %in% chains) {
      seuratObj <- .combine_matrices(formatted_metadata, chain1 = 'TRA', chain2 = 'TRB', matSuffix = '_fl', seuratObj = seuratObj)
      seuratObj <- .combine_matrices(formatted_metadata, chain1 = 'TRA', chain2 = 'TRB', matSuffix = '_cdr3', seuratObj = seuratObj)
    }

    if ('TRG' %in% chains && 'TRD' %in% chains) {
      seuratObj <- .combine_matrices(formatted_metadata, chain1 = 'TRG', chain2 = 'TRD', matSuffix = '_fl', seuratObj = seuratObj)
      seuratObj <- .combine_matrices(formatted_metadata, chain1 = 'TRG', chain2 = 'TRD', matSuffix = '_cdr3', seuratObj = seuratObj)
    }
  }

  return(seuratObj)
}

.combine_matrices <- function(metadata, chain1, chain2, matSuffix, seuratObj) {
  print(paste0('Calculating joint distance matrix for: ', chain1, ' and ', chain2))

  assayName <- paste0(chain1, '-', chain2, matSuffix)
  cloneIdxField1 <- paste0(chain1, '-CloneIdx')
  cloneIdxField2 <- paste0(chain2, '-CloneIdx')
  cloneIdxFieldBoth <- paste0(chain1, '_', chain2, '-CloneIdx')
  cloneValidField1 <- paste0(chain1, '_ValidForClustering')
  cloneValidField2 <- paste0(chain2, '_ValidForClustering')

  for (fn in c(cloneIdxField1, cloneIdxField2, cloneIdxFieldBoth, cloneValidField1, cloneValidField2)) {
    if (! fn %in% names(metadata)) {
      stop(paste0('Missing field in metadata: ', fn))
    }
  }

  sel <- !is.na(metadata[[cloneValidField1]]) & metadata[[cloneValidField1]] & !is.na(metadata[[cloneValidField2]]) & metadata[[cloneValidField2]]
  if (sum(sel) == 0) {
    warning('ERROR: No passing dual chains, skipping')
    return(seuratObj)
  }

  print(paste0(assayName, ', total passing cells: ', sum(sel)))
  cloneMapping <- metadata[sel, c(cloneIdxFieldBoth, cloneIdxField1, cloneIdxField2)] %>%
    unique()

  if (any(duplicated(cloneMapping[[cloneIdxFieldBoth]]))) {
    stop('There were duplicated clone IDs')
  }
  rownames(cloneMapping) <- cloneMapping[[cloneIdxFieldBoth]]

  mat1 <- Seurat::GetAssayData(seuratObj, assay = paste0(chain1, matSuffix), layer = 'counts')
  mat2 <- Seurat::GetAssayData(seuratObj, assay = paste0(chain2, matSuffix), layer = 'counts')

  # NOTE: due to the cloneSize filter, this can be a superset of clones in the distance matrix
  validClones1 <- cloneMapping[[cloneIdxField1]]
  validClones1 <- intersect(validClones1, rownames(mat1))
  print(paste0('Total valid chain 1 clones: ', length(validClones1)))

  validClones2 <- cloneMapping[[cloneIdxField2]]
  validClones2 <- intersect(validClones2, rownames(mat2))
  print(paste0('Total valid chain 2 clones: ', length(validClones2)))

  translation1 <- data.frame(x = rownames(mat1)) %>%
    dplyr::inner_join(cloneMapping, by = c('x' = cloneIdxField1))
  rownames(translation1) <- translation1[[cloneIdxFieldBoth]]

  translation2 <- data.frame(x = rownames(mat2)) %>%
    dplyr::inner_join(cloneMapping, by = c('x' = cloneIdxField2))
  rownames(translation2) <- translation2[[cloneIdxFieldBoth]]

  sharedClones <- intersect(translation1[[cloneIdxFieldBoth]], translation2[[cloneIdxFieldBoth]])

  # Now find the intersect across chain1/2:
  translation1 <- translation1[sharedClones,]
  translation2 <- translation2[sharedClones,]

  mat1 <- mat1[translation1$x,]
  rownames(mat1) <- sharedClones

  mat2 <- mat2[translation2$x,]
  rownames(mat2) <- sharedClones

  if (any(colnames(mat1) != colnames(mat2))) {
    stop('Matrix 1 and 2 column names did not match!')
  }

  mat <- Seurat::as.sparse(as.matrix(mat1) + as.matrix(mat2))
  colnames(mat) <- colnames(mat1)
  rownames(mat) <- rownames(mat1)
  seuratObj[[assayName]] <- Seurat::CreateAssayObject(counts = mat)

  # Now make the raw distance matrix:
  mat1 <- seuratObj@misc$TCR_Distances[[paste0(chain1, matSuffix)]]
  mat1 <- mat1[translation1$x, translation1$x]

  mat2 <- seuratObj@misc$TCR_Distances[[paste0(chain2, matSuffix)]]
  mat2 <- mat2[translation2$x, translation2$x]

  mat <- Seurat::as.sparse(as.matrix(mat1) + as.matrix(mat2))
  colnames(mat) <- sharedClones
  rownames(mat) <- sharedClones

  seuratObj@misc$TCR_Distances[[assayName]] <- mat

  return(seuratObj)
}

# NOTE: this will probbaly not work correctly with spike-in values. We could consider subsetting metadata[!metadata$IsSpikeInClone,].
.restore_matrix_to_seurat_obj <- function(metadata, mat, cloneFieldName, seuratObj, assayName, assayMeta) {
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
      new_mat <- cbind(new_mat, stats::setNames(rep(nrow(mat), x = NA), rownames(mat)))
    }
    else if (colName %in% colnames(mat)) {
      new_mat <- cbind(new_mat, mat[,colName])
    } else {
      new_mat <- cbind(new_mat, stats::setNames(rep(nrow(mat), x = NA), rownames(mat)))
    }
  }

  # Treat these like cells
  colnames(new_mat) <- rownames(metadata)

  # These are features/genes, and must be unique:
  if (any(duplicated(rownames(new_mat)))) {
    dups <- unique(rownames(new_mat)[duplicated(rownames(new_mat))])
    stop(paste0('There were duplicated rownames on the final matrix: ', paste0(dups, collapse = ',')))
  }

  if (any(duplicated(colnames(new_mat)))) {
    dups <- unique(colnames(new_mat)[duplicated(colnames(new_mat))])
    stop(paste0('There were duplicated colnames on the final matrix: ', paste0(dups, collapse = ',')))
  }

  if (is.null(rownames(assayMeta))) {
    stop('Missing rownames on assayMeta')
  }

  if (any(rownames(assayMeta) != rownames(new_mat))) {
    stop('Rownames on new_mat do not match assayMeta')
  }

  if (is.null(seuratObj)) {
    seuratObj <- Seurat::CreateSeuratObject(counts = Seurat::as.sparse(new_mat), assay = assayName, meta.data = metadata)
    seuratObj@misc$TCR_Distances <- list()
  } else {
    seuratObj[[assayName]] <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(new_mat))
  }

  seuratObj[[assayName]] <- Seurat::AddMetaData(seuratObj[[assayName]], metadata = assayMeta)
  seuratObj@misc$TCR_Distances[[assayName]] <- mat

  return(seuratObj)
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

  # Always calculate cloneSize
  metadata <- metadata |>
    dplyr::group_by(dplyr::across(dplyr::all_of(cloneIdxCol))) |>
    dplyr::mutate(count = dplyr::n()) |>
    dplyr::ungroup()

  # Filter out unique/rare clones
  if (minimumCloneSize > 1) {
    before_filter <- nrow(metadata)
    if ('IsSpikeInClone' %in% names(metadata)) {
      metadata <- metadata |>
        dplyr::filter(IsSpikeInClone | count >= minimumCloneSize) |>
        dplyr::ungroup()
    } else {
      metadata <- metadata |>
        dplyr::filter(count >= minimumCloneSize) |>
        dplyr::ungroup()
    }

    message(paste0('Rows remaining after filtering clones with cloneSize less than ', minimumCloneSize, ': ', nrow(metadata), ' (total dropped: ', (before_filter - nrow(metadata)), ')'))
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
    dplyr::summarize() |>
    unique.data.frame()

  message(paste0('Unique metadata rows after grouping: ', nrow(metadata)))

  # Reformat data for tcrDist3, iterate over chains specified:
  formatted_data <- data.frame(CloneId = metadata[[cloneIdxCol]], count = metadata$count)
  if (any(duplicated(formatted_data$CloneId))) {
    dupes <- formatted_data$CloneId[duplicated(formatted_data$CloneId)]
    stop(paste0('There were duplicated CloneId values in .filter_and_group_for_tcrdist3(): ', paste0(dupes, collapse = ',')))
  }

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

  rownames(formatted_data) <- formatted_data$CloneId

  return(formatted_data)
}

.convert_chain_for_python <- function(chain) {
  return(switch(chain,
                'TRA' = 'alpha',
                'TRB' = 'beta',
                'TRG' = 'gamma',
                'TRD' = 'delta',
                stop(paste0('Unknown chain: ', chain))
  ))
}
