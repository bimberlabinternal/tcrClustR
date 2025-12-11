#' @include Utils.R
#' @include FormatMetadata.R
#'
utils::globalVariables(
  names = c('SubjectId', 'TRA_V', 'TRA_J', 'TRB_V', 'TRB_J', 'CloneNames', 'count'),
  package = 'tcrClustR',
  add = TRUE
)

#' @title RunTcrdist3
#' @description This function runs the tcrdist3 pipeline on a set of CDR3 sequences in a Seurat Object.
#'
#' @param seuratObj The seuratObj containing formatted metadata
#' @param organism Organism to use for tcrdist3. Default is 'human'.
#' @param chains Vector of TCR chains to include in the analysis. Default is c("TRA", "TRB").
#' @param minimumCloneSize Minimum number of clones per subject to include in the analysis. Default is 2.
#' @param calculateChainPairs Boolean controlling whether to compute joint multi-chain distance matrices for observed chain combinations.
#' @param debugTcrdist3 Boolean controlling whether to run tcrdist3 in debug mode.
#' @param verbose Boolean controlling whether to display processing steps.
#' @import Matrix
#' @importFrom methods as
#' @importFrom methods is
RunTcrdist3 <- function(seuratObj,
                        organism = 'human',
                        chains = c("TRA", "TRB"),
                        minimumCloneSize = 2,
                        calculateChainPairs = FALSE,
                        debugTcrdist3 = FALSE,
                        verbose = FALSE
) {
  pythonExecutable <- .get_python_executable(verbose)

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

  for (chain in chains) {
    message(paste0('Preparing chain: ', chain))

    script <- tempfile(fileext = ".py")
    template <- readr::read_file(system.file("scripts/tcrdist3TcrDistances.py", package = "tcrClustR"))
    readr::write_file(template, script)

    input_data_file <- tempfile(fileext = ".csv")
    input_data_file <- gsub(input_data_file, pattern = '\\\\', replacement = '/')
    dat <- .filter_and_group_for_tcrdist3(
      metadata = seuratObj@meta.data,
      chains = chain,
      minimumCloneSize = minimumCloneSize
    )

    utils::write.table(dat, file = input_data_file, sep = ',', row.names = FALSE, quote = TRUE)

    #format and write the python function to the end of the script
    rdsOutputPath <- gsub(tempdir(), pattern = '\\\\', replacement = '/')
    command <- paste0("writeTcrDistances(csv_path = '", input_data_file,
                      "', organism = '", organism,
                      "', chainsString = '", .convert_chain_for_python(chain),
                      "', db_file = 'combo_xcr_2024-03-05.tsv",
                      "', rds_output_path = '", rdsOutputPath,
                      "', debug ='", ifelse(!is.null(debugTcrdist3) && debugTcrdist3, yes = 'True', no = 'False'),
                      "')")
    readr::write_file(command, script, append = TRUE)
    pythonOutput <- system2(pythonExecutable, script, stdout = TRUE, stderr = TRUE)
    if (length(pythonOutput) > 0) {
      message(paste(pythonOutput, collapse = "\n"))
    }
    
    #ensure cell_df vs clone_df order validation occured in the python layer
    chain_tcrdist3 <- .convert_chain_for_python(chain)
    chain_cdr3_id <- tolower(strsplit(chain_tcrdist3, split = "")[[1]][1])

    unlink(input_data_file)
    unlink(script)

    #read the RDS files
    
    #read clone_ids from Python (dimnames set in R since rpy2 dimnames assignment fails)
    cloneIdsFile <- paste0(rdsOutputPath, "/clone_ids.rds")
    if (!file.exists(cloneIdsFile)) {
      stop(paste0("Clone IDs RDS file not found: ", cloneIdsFile))
    }
    python_clone_ids <- readRDS(cloneIdsFile)
    
    #validate clone_ids match input data
    if (length(python_clone_ids) != nrow(dat)) {
      stop(paste0('Clone IDs length (', length(python_clone_ids), ') does not match input data rows (', nrow(dat), ')'))
    }
    if (!identical(as.character(python_clone_ids), as.character(dat$CloneId))) {
      stop('Clone IDs from Python do not match input CloneId values (order or content mismatch)')
    }
    
    rdsFile <- paste0(rdsOutputPath, "/pw_", chain_tcrdist3, ".rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise 'full' tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_full_length <- readRDS(rdsFile)

    if (nrow(distanceMatrix_full_length) != length(python_clone_ids) || ncol(distanceMatrix_full_length) != length(python_clone_ids)) {
      stop('Expected tcrdist3 fl matrix to have the same dimensions as clone_ids')
    }
    colnames(distanceMatrix_full_length) <- as.character(python_clone_ids)
    rownames(distanceMatrix_full_length) <- as.character(python_clone_ids)

    unlink(rdsFile)

    #process the CDR3 only TCR distance matrix
    #grab the first letter of the chain (distance matrices for the cdr3 are stored as "pw_cdr3_b_aa.rds" for a beta chain)
    rdsFile <- paste0(rdsOutputPath, "/pw_cdr3_", chain_cdr3_id, "_aa.rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise CDR3 tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_CDR3 <- readRDS(rdsFile)
    if (nrow(distanceMatrix_CDR3) != length(python_clone_ids) || ncol(distanceMatrix_CDR3) != length(python_clone_ids)) {
      stop('Expected tcrdist3 cdr3 matrix to have the same dimensions as clone_ids')
    }

    colnames(distanceMatrix_CDR3) <- as.character(python_clone_ids)
    rownames(distanceMatrix_CDR3) <- as.character(python_clone_ids)
    unlink(rdsFile)
    unlink(cloneIdsFile)

    fieldName <- paste0(chain, '_CloneIdx')
    seuratObj <- .validate_distance_mat_and_store(seuratObj, distanceMatrix_full_length, fieldName, assayName = paste0(chain, '_fl'), assayMeta = dat)
    seuratObj <- .validate_distance_mat_and_store(seuratObj, distanceMatrix_CDR3, fieldName, assayName = paste0(chain, '_cdr3'), assayMeta = dat)

    rm(distanceMatrix_CDR3)
    rm(distanceMatrix_full_length)
  }

  if (calculateChainPairs) {
    message('Calculating joint-chain distances')

    # The idea below is to subset to clones that are valid in both chains:
    if ('TRA' %in% chains && 'TRB' %in% chains) {
      seuratObj <- .combine_matrices(seuratObj, chain1 = 'TRA', chain2 = 'TRB', matSuffix = '_fl')
      seuratObj <- .combine_matrices(seuratObj, chain1 = 'TRA', chain2 = 'TRB', matSuffix = '_cdr3')
    }

    if ('TRG' %in% chains && 'TRD' %in% chains) {
      seuratObj <- .combine_matrices(seuratObj, chain1 = 'TRG', chain2 = 'TRD', matSuffix = '_fl')
      seuratObj <- .combine_matrices(seuratObj, chain1 = 'TRG', chain2 = 'TRD', matSuffix = '_cdr3')
    }
  }

  return(seuratObj)
}

.combine_matrices <- function(seuratObj, chain1, chain2, matSuffix) {
  message(paste0('Calculating joint distance matrix for: ', chain1, ' and ', chain2))

  assayName <- paste0(chain1, ASSAY_MULTI_CHAIN_DELIM, chain2, matSuffix)
  cloneIdxField1 <- paste0(chain1, '_CloneIdx')
  cloneIdxField2 <- paste0(chain2, '_CloneIdx')
  cloneIdxFieldBoth <- paste0(chain1, ASSAY_MULTI_CHAIN_DELIM, chain2, '_CloneIdx')
  cloneValidField1 <- paste0(chain1, '_ValidForClustering')
  cloneValidField2 <- paste0(chain2, '_ValidForClustering')

  for (fn in c(cloneIdxField1, cloneIdxField2, cloneIdxFieldBoth, cloneValidField1, cloneValidField2)) {
    if (! fn %in% names(seuratObj@meta.data)) {
      stop(paste0('Missing field in metadata: ', fn))
    }
  }

  sel <- !is.na(seuratObj@meta.data[[cloneValidField1]]) & seuratObj@meta.data[[cloneValidField1]] & !is.na(seuratObj@meta.data[[cloneValidField2]]) & seuratObj@meta.data[[cloneValidField2]]
  if (sum(sel) == 0) {
    warning('ERROR: No passing dual chains, skipping')
    return(seuratObj)
  }

  message(paste0(assayName, ', total passing cells: ', sum(sel)))
  cloneMapping <- seuratObj@meta.data[sel, c(cloneIdxFieldBoth, cloneIdxField1, cloneIdxField2)] %>%
    unique()

  if (any(duplicated(cloneMapping[[cloneIdxFieldBoth]]))) {
    stop('There were duplicated clone IDs')
  }
  rownames(cloneMapping) <- cloneMapping[[cloneIdxFieldBoth]]

  mat1 <- Seurat::GetAssayData(seuratObj@misc$TCR_Distances[[paste0(chain1, matSuffix)]], layer = 'counts')
  mat2 <- Seurat::GetAssayData(seuratObj@misc$TCR_Distances[[paste0(chain2, matSuffix)]], layer = 'counts')

  # NOTE: due to the cloneSize filter, this can be a superset of clones in the distance matrix
  validClones1 <- cloneMapping[[cloneIdxField1]]
  validClones1 <- intersect(validClones1, rownames(mat1))
  message(paste0('Total valid chain 1 clones: ', length(validClones1)))

  validClones2 <- cloneMapping[[cloneIdxField2]]
  validClones2 <- intersect(validClones2, rownames(mat2))
  message(paste0('Total valid chain 2 clones: ', length(validClones2)))

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

  # Now make the raw distance matrix:
  mat1 <- mat1[translation1$x, translation1$x]
  mat2 <- mat2[translation2$x, translation2$x]

  mat <- Seurat::as.sparse(as.matrix(mat1) + as.matrix(mat2))
  colnames(mat) <- sharedClones
  rownames(mat) <- sharedClones

  seuratObj@misc$TCR_Distances[[assayName]] <- Seurat::CreateAssayObject(counts = mat)

  return(seuratObj)
}

.validate_distance_mat_and_store <- function(seuratObj, distanceMat, cloneFieldName, assayName, assayMeta) {
  if (!cloneFieldName %in% names(seuratObj@meta.data)) {
    stop(paste0('Missing field: ', cloneFieldName, ', fields present: ', paste0(sort(names(seuratObj@meta.data)), collapse = ',')))
  }

  # NOTE: these can be duplicated:
  cloneNames <- seuratObj@meta.data[[cloneFieldName]]

  if (!all(colnames(distanceMat) %in% cloneNames)) {
    stop('The metadata cloneNames and matrix colnames did not match')
  }

  if (!all(rownames(distanceMat) %in% cloneNames)) {
    stop('The metadata cloneNames and matrix rownames did not match')
  }

  seuratObj@misc$TCR_Distances[[assayName]] <- Seurat::CreateAssayObject(counts = Seurat::as.sparse(distanceMat))
  seuratObj@misc$TCR_Distances[[assayName]] <- Seurat::AddMetaData(seuratObj@misc$TCR_Distances[[assayName]], metadata = assayMeta)

  return(seuratObj)
}

# This is an internal method that expects the dataframe produced by .FormatMetadata. This dataframe should contain columns to uniquely identify each
.filter_and_group_for_tcrdist3 <- function(metadata, chains, minimumCloneSize = 1) {
  initialRows <- nrow(metadata)

  chainId <- .get_chain_field_prefix(chains)
  cloneIdxCol <- paste0(chainId, '_CloneIdx')
  if (! cloneIdxCol %in% names(metadata)) {
    print(names(metadata))
    stop(paste0('Metadata missing column: ', cloneIdxCol))
  }

  validCol <- paste0(chainId, '_ValidForClustering')
  if (! validCol %in% names(metadata)) {
    print(names(metadata))
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

  # Drop the column if unused:
  if ('IsSpikeInClone' %in% names(metadata) && length(unique(metadata$IsSpikeInClone)) == 1) {
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
