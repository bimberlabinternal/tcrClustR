utils::globalVariables(
  names = c('TRA', 'TRA_V', 'TRA_J', 'TRB', 'TRB_V', 'TRB_J', 'CloneNames', 'count', 'Cluster', 'DistanceSum', 'gene_segments', 'IsSpikeInClone'),
  package = 'tcrClustR',
  add = TRUE
)

#' Format metadata for processing
#' @description This function formats a seurat object's metadata (with TCR information appended) for tcrDist3 distance caluclations.
#'
#' @param metadata Data frame containing metadata.
#' @param chains TCR chains to include in the analysis. TRA/TRB supported and tested, but others likely work.
#' @param organism Organism to use for tcrDist3. Default is 'human'.
#' @param calculateChainPairs If true, this will prepare the columns needed for A/B and/or  G/D (depending on values of chains
#' @param spikeInDataframe Data frame containing spike-in data. Default is NULL. See examples for formatting requirements.
#' @param verbose Boolean controlling whether to display processing steps. Default is FALSE.
#' @return a properly formatted metadata dataframe.
#' @examples
#' \dontrun{
#' spikeInDataframe <- data.frame(CloneNames = rep(1:3),
#'                                  TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
#'                                  TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
#'                                  TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
#'                                  TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
#'                                  TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
#'                                  TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY"))
#' }

.FormatMetadata <- function(metadata,
                                      chains = c("TRA", "TRB"),
                                      organism = 'human',
                                      calculateChainPairs = TRUE,
                                      spikeInDataframe = NULL,
                                      verbose = FALSE
) {
  if (verbose) {
    message("Initial metadata dimensions: ", nrow(metadata), " rows, ", ncol(metadata), " columns")
    message("Requested chains: ", paste(chains, collapse = ", "))
  }

  unknownChains <- chains[!chains %in% c('TRA', 'TRB', 'TRD', 'TRG')]
  if (length(unknownChains) > 0) {
    stop(paste0("Unknown chains: ", paste0(unknownChains, collapse = ',')))
  }

  # Add the corresponding pair, as needed:
  if (calculateChainPairs) {
    if ('TRA' %in% chains || 'TRB' %in% chains) {
      chains <- unique(c(chains, 'TRA', 'TRB'))
    }

    if ('TRG' %in% chains || 'TRD' %in% chains) {
      chains <- unique(c(chains, 'TRG', 'TRD'))
    }
  }

  # Check spikeInDataframe's formatting
  metadata$IsSpikeInClone <- FALSE
  metadata$RowName <- rownames(metadata)
  if (!is.null(spikeInDataframe)) {
    # Check that the spikeInDataframe has columns that match the chains requested
    for (chain in chains) {
      if (! .has_chain_columns(spikeInDataframe, chain)) {
        stop(paste0("The spikeInDataframe must have the columns '", chain, "' (the CDR3), '", paste0(chain, '_V'), "', and '", paste0(chain, '_J'), "'"))
      }
    }

    # Check that the spikeInDataframe has the columns 'CloneNames'
    if (!"CloneNames" %in% colnames(spikeInDataframe)) {
      stop("The spikeInDataframe must have the column 'CloneNames'")
    }

    rownames(spikeInDataframe) <- paste0("SpikeIn_", seq_len(nrow(spikeInDataframe)))

    # Bind the spikeInDataframe to the metadata
    spikeInDataframe$IsSpikeInClone <- TRUE

    newRowNames <- c(rownames(metadata), rownames(spikeInDataframe))
    spikeInDataframe$RowName <- rownames(spikeInDataframe)

    metadata <- plyr::rbind.fill(metadata, spikeInDataframe)
    rownames(metadata) <- newRowNames
    if (any(rownames(metadata) != metadata$RowName)) {
      stop('Row names did not match after plyr::rbind.fill')
    }
  }

  metadata <- .flag_valid_rows(metadata = metadata, chains = chains, organism = organism, verbose = verbose)
  if (nrow(metadata) == 0) {
    stop("No data remaining after filtering. Check your gene segment names and database compatibility.")
  }

  for (chain in chains) {
    colName <- paste0(chain, '_ValidForClustering')
    if (! colName %in% names(metadata)) {
      stop(paste0('Missing column: ', colName))
    }

    metadata$IsValidForChain <- metadata[[colName]]

    tcr_grouping_columns <- c(c(chain, paste0(chain, "_V"), paste0(chain, "_J")))
    metadata <- metadata |>
      tibble::rownames_to_column('RowNames_') |>
      dplyr::group_by(dplyr::across(dplyr::all_of(tcr_grouping_columns))) |>
      dplyr::mutate(`_CloneIdx_` = paste0(chain, '-', 'Clone-', dplyr::cur_group_id()), `_CloneSize_` = dplyr::n()) |>
      dplyr::ungroup() %>%
      tibble::column_to_rownames('RowNames_')

    metadata[['_CloneIdx_']][!metadata$IsValidForChain] <- NA

    # There might be a more elegant solution to this:
    names(metadata)[names(metadata) == '_CloneIdx_'] <- paste0(chain, '_CloneIdx')
    names(metadata)[names(metadata) == '_CloneSize_'] <- paste0(chain, '_CloneSize')
    metadata$IsValidForChain <- NULL
  }

  # Repeat for A/B and G/D:
  if (calculateChainPairs) {
    toTest <- c()
    if ('TRA' %in% chains && 'TRB' %in% chains) {
      toTest <- c(toTest, .get_chain_field_prefix(c('TRA', 'TRB'), delim = '-'))
    }

    if ('TRG' %in% chains && 'TRD' %in% chains) {
      toTest <- c(toTest, .get_chain_field_prefix(c('TRG', 'TRD'), delim = '-'))
    }

    for (chainId in toTest) {
      chains <- unlist(strsplit(chainId, split = '-'))
      col1 <- paste0(chains[1], '_ValidForClustering')
      if (!col1 %in% names(metadata)) {
        stop(paste0('Missing column: ', col1))
      }

      col2 <- paste0(chains[2], '_ValidForClustering')
      if (!col2 %in% names(metadata)) {
        stop(paste0('Missing column: ', col2))
      }

      metadata$IsValidForChain <- metadata[[col1]] & metadata[[col2]]
      tcr_grouping_columns <- c(chains[1], paste0(chains[1], "_V"), paste0(chains[1], "_J"), chains[2], paste0(chains[2], "_V"), paste0(chains[2], "_J"))
      metadata <- metadata |>
        tibble::rownames_to_column('RowNames_') |>
        dplyr::group_by(dplyr::across(dplyr::all_of(tcr_grouping_columns))) |>
        dplyr::mutate(`_CloneIdx_` = paste0(chainId, '-', 'Clone-', dplyr::cur_group_id()), `_CloneSize_` = dplyr::n()) |>
        dplyr::ungroup() %>%
        tibble::column_to_rownames('RowNames_')

      metadata[['_CloneIdx_']][!metadata$IsValidForChain] <- NA

      # NOTE: Seurat doesnt handle hypthens in column names well, so substitute:
      chainIdForColNames <- gsub(chainId, pattern = '-', replacement = '_')
      names(metadata)[names(metadata) == '_CloneIdx_'] <- paste0(chainIdForColNames, '_CloneIdx')
      names(metadata)[names(metadata) == '_CloneSize_'] <- paste0(chainIdForColNames, '_CloneSize')
      metadata$IsValidForChain <- NULL
    }
  }

  if (any(rownames(metadata) != metadata$RowName)) {
    print(head(rownames(metadata)))
    stop('Row names did not match after .FormatMetadata')
  }

  return(metadata)
}

.add_gene_suffix <- function(df, colName, suffix = '*01') {
  if (!colName %in% names(df)) {
    stop('Missing df column: ', colName)
  }

  dat <- df[[colName]]
  sel <- !grepl(dat, pattern = '\\*')
  if (sum(sel) > 0) {
    dat <- as.character(dat)
    dat[sel] <- paste0(dat[sel], suffix)
  }

  return(dat)
}

.reverse_translate_cdr3 <- function(cdr3_aa_seq) {
  # Handle DUMMY or invalid sequences
  if (is.na(cdr3_aa_seq) || length(cdr3_aa_seq) == 0 || grepl("DUMMY", cdr3_aa_seq, fixed = TRUE)) {
    return("AAA")  # Return dummy codon sequence
  }

  #codon table
  codon_table <- list(
    A = c("GCT", "GCC", "GCA", "GCG"),
    C = c("TGT", "TGC"),
    D = c("GAT", "GAC"),
    E = c("GAA", "GAG"),
    `F` = c("TTT", "TTC"),
    G = c("GGT", "GGC", "GGA", "GGG"),
    H = c("CAT", "CAC"),
    I = c("ATT", "ATC", "ATA"),
    K = c("AAA", "AAG"),
    L = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    M = c("ATG"),
    N = c("AAT", "AAC"),
    P = c("CCT", "CCC", "CCA", "CCG"),
    Q = c("CAA", "CAG"),
    R = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    S = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    `T` = c("ACT", "ACC", "ACA", "ACG"),
    V = c("GTT", "GTC", "GTA", "GTG"),
    W = c("TGG"),
    Y = c("TAT", "TAC"),
    "*" = c("TAA", "TAG", "TGA")  #stop codons (unused in this implementation, but could be useful for other applications)
  )

  #sample a reverse translation at random for each AA
  cdr3_nuc_seq <- sapply(strsplit(as.character(cdr3_aa_seq), NULL)[[1]], function(aa) {
    # Check if amino acid exists in codon table
    if (aa %in% names(codon_table) && length(codon_table[[aa]]) > 0) {
      # TODO: GW, do you see problems with this?
      # NOTE: arbitrarily take the first codon when multiple are available, which ensures we can uniquely group rows by CDR3/V/J:
      codon_table[[aa]][1]
    } else {
      # Return AAA for unknown amino acids
      warning("Unknown amino acid '", aa, "' in CDR3 sequence '", cdr3_aa_seq, "'. Using AAA.")
      "AAA"
    }
  })

  #paste the sampled codons together
  return(paste(cdr3_nuc_seq, collapse = ""))
}

.PullTcrdist3Db <- function(organism = 'human',
                            outputFilePath,
                            verbose = FALSE) {

  #validate Python modules are available before proceeding
  pythonExecutable <- .get_python_executable(verbose)
  if (verbose) message("Validating Python dependencies...")

  # NOTE: The package is installed as 'tcrdist3' but imported as 'tcrdist'
  required_modules <- c("tcrdist", "pandas")
  missing_modules <- character(0)

  for (mod in required_modules) {
    check_result <- system2(pythonExecutable,
                            c("-c", shQuote(paste0("import ", mod))),
                            stdout = FALSE,
                            stderr = FALSE)
    if (check_result != 0) {
      missing_modules <- c(missing_modules, mod)
    }
  }

  if (length(missing_modules) > 0) {
    stop("Missing required Python modules: ", paste(missing_modules, collapse = ", "),
         "\n\nPlease run SetupPythonEnvironment() to install dependencies, or manually install with:\n",
         pythonExecutable, " -m pip install ", paste(missing_modules, collapse = " "),
         "\n\nFor tcrdist3, use: ", pythonExecutable, " -m pip install git+https://github.com/bimberlabinternal/tcrdist3.git0.3")
  }

  outputFilePath <- R.utils::getAbsolutePath(outputFilePath)
  template <- readr::read_file(system.file("scripts/PullTcrdist3Db.py", package = "tcrClustR"))
  script <- tempfile(fileext = ".py")
  readr::write_file(template, script)

  command <- paste0("PullTcrdist3Db(organism = '", organism,
                    "', outputFilePath = '", outputFilePath,
                    "')")
  readr::write_file(command, script, append = TRUE)

  # Execute with proper error capture
  if (verbose) {
    message("Python executable: ", pythonExecutable)
    message("Python script: ", script)
    message("Output path: ", outputFilePath)
  }

  # Capture both stdout and stderr
  stdout_file <- tempfile(fileext = ".stdout")
  stderr_file <- tempfile(fileext = ".stderr")

  exit_code <- system2(pythonExecutable,
                       script,
                       stdout = stdout_file,
                       stderr = stderr_file)

  # Read captured output
  stdout_content <- if (file.exists(stdout_file)) readLines(stdout_file, warn = FALSE) else character(0)
  stderr_content <- if (file.exists(stderr_file)) readLines(stderr_file, warn = FALSE) else character(0)

  # Clean up temp files
  unlink(c(stdout_file, stderr_file))

  # Display output if verbose or if there was an error
  if (verbose || exit_code != 0) {
    if (length(stdout_content) > 0) {
      message("Python stdout:")
      cat(paste(stdout_content, collapse = "\n"), "\n")
    }
    if (length(stderr_content) > 0) {
      message("Python stderr:")
      cat(paste(stderr_content, collapse = "\n"), "\n")
    }
  }

  #check that the gene segments file is created
  if (!file.exists(outputFilePath)) {
    error_msg <- paste0(
      "tcrdist3_gene_segments.txt generation failed.\n\n",
      "Python executable: ", pythonExecutable, "\n",
      "Exit code: ", exit_code, "\n\n"
    )

    if (length(stderr_content) > 0) {
      error_msg <- paste0(error_msg, "Python Error Output:\n", paste(stderr_content, collapse = "\n"), "\n\n")
    }

    if (length(stdout_content) > 0) {
      error_msg <- paste0(error_msg, "Python Standard Output:\n", paste(stdout_content, collapse = "\n"), "\n\n")
    }

    error_msg <- paste0(error_msg,
                        "Troubleshooting steps:\n",
                        "1. Run SetupPythonEnvironment() to validate and install Python dependencies\n",
                        "2. Verify tcrdist3 is installed: ", pythonExecutable, " -c 'import tcrdist'\n",
                        "3. Check Python script: ", script, "\n"
    )

    stop(error_msg)
  }

  if (verbose) {
    message("Successfully generated gene segments file: ", outputFilePath)
  }
}

.has_chain_columns <- function(df, chain) {
  vCol <- paste0(chain, '_V')
  jCol <- paste0(chain, '_J')

  return(all(c(chain, vCol, jCol) %in% colnames(df)))
}

.flag_valid_rows <- function(metadata, chains, organism, verbose) {
  segmentTempFile <- tempfile(fileext = '.txt')
  .PullTcrdist3Db(organism = organism,
                  outputFilePath = segmentTempFile,
                  verbose = verbose
  )

  gene_segments_in_db <- readr::read_csv(segmentTempFile, show_col_types = FALSE) |>
    dplyr::mutate(gene_segments = gsub(pattern = "\\*[0-9]+$", replacement = "", gene_segments)) |>
    unlist() |>
    unique()

  gene_segments_in_db <- as.character(gene_segments_in_db)
  if (verbose) {
    message("Found ", length(gene_segments_in_db), " gene segments in tcrdist3 database")
  }

  unlink(segmentTempFile)

  initial_rows <- nrow(metadata)
  if (verbose) {
    message("Starting metadata cleaning with ", nrow(metadata), " rows")
  }

  for (chain in chains) {
    vCol <- paste0(chain, '_V')
    jCol <- paste0(chain, '_J')
    chainColumns <- c(chains, vCol, jCol)
    validCol <- paste0(chain, '_ValidForClustering')

    if (verbose) {
      message("Processing chain: ", chain)
    }

    if (any(!chainColumns %in% colnames(metadata))) {
      missingCols <- chainColumns[! chainColumns %in% colnames(metadata)]
      stop(paste0('Missing column(s): ', paste0(missingCols, collapse = ',')))
    }

    metadata[[validCol]] <- !is.na(metadata[[chain]]) &
      !is.na(metadata[[vCol]]) &
      !is.na(metadata[[jCol]]) &
      metadata[[chain]] != '' &
      metadata[[vCol]] != '' &
      metadata[[jCol]] != ''

    after_na_filter <- sum(metadata[[validCol]])
    if (verbose) {
      message("After filtering NA values in ", chain, ": ", after_na_filter, " rows (filtered ", initial_rows - after_na_filter, " rows)")
    }

    # Filter rows with commas (multiple segments detected in a cell) in the requested chains
    comma_filter <- grepl(",", metadata[[chain]]) | grepl(",", metadata[[vCol]]) | grepl(",", metadata[[jCol]])
    metadata[[validCol]][comma_filter] <- FALSE
    if (verbose) {
      message("Dropping ", sum(comma_filter), ' rows for ', chain, ", remaining valid rows: ", sum(metadata[[validCol]]))
    }

    # Flag gene segments not found in conga's database
    toTest <- as.character(metadata[[vCol]])
    unknown_v_segments <- metadata[[validCol]] & !is.na(toTest) & !(toTest %in% gene_segments_in_db)
    if (sum(unknown_v_segments) > 0) {
      unk <- sort(unique(toTest[unknown_v_segments]))
      warning('The following ', length(unk), ' ', vCol, ' values were not found in the DB: ', paste0(unk, collapse = ','), '. ', paste0("Run tcrClustR:::.PullTcrdist3Db(organism = '", organism, "', outputFilePath = '...') to obtain the list of known segments."))
      metadata[[validCol]][unknown_v_segments] <- FALSE
    }

    toTest <- as.character(metadata[[jCol]])
    unknown_j_segments <- metadata[[validCol]] & !is.na(toTest) & !(toTest %in% gene_segments_in_db)
    if (sum(unknown_j_segments) > 0) {
      unk <- sort(unique(toTest[unknown_j_segments]))
      warning('The following ', length(unk), ' ', jCol, ' values were not found in the DB: ', paste0(unk, collapse = ','), '. ', paste0("Run tcrClustR:::.PullTcrdist3Db(organism = '", organism, "', outputFilePath = '...') to obtain the list of known segments."))
      metadata[[validCol]][unknown_j_segments] <- FALSE
    }

    if (verbose) {
      message(paste0("Finished metadata cleaning for ", chain, ". Initial rows: ", nrow(metadata), ", remaining after filters: ", sum(metadata[[validCol]])))
    }
  }

  return(metadata)
}