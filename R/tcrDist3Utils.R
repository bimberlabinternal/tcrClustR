
utils::globalVariables(
  names = c('SubjectId', 'TRA', 'TRA_V', 'TRA_J', 'TRB', 'TRB_V', 'TRB_J', 'CloneNames', 'count'),
  package = 'tcrClustR',
  add = TRUE
)

#' Format metadata for tcrDist3
#' @description This function formats a seurat object's metadata (with TCR information appended) for tcrDist3 distance caluclations.
#'
#' @param metadata Data frame containing metadata.
#' @param chains TCR chains to include in the analysis. TRA/TRB supported and tested, but others likely work.
#' @param organism Organism to use for tcrDist3. Default is 'human'.
#' @param cleanMetadata Boolean controlling whether to clean the metadata by removing rows with NA values or commas in the specified chains.
#' @param summarizeClones Boolean controlling whether to summarize clones by SubjectId, TRA, TRB, TRA_V, TRA_J, TRB_V, and TRB_J.
#' @param imputeCloneNames Boolean controlling whether to impute clone names if they are missing.
#' @param writeUnannotatedGeneSegmentsToFile Boolean controlling whether to write unannotated gene segments to a file (filtered_(chain)_gene_segments.csv).
#' @param outputCsv Path to the output CSV file.
#' @param minimumClonesPerSubject Minimum number of clones per subject to include in the analysis. Default is 2.
#' @param spikeInDataframe Data frame containing spike-in data. Default is NULL. See examples for formatting requirements.
#' @param pythonExecutable Path to the python executable. Default is NULL, but imputes to reticulate::py_exe().
#' @param verbose Boolean controlling whether to display processing steps. Default is FALSE.
#' @return a properly formatted metadata dataframe.
#' @export
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

#TODO: flesh out examples demonstrating formatting requirements for spikeInDataframe

FormatMetadataForTcrDist3 <- function(metadata,
                                      outputCsv = './tcrDist3Input.csv',
                                      chains = c("TRA", "TRB"),
                                      organism = 'human',
                                      cleanMetadata = T,
                                      summarizeClones = T,
                                      imputeCloneNames = T,
                                      minimumClonesPerSubject = 100,
                                      writeUnannotatedGeneSegmentsToFile = T,
                                      spikeInDataframe = NULL,
                                      pythonExecutable = NULL,
                                      verbose = FALSE
) {
  if (verbose) {
    message("Initial metadata dimensions: ", nrow(metadata), " rows, ", ncol(metadata), " columns")
    message("Requested chains: ", paste(chains, collapse = ", "))
  }
  #check spikeInDataframe's formatting
  if (!is.null(spikeInDataframe)) {
    #check that the spikeInDataframe has columns that match the chains requested
    if ("TRA" %in% chains) {
      if (!all(c("TRA_V", "TRA_J", "TRA") %in% colnames(spikeInDataframe))) {
        stop("The spikeInDataframe must have the columns 'TRA_V', 'TRA_J', and 'CDR3' for TRA chains.")
      }
    } else if ("TRB" %in% chains) {
      if (!all(c("TRB_V", "TRB_J", "TRB") %in% colnames(spikeInDataframe))) {
        stop("The spikeInDataframe must have the columns 'TRB_V', 'TRB_J', and 'CDR3' for TRB chains.")
      }
    } else if ("TRG" %in% chains) {
      if (!all(c("TRG_V", "TRG_J", "TRG") %in% colnames(spikeInDataframe))) {
        stop("The spikeInDataframe must have the columns 'TRG_V', 'TRG_J', and 'CDR3' for TRG chains.")
      }
    } else if ("TRD" %in% chains) {
      if (!all(c("TRD_V", "TRD_J", "TRD") %in% colnames(spikeInDataframe))) {
        stop("The spikeInDataframe must have the columns 'TRD_V', 'TRD_J', and 'CDR3' for TRD chains.")
      }
    } else {
      stop(paste0("Chain ", chain, " is not supported."))
    }
    #check that the spikeInDataframe has the columns 'CloneNames' and impute a SubjectId if missing
    if (!"CloneNames" %in% colnames(spikeInDataframe)) {
      stop("The spikeInDataframe must have the column 'CloneNames'.")
    }
    if (!"SubjectId" %in% colnames(spikeInDataframe)) {
      spikeInDataframe$SubjectId <- paste0("spikeIn_", seq_len(nrow(spikeInDataframe)))
    }
    #force spikeInDataframe to exceed the minimum number of clones per subject
    if (minimumClonesPerSubject > 1) {
      spikeInDataframe <- do.call("rbind",
                                  replicate(minimumClonesPerSubject,
                                            spikeInDataframe,
                                            simplify = FALSE))
    }
    #bind the spikeInDataframe to the metadata
    metadata <- plyr::rbind.fill(metadata, spikeInDataframe)
  }
  if (cleanMetadata) {
    if (verbose) message("Starting metadata cleaning with ", nrow(metadata), " rows")
    for (chain in chains) {
      if (verbose) message("Processing chain: ", chain)

      #check if chain column exists
      if (!chain %in% colnames(metadata)) {
        if (verbose) {
          message("WARNING - Chain column ", chain, " not found in metadata!")
          message("Available columns: ", paste(colnames(metadata), collapse = ", "))
        }
        next
      }

      initial_rows <- nrow(metadata)

      #filter rows with NA values in the requested chains
      metadata <- metadata[!is.na(metadata[[chain]]), ]
      after_na_filter <- nrow(metadata)
      if (verbose) message("After filtering NA values in ", chain, ": ", after_na_filter, " rows (removed ", initial_rows - after_na_filter, " rows)")

      #filter rows with commas (multiple segments detected in a cell) in the requested chains
      metadata <- metadata[!grepl(",", metadata[[chain]]), ]
      after_comma_filter <- nrow(metadata)
      if (verbose) message("After filtering commas in ", chain, ": ", after_comma_filter, " rows (removed ", after_na_filter - after_comma_filter, " rows)")

      .PullTcrdist3Db(organism = organism,
                      outputFilePath = file.path(dirname(outputCsv), 'tcrdist3_gene_segments.txt'),
                      verbose = verbose)
      gene_segments_in_db <- readr::read_csv(file.path(dirname(outputCsv), 'tcrdist3_gene_segments.txt'), show_col_types = FALSE) |>
        dplyr::mutate(`gene_segments` = gsub("\\*[0-9]+$", "", `gene_segments`)) |>
        unlist() |>
        unique()
      if (verbose) {
        message("Found ", length(gene_segments_in_db), " gene segments in tcrdist3 database")
        message("First 10 database gene segments: ", paste(head(gene_segments_in_db, 10), collapse = ", "))
      }

      #remove gene segments not found in conga's database
      #TODO: add message flag that trips when TCRDist3 detects an unannotated gene segment
      if (chain == "TRA") {
        if (writeUnannotatedGeneSegmentsToFile) {
          if (any(!(metadata$TRA_V %in% gene_segments_in_db) | any(!metadata$TRA_J %in% gene_segments_in_db))) {
            message("TRA gene segments present in the data, but not found in conga database!")
            #store filtered gene segments
            filtered_genes <- metadata |>
              dplyr::filter(!is.na(TRA_V)) |>
              dplyr::filter(!is.na(TRA_J)) |>
              dplyr::mutate(
                VALID_V = dplyr::case_when(TRA_V %in% gene_segments_in_db ~ "valid",
                                           TRUE ~ "invalid"),
                VALID_J = dplyr::case_when(TRA_J %in% gene_segments_in_db ~ "valid",
                                           TRUE ~ "invalid" )) |>
              #write only the invalid V and J segments, so that the valid ones appear as "valid" in the text file
              dplyr::mutate(TRA_V = dplyr::case_when(VALID_V == "valid" ~ "valid",
                                                     TRUE ~ TRA_V)) |>
              dplyr::mutate(TRA_J = dplyr::case_when(VALID_J == "valid" ~ "valid",
                                                     TRUE ~ TRA_J)) |>
              dplyr::filter(TRA_V != "valid" | TRA_J != "valid") |>
              dplyr::select(TRA_V, TRA_J) |>
              unique.data.frame()
            if (verbose) message("Writing TRA segments present in the data, but missing in tcrdist3 database to file: ", R.utils::getAbsolutePath(file.path(dirname(outputCsv),'filtered_TRA_gene_segments.csv')))
            utils::write.csv(filtered_genes, file = file.path(dirname(outputCsv),'filtered_TRA_gene_segments.csv'), row.names = FALSE)
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRA_V %in% gene_segments_in_db) |>
          dplyr::filter(TRA_J %in% gene_segments_in_db)
      } else if (chain == "TRB") {
        if (verbose) message("Processing TRB chain filtering")

        #check for required TRB columns
        if (!"TRB_V" %in% colnames(metadata)) {
          if (verbose) {
            message("ERROR - TRB_V column not found in metadata!")
            message("Available columns: ", paste(colnames(metadata), collapse = ", "))
          }
        }
        if (!"TRB_J" %in% colnames(metadata)) {
          if (verbose) {
            message("ERROR - TRB_J column not found in metadata!")
            message("Available columns: ", paste(colnames(metadata), collapse = ", "))
          }
        }

        #check unique values in TRB_V and TRB_J columns
        if ("TRB_V" %in% colnames(metadata)) {
          unique_v_genes <- unique(metadata$TRB_V[!is.na(metadata$TRB_V)])
          if (verbose) {
            message("Found ", length(unique_v_genes), " unique TRB_V genes")
            message("First 10 TRB_V genes: ", paste(head(unique_v_genes, 10), collapse = ", "))
          }
          v_genes_in_db <- sum(unique_v_genes %in% gene_segments_in_db)
          if (verbose) message(v_genes_in_db, " out of ", length(unique_v_genes), " TRB_V genes found in database")
        }

        if ("TRB_J" %in% colnames(metadata)) {
          unique_j_genes <- unique(metadata$TRB_J[!is.na(metadata$TRB_J)])
          if (verbose) {
            message("Found ", length(unique_j_genes), " unique TRB_J genes")
            message("First 10 TRB_J genes: ", paste(head(unique_j_genes, 10), collapse = ", "))
          }
          j_genes_in_db <- sum(unique_j_genes %in% gene_segments_in_db)
          if (verbose) message(j_genes_in_db, " out of ", length(unique_j_genes), " TRB_J genes found in database")
        }

        before_v_filter <- nrow(metadata)
        if (writeUnannotatedGeneSegmentsToFile) {
          if (any(!(metadata$TRB_V %in% gene_segments_in_db) | any(!metadata$TRB_J %in% gene_segments_in_db))) {
            message("TRB gene segments present in the data, but not found in conga database!")
            #store filtered gene segments
            filtered_genes <- metadata |>
              dplyr::filter(!is.na(TRB_V)) |>
              dplyr::filter(!is.na(TRB_J)) |>
              dplyr::mutate(
                VALID_V = dplyr::case_when(TRB_V %in% gene_segments_in_db ~ "valid",
                                           TRUE ~ "invalid"),
                VALID_J = dplyr::case_when(TRB_J %in% gene_segments_in_db ~ "valid",
                                           TRUE ~ "invalid" )) |>
              #NA the valid V and J segments, so that the valid ones appear as "valid" in the text file
              dplyr::mutate(TRB_V = dplyr::case_when(VALID_V == "valid" ~ "valid",
                                                     TRUE ~ TRB_V)) |>
              dplyr::mutate(TRB_J = dplyr::case_when(VALID_J == "valid" ~ "valid",
                                                     TRUE ~ TRB_J)) |>
              dplyr::filter(TRB_V != "valid" | TRB_J != "valid") |>
              dplyr::select(TRB_V, TRB_J) |>
              unique.data.frame()
            if (verbose) message("Writing TRB segments present in the data, but missing in tcrdist3 database to file: ", R.utils::getAbsolutePath(file.path(dirname(outputCsv),'filtered_TRB_gene_segments.csv')))
            utils::write.csv(filtered_genes, file = file.path(dirname(outputCsv), 'filtered_TRB_gene_segments.csv'), row.names = FALSE)
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRB_V %in% gene_segments_in_db) |>
          dplyr::filter(TRB_J %in% gene_segments_in_db)
        after_gene_filter <- nrow(metadata)
        if (verbose) {
          message("After filtering TRB gene segments: ", after_gene_filter, " rows (removed ", before_v_filter - after_gene_filter, " rows)")
          message("Metadata after TRB filtering has ", nrow(metadata), " rows")
        }

      } else if (chain == "TRG") {
        if (writeUnannotatedGeneSegmentsToFile) {
          if (any(!(metadata$TRG_V %in% gene_segments_in_db) | any(!metadata$TRG_J %in% gene_segments_in_db))) {
            message("TRG gene segments present in the data, but not found in conga database!")
            #store filtered gene segments
            filtered_genes <- metadata |>
              dplyr::filter(!is.na(TRG_V)) |>
              dplyr::filter(!is.na(TRG_J)) |>
              dplyr::mutate(
                VALID_V = dplyr::case_when(TRG_V %in% gene_segments_in_db ~ "valid",
                                           TRUE ~ "invalid"),
                VALID_J = dplyr::case_when(TRG_J %in% gene_segments_in_db ~ "valid",
                                           TRUE ~ "invalid" )) |>
              #NA the valid V and J segments, so that the valid ones appear as "valid" in the text file
              dplyr::mutate(TRG_V = dplyr::case_when(VALID_V == "valid" ~ "valid",
                                                     TRUE ~ TRG_V)) |>
              dplyr::mutate(TRG_J = dplyr::case_when(VALID_J == "valid" ~ "valid",
                                                     TRUE ~ TRG_J)) |>
              dplyr::filter(TRG_V != "valid" | TRG_J != "valid") |>
              dplyr::select(TRG_V, TRG_J) |>
              unique.data.frame()
            if (verbose) message("Writing TRG segments present in the data, but missing in tcrdist3 database to file: ", R.utils::getAbsolutePath(file.path(dirname(outputCsv),'filtered_TRG_gene_segments.csv')))
            utils::write.csv(filtered_genes, file = file.path(dirname(outputCsv),'filtered_TRG_gene_segments.csv'), row.names = FALSE)
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRG_V %in% gene_segments_in_db) |>
          dplyr::filter(TRG_J %in% gene_segments_in_db)
      } else if (chain == "TRD") {
        if (writeUnannotatedGeneSegmentsToFile) {
          if (any(!(metadata$TRD_V %in% gene_segments_in_db) | any(!metadata$TRD_J %in% gene_segments_in_db))) {
            message("TRD gene segments present in the data, but not found in conga database!")
            #store filtered gene segments
            filtered_genes <- metadata |>
              dplyr::filter(!is.na(TRD_V)) |>
              dplyr::filter(!is.na(TRD_J)) |>
              dplyr::mutate(
                VALID_V = dplyr::case_when(TRD_V %in% gene_segments_in_db ~ "valid",
                                           TRUE ~ "invalid"),
                VALID_J = dplyr::case_when(TRD_J %in% gene_segments_in_db ~ "valid",
                                           TRUE ~ "invalid" )) |>
              #NA the valid V and J segments, so that the valid ones appear as "valid" in the text file
              dplyr::mutate(TRD_V = dplyr::case_when(VALID_V == "valid" ~ "valid",
                                                     TRUE ~ TRD_V)) |>
              dplyr::mutate(TRD_J = dplyr::case_when(VALID_J == "valid" ~ "valid",
                                                     TRUE ~ TRD_J)) |>
              dplyr::filter(TRD_V != "valid" | TRD_J != "valid") |>
              dplyr::select(TRD_V, TRD_J) |>
              unique.data.frame()
            if (verbose) message("Writing TRD segments present in the data, but missing in tcrdist3 database to file: ", R.utils::getAbsolutePath(file.path(dirname(outputCsv),'filtered_TRD_gene_segments.csv')))
            utils::write.csv(filtered_genes, file = file.path(dirname(outputCsv),'filtered_TRD_gene_segments.csv'), row.names = FALSE)
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRD_V %in% gene_segments_in_db) |>
          dplyr::filter(TRD_J %in% gene_segments_in_db)
      } else {
        stop(paste0("Chain ", chain, " is not supported."))
      }
    }
    if (verbose) message("Finished metadata cleaning with ", nrow(metadata), " rows")
  }
  #impute clone names if asked
  if (imputeCloneNames) {
    if (verbose) message("Starting clone name imputation with ", nrow(metadata), " rows")

    if (nrow(metadata) == 0) {
      if (verbose) message("ERROR - metadata has 0 rows when trying to impute clone names!")
      stop("No data remaining after filtering. Check your gene segment names and database compatibility.")
    }

    if (!"CloneNames" %in% colnames(metadata)) {
      #initialize the CloneNames metadata column
      if (verbose) message("Creating CloneNames column")
      metadata$CloneNames <- "undefined_clone"
    }
    #assume that clone names are set by Rdiscvr, but if they're NA (like for the tests) we need to impute them

    if (!is.null(spikeInDataframe)){
      #if a user submits a spike-in dataframe, the subject IDs will need to be converted to a character column to merge with SubjectId == "SpikeIn"
      metadata$SubjectId <- as.character(metadata$SubjectId)
    }

    # Check if SubjectId column exists, if not create it
    if (!"SubjectId" %in% colnames(metadata)) {
      if (verbose) message("SubjectId column not found, creating default SubjectId")
      metadata$SubjectId <- "DefaultSubject"
    }

    # Check if TRA column exists, if not create dummy column
    if (!"TRA" %in% colnames(metadata)) {
      if (verbose) message("TRA column not found, creating dummy TRA column")
      metadata$TRA <- "DUMMY_TRA"
    }

    # Check if TRB column exists, if not create dummy column
    if (!"TRB" %in% colnames(metadata)) {
      if (verbose) message("TRB column not found, creating dummy TRB column")
      metadata$TRB <- "DUMMY_TRB"
    }

    # Check if TRA_V column exists, if not create dummy column
    if (!"TRA_V" %in% colnames(metadata)) {
      if (verbose) message("TRA_V column not found, creating dummy TRA_V column")
      metadata$TRA_V <- "DUMMY_TRA_V"
    }

    # Check if TRA_J column exists, if not create dummy column
    if (!"TRA_J" %in% colnames(metadata)) {
      if (verbose) message("TRA_J column not found, creating dummy TRA_J column")
      metadata$TRA_J <- "DUMMY_TRA_J"
    }

    #construct the grouping columns based on the supplied chains
    grouping_columns <- c("SubjectId")
    if ("TRA" %in% chains) {
      grouping_columns <- c(grouping_columns, "TRA", "TRA_V", "TRA_J")
    }
    if ("TRB" %in% chains) {
      grouping_columns <- c(grouping_columns, "TRB", "TRB_V", "TRB_J")
    }
    if ("TRG" %in% chains) {
      grouping_columns <- c(grouping_columns, "TRG", "TRG_V", "TRG_J")
    }
    if ("TRD" %in% chains) {
      grouping_columns <- c(grouping_columns, "TRD", "TRD_V", "TRD_J")
    }

    if (verbose) {
      print(head(metadata))
    }

    metadata <- metadata |>
      #if a user submits a spike-in dataframe, the clones will be missing a subject Id
      dplyr::mutate(SubjectId = dplyr::case_when(
        is.na(CloneNames) & is.na(SubjectId) ~ "SpikeIn",
        TRUE ~ as.character(SubjectId)
      )) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_columns))) |>
      dplyr::mutate(CloneNames =
                      dplyr::case_when( is.na(CloneNames) ~ paste0(SubjectId, "_", dplyr::cur_group_id()),
                                        TRUE ~ as.character(CloneNames)))
  }

  #TODO: this implementation only works for TRA+TRB, need to fix eventually.

  if (summarizeClones) {
    #TODO: figure out if we need to index clones jointly (across both chains)
    #or singly (TRAs would have a clone ID and TRBs would have their own clone ID)
    metadata <- metadata |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grouping_columns))) |>
      dplyr::reframe(count = dplyr::n(), CloneNames) |>
      unique.data.frame()
    #filter out unique/rare clones
    #TODO: check for parameter nesting between summarizeClones and imputeCloneNames appropriately
    if (minimumClonesPerSubject > 1) {
      metadata <- metadata |>
        dplyr::filter(count >= minimumClonesPerSubject)
    }
  }

  #unique-ify the metadata prior to formatting, since the random sampling in the reverse translation function can cause duplicates
  metadata <- metadata |> unique.data.frame()
  if (verbose) message("Final metadata dimensions after cleaning and summarizing: ", nrow(metadata), " rows, ", ncol(metadata), " columns")
  #reformat data for tcrDist3, iterate over chains specified:

  formatted_data <- data.frame( subject = metadata$SubjectId,
                                epitope = rep("dummy_epitope", nrow(metadata)),  #epitope information is not available in seuratMetadata.csv
                                count = if("count" %in% colnames(metadata)) metadata$count else rep(1, nrow(metadata))
                                )
  if (verbose) {
    print(head(formatted_data))
    message("Dimensions of formatted_data: ", paste(dim(formatted_data), collapse = " x "))
    message("Dimensions of metadata: ", paste(dim(metadata), collapse = " x "))
  }
  if ("TRA" %in% chains) {
    if (verbose) message("Formatting TRA chains for tcrDist3")
    formatted_data <- cbind(formatted_data,
                            data.frame(
                              v_a_gene = paste0(metadata$TRA_V, "*01"),
                              j_a_gene = paste0(metadata$TRA_J, "*01"),
                              cdr3_a_aa = metadata$TRA,
                              cdr3_a_nucseq = sapply(metadata$TRA, .reverse_translate_cdr3)
                            )
    )
  }
  if ("TRB" %in% chains) {
    if (verbose) message("Formatting TRB chains for tcrDist3")
    formatted_data <- cbind(formatted_data,
                            data.frame(
                              v_b_gene = paste0(metadata$TRB_V, "*01"),
                              j_b_gene = paste0(metadata$TRB_J, "*01"),
                              cdr3_b_aa = metadata$TRB,
                              cdr3_b_nucseq = sapply(metadata$TRB, .reverse_translate_cdr3)
                            )
    )
  }
  if ("TRG" %in% chains) {
    if (verbose) message("Formatting TRG chains for tcrDist3")
    formatted_data <- cbind(formatted_data,
                            data.frame(
                              v_g_gene = paste0(metadata$TRG_V, "*01"),
                              j_g_gene = paste0(metadata$TRG_J, "*01"),
                              cdr3_g_aa = metadata$TRG,
                              cdr3_g_nucseq = sapply(metadata$TRG, .reverse_translate_cdr3)
                            )
    )
  }
  if ("TRD" %in% chains) {
    if (verbose) message("Formatting TRD chains for tcrDist3")
    formatted_data <- cbind(formatted_data,
                            data.frame(
                              v_d_gene = paste0(metadata$TRD_V, "*01"),
                              j_d_gene = paste0(metadata$TRD_J, "*01"),
                              cdr3_d_aa = metadata$TRD,
                              cdr3_d_nucseq = sapply(metadata$TRD, .reverse_translate_cdr3)
                            )
    )
  }
  #write the formatted data to the output CSV file
  utils::write.csv(formatted_data, outputCsv, row.names = FALSE)
  return(formatted_data)
}

.reverse_translate_cdr3 <- function(cdr3_aa_seq) {
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
    sample(codon_table[[aa]], 1)
  })

  #paste the sampled codons together
  return(paste(cdr3_nuc_seq, collapse = ""))
}

.PullTcrdist3Db <- function(organism = 'human',
                            outputFilePath = './tcrdist3_gene_segments',
                            pythonExecutable = NULL,
                            verbose = FALSE) {
  if (is.null(pythonExecutable)) {
    pythonExecutable <- reticulate::py_exe()
    #fallback if reticulate fails
    if (is.null(pythonExecutable) || pythonExecutable == "") {
      pythonExecutable <- Sys.which("python3")
    }
  }
  outputFilePath <- R.utils::getAbsolutePath(outputFilePath)
  template <- readr::read_file(system.file("scripts/PullTcrdist3Db.py", package = "tcrClustR"))
  script <- tempfile(fileext = ".py")
  readr::write_file(template, script)
  #format and write the Python function call to the script
  command <- paste0("PullTcrdist3Db(organism = '", organism,
                    "', outputFilePath = '", outputFilePath,
                    "')")
  readr::write_file(command, script, append = TRUE)
  #add execution permissions to script and parent directory
  Sys.chmod(script, mode = "755")
  system(paste("chmod 755", dirname(script)))
  #execute
  if (verbose) message("Python executable: ", pythonExecutable)
  result <- system2(pythonExecutable, script, stdout = TRUE, stderr = TRUE)
  if (verbose) cat(result) 
  #check that the gene segments file is created
  if (!file.exists(outputFilePath)) {
    stop("tcrdist3_gene_segments.txt generation failed. Check Python script execution.")
  }
}

#' Create a ComplexHeatmap for a single TCR distance assay
#'
#' @description
#' Internal utility to generate a ComplexHeatmap for a single TCR distance assay.
#'
#' @param seuratObj_TCR A Seurat object containing TCR distance assays.
#' @param assay Character string specifying the assay name to plot.
#' @param cluster_info Factor vector of cluster assignments for each cell in the assay.
#' @param cluster_colors Named character vector of colors corresponding to cluster levels.
#' @param annotate_clusters Boolean specifying whether to display clustering information.
#'
#' @return A ComplexHeatmap object ready for drawing.
#' @importFrom ComplexHeatmap Heatmap draw HeatmapAnnotation rowAnnotation
#' @keywords internal
#' @examples
#' \dontrun{
#' .TCRDistanceHeatmap(
#'   seuratObj_TCR = seuratObj,
#'   assay = "TCRAssay",
#'   cluster_info = info_df,
#'   cluster_colors = color_list,
#'   annotate_clusters = TRUE
#' )
#' }
.TCRDistanceHeatmap <- function(
    seuratObj_TCR = NULL,
    assay = NULL,
    cluster_info = NULL,
    cluster_colors = NULL,
    annotate_clusters = TRUE
) {
  if (is.null(seuratObj_TCR)) {
    stop("seuratObj_TCR must not be NULL.")
  }
  if (is.null(assay)) {
    stop("assay must not be NULL.")
  }
  if (is.null(cluster_info)) {
    stop("cluster_info must not be NULL.")
  }
  if (is.null(cluster_colors)) {
    stop("cluster_colors must not be NULL.")
  }
  m <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))
  cluster_info <- factor(cluster_info)

  if (annotate_clusters) {
    col_annotation <- ComplexHeatmap::HeatmapAnnotation(
      cluster = cluster_info,
      col = list(cluster = cluster_colors),
      show_annotation_name = FALSE,
      annotation_name_side = "left",
      which = "column",
      show_legend = TRUE
    )

    row_annotation <- ComplexHeatmap::rowAnnotation(
      cluster = cluster_info,
      col     = list(cluster = cluster_colors),
      show_annotation_name = FALSE,
      show_legend          = FALSE
    )

    ComplexHeatmap::Heatmap(
      m,
      name               = assay,
      column_title       = assay,
      border_gp          = gpar(col = "black", lty = 2),
      top_annotation     = col_annotation,
      left_annotation    = row_annotation,
      use_raster         = TRUE,
      cluster_columns    = TRUE,
      cluster_rows       = TRUE,
      row_dend_side      = "left",
      column_dend_side   = "top",
      column_split       = cluster_info,
      row_split          = cluster_info,
      show_heatmap_legend = TRUE,
      show_column_names    = FALSE,
      show_row_names       = FALSE
    )
  } else {
    ComplexHeatmap::Heatmap(
      m,
      name               = assay,
      column_title       = assay,
      border_gp          = gpar(col = "black", lty = 2),
      use_raster         = TRUE,
      cluster_columns    = FALSE,
      cluster_rows       = FALSE,
      show_heatmap_legend = TRUE,
      show_column_names    = FALSE,
      show_row_names       = FALSE
    )
  }
}

#' Generate and display heatmaps of TCR similarity for each assay
#'
#' @description
#' Generates and displays ComplexHeatmaps of TCR similarity for each assay in a Seurat object.
#'
#' @param seuratObj_TCR A Seurat object containing one or more TCR distance assays.
#' @param assayList Character vector of assay names to include. Default is NULL (all assays).
#' @param resolution Numeric clustering resolution parameter matching metadata column suffix.
#' @param annotate_clusters Boolean specifying whether to display clustering information.
#'
#' @return Returns the patchwork object containing ComplexHeatmaps of the TCR distance data.
#' @export
#' @examples
#' \dontrun{
#' TCRDistanceHeatmaps(
#'   seuratObj_TCR = seuratObj,
#'   assayList = NULL,
#'   resolution = 0.1,
#'   annotate_clusters = TRUE
#' )
#' }
TCRDistanceHeatmaps <- function(
    seuratObj_TCR = NULL,
    assayList = NULL,
    resolution = 0.1,
    annotate_clusters = TRUE
) {
  if (is.null(seuratObj_TCR)) {
    stop("Please provide a Seurat Object with TCR distance assays.")
  }

  assays_to_use <- if (is.null(assayList)) {
    SeuratObject::Assays(seuratObj_TCR)
  } else {
    assayList
  }

  heatmaps <- list()

  for (assay in assays_to_use) {
    distance_matrix <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))

    # Find metadata column from clustering pipeline
    cluster_column <- paste0("TcrClustR_", assay, "_", resolution)
    if (!(cluster_column %in% colnames(seuratObj_TCR@meta.data))) {
      message(paste("Skipping", assay, "- no", cluster_column))
      next
    }

    # Cluster metadata is an array of assay_size x n_assay, so iterate and slice
    # the correct portion of that array
    full_cluster_info <- seuratObj_TCR@meta.data[[cluster_column]]
    n_cells_in_assay <- ncol(distance_matrix)
    assay_start_index <- 1
    for (a in SeuratObject::Assays(seuratObj_TCR)) {
      if (a == assay) break
      assay_start_index <- assay_start_index + ncol(Seurat::GetAssayData(seuratObj_TCR, assay = a, layer = "counts"))
    }
    assay_end_index <- assay_start_index + n_cells_in_assay - 1

    cluster_info <- full_cluster_info[assay_start_index:assay_end_index]
    cluster_info <- as.factor(cluster_info)
    cluster_levels <- levels(cluster_info)
    cluster_colors <- setNames(
      if (length(cluster_levels) <= 8) {
        RColorBrewer::brewer.pal(length(cluster_levels), "Set2")
      } else if (length(cluster_levels) <= 12) {
        RColorBrewer::brewer.pal(length(cluster_levels), "Set3")
      } else {
        grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Set3"))(length(cluster_levels))
      },
      cluster_levels
    )

    # Get a ComplexHeatmap
    heatmap_obj <- .TCRDistanceHeatmap(seuratObj_TCR, assay, cluster_info, cluster_colors, annotate_clusters)
    drawn_heatmap <- draw(heatmap_obj, merge_legend = FALSE, heatmap_legend_side = "right", annotation_legend_side = "right", newpage = FALSE)

    if (!is.null(drawn_heatmap)) {
      heatmaps[[assay]] <- drawn_heatmap
    }
  }

  # Composite with patchwork
  combined_heatmaps <- patchwork::wrap_plots(lapply(heatmaps, function(hm) grid.grabExpr(draw(hm))), ncol = 1)

  final_plot <- combined_heatmaps +
    patchwork::plot_annotation(title = "TCR Similarity",
                               theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

  print(final_plot)

  return(final_plot)
}

#' Plot histograms of summed TCR distances by cluster for each assay
#'
#' @description
#' Plots histograms of summed TCR distances by cluster for each assay in a Seurat object.
#'
#' @param seuratObj_TCR A Seurat object with TCR distance assay data.
#' @param assayList Character vector of assays to plot. Default is NULL (all assays).
#' @param resolution Numeric clustering resolution matching metadata column suffix.
#'
#' @return Returns the patchwork object containing the TCR distance histograms.
#' @export
#' @examples
#' \dontrun{
#' TCRDistanceHistograms(
#'   seuratObj_TCR = seuratObj,
#'   assayList = NULL,
#'   resolution = 0.1
#' )
#' }
TCRDistanceHistograms <- function(
    seuratObj_TCR = NULL,
    assayList  = NULL,
    resolution = 0.1
) {
  if (is.null(seuratObj_TCR)) {
    stop("Please provide a Seurat object with TCR distance assays.")
  }

  assays <- if (is.null(assayList)) {
    SeuratObject::Assays(seuratObj_TCR)
  } else {
    assayList
  }

  # precompute how many cells each assay has
  cell_counts <- setNames(
    vapply(assays, function(a) {
      ncol(Seurat::GetAssayData(seuratObj_TCR, assay = a, layer = "counts"))
    }, integer(1)),
    assays
  )

  # compute start/end indices for slicing the metadata vector
  starts <- cumsum(c(1, head(cell_counts, -1)))
  ends   <- cumsum(cell_counts)

  plots <- list()
  for (i in seq_along(assays)) {
    assay <- assays[i]
    cluster_col <- paste0("TcrClustR_", assay, "_", resolution)

    if (!(cluster_col %in% colnames(seuratObj_TCR@meta.data))) {
      message("Skipping ", assay, ": no metadata column ", cluster_col)
      next
    }

    dist_mat  <- as.matrix(Seurat::GetAssayData(
      seuratObj_TCR, assay = assay, layer = "counts"
    ))

    # slice the full clustering vector down to this assay
    full_info    <- seuratObj_TCR@meta.data[[cluster_col]]
    cluster_info <- factor(full_info[starts[i]:ends[i]])

    # Build palette
    n_clust <- length(levels(cluster_info))
    pal     <- RColorBrewer::brewer.pal(min(n_clust, 8), "Set2")
    cl_cols <- setNames(pal, levels(cluster_info))

    df <- data.frame(
      DistanceSum = rowSums(dist_mat),
      Cluster     = cluster_info
    )

    p <- ggplot(df, aes(x = DistanceSum, fill = Cluster)) +
      geom_histogram(bins = 50, color = "black") +
      scale_fill_manual(values = cl_cols) +
      facet_wrap(~ Cluster, nrow = 1, scales = "free_y") +
      labs(
        title = assay,
        x     = "Summed Distances",
        y     = "# of cells"
      ) +
      theme_minimal() +
      theme(
        plot.title      = element_text(size = 14, face = "bold", hjust = 0),
        axis.title      = element_text(size = 12),
        axis.text       = element_text(size = 10),
        strip.text      = element_text(face = "bold"),
        legend.position = "none"
      )

    plots[[assay]] <- p
  }

  if (length(plots) == 0) {
    stop("No assays with the requested clustering metadata were found.")
  }

  combined <- patchwork::wrap_plots(plots, ncol = 1)
  print(combined)

  return(combined)
}
