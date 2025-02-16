
utils::globalVariables(
  names = c('id', 'unsequencedTCRs', 'TRA_V', 'TRA_J', 'organism', 'TRB_V', 'TRB_J', 'TRG_V', 'TRG_J', 'TRD_V', 'TRD_J', 'gene_segments', 'CloneNames', 'count'),
  package = 'tcrClustR',
  add = TRUE
)

#' @title FormatMetadataForConga
#' @description This function formats metadata for the conga pipeline.
#'
#' @param metadata Data frame containing metadata.
#' @param outputCsv Path to the output CSV file.
#' @param organism Organism to use for conga's TCR db. Default is 'human'.
#' @param chains TCR chains to include in the analysis. TRA/TRB supported and tested, but others likely work
#' @param cleanMetadata Boolean controlling whether to clean the metadata by removing rows with NA values or commas in the specified chains.
#' @param minimumClonesPerSubject Minimum number of clones per subject to include in the analysis. Default is 2.
#' @param writeUnannotatedGeneSegmentsToFile Boolean controlling whether to write unannotated gene segments to a file (filtered_(chain)_gene_segments.csv).
#' @param spikeInDataframe Data frame containing spike-in data. Default is NULL.
#' @export

FormatMetadataForConga <- function(metadata,
                                   outputCsv = './congaInput.csv',
                                   organism = 'human',
                                   chains = NULL,
                                   cleanMetadata = T,
                                   minimumClonesPerSubject = 100,
                                   writeUnannotatedGeneSegmentsToFile = T,
                                   spikeInDataframe = NULL) {
  #determine gene segments from chains
  gene_segments_and_chains <- c()
  for (chain in chains) {
    if (chain == "TRA") {
      gene_segments_and_chains <- c(gene_segments_and_chains, "TRA_V", "TRA_J", "TRA")
    } else if (chain == "TRB") {
      gene_segments_and_chains <- c(gene_segments_and_chains, "TRB_V", "TRB_J", "TRB")
    } else if (chain == "TRD") {
      gene_segments_and_chains <- c(gene_segments_and_chains, "TRD_V", "TRD_J", "TRD")
    } else if (chain == "TRG") {
      gene_segments_and_chains <- c(gene_segments_and_chains, "TRG_V", "TRG_J", "TRG")
    } else {
      stop(paste0("Chain ", chain, " is not supported."))
    }
  }
  if (cleanMetadata) {
    for (chain in chains) {
      #filter rows with NA values in the requested chains
      metadata <- metadata[!is.na(metadata[[chain]]), ]
      #filter rows with commas (multiple segments detected in a cell) in the requested chains
      metadata <- metadata[!grepl(",", metadata[[chain]]), ]
      #filter gene segments with commas and NA values
      metadata <- metadata[!grepl(",", metadata[[paste0(chain, "_V")]]), ]
      metadata <- metadata[!grepl(",", metadata[[paste0(chain, "_J")]]), ]
      metadata <- metadata[!is.na(metadata[[paste0(chain, "_V")]]), ]
      metadata <- metadata[!is.na(metadata[[paste0(chain, "_J")]]), ]

      .PullCongaDb(organism = organism,
                   outputFilePath = './conga_gene_segments.txt')
      gene_segments_in_db <- readr::read_csv('./conga_gene_segments.txt', show_col_types = FALSE) |>
        dplyr::mutate(`id` = gsub("\\*[0-9]+$", "", `id`)) |>
        dplyr::pull(id) |>
        unlist() |>
        unique()

      #remove gene segments not found in conga's database
      #TODO: store the gene segments in the data that are not found in the database
      if (chain == "TRA") {
        if (writeUnannotatedGeneSegmentsToFile) {
          #store filtered gene segments
          if (any(!(metadata$TRA_V %in% gene_segments_in_db) | any(!metadata$TRA_J %in% gene_segments_in_db))) {
            message("TRA gene segments present in the data, but not found in conga database!")
            print(paste0("Writing TRA segments present in the data, but missing in conga database to file: ", R.utils::getAbsolutePath('./filtered_TRA_gene_segments.csv')))
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
            utils::write.csv(filtered_genes, file = './filtered_TRA_gene_segments.csv', row.names = FALSE)
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRA_V %in% gene_segments_in_db) |>
          dplyr::filter(TRA_J %in% gene_segments_in_db)
      } else if (chain == "TRB") {
        if (writeUnannotatedGeneSegmentsToFile) {
          #store filtered gene segments
          if (any(!(metadata$TRB_V %in% gene_segments_in_db) | any(!metadata$TRB_J %in% gene_segments_in_db))) {
            message("TRB gene segments present in the data, but not found in conga database!")
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
            utils::write.csv(filtered_genes, file = './filtered_TRB_gene_segments.csv', row.names = FALSE)
            print(paste0("Writing TRB segments present in the data, but missing in conga database to file: ", R.utils::getAbsolutePath('./filtered_TRB_gene_segments.csv')))
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRA_V %in% gene_segments_in_db) |>
          dplyr::filter(TRA_J %in% gene_segments_in_db)
      } else if (chain == "TRG") {
        if (writeUnannotatedGeneSegmentsToFile) {
          #store filtered gene segments
          if (any(!(metadata$TRG_V %in% gene_segments_in_db) | any(!metadata$TRG_J %in% gene_segments_in_db))) {
            message("TRG gene segments present in the data, but not found in conga database!")
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
            print(paste0("Writing TRG segments present in the data, but missing in conga database to file: ", R.utils::getAbsolutePath('./filtered_TRG_gene_segments.csv')))
            utils::write.csv(filtered_genes, file = './filtered_TRG_gene_segments.csv', row.names = FALSE)
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRG_V %in% gene_segments_in_db) |>
          dplyr::filter(TRG_J %in% gene_segments_in_db)
      } else if (chain == "TRD") {
        if (writeUnannotatedGeneSegmentsToFile) {
          #store filtered gene segments
          if (any(!(metadata$TRD_V %in% gene_segments_in_db) | any(!metadata$TRD_J %in% gene_segments_in_db))) {
            message("TRD gene segments present in the data, but not found in conga database!")
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
            print(paste0("Writing TRD segments present in the data, but missing in conga database to file: ", R.utils::getAbsolutePath('./filtered_TRD_gene_segments.csv')))
            utils::write.csv(filtered_genes, file = './filtered_TRD_gene_segments.csv', row.names = FALSE)
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRD_V %in% gene_segments_in_db) |>
          dplyr::filter(TRD_J %in% gene_segments_in_db)
      } else {
        stop(paste0("Chain ", chain, " is not supported."))
      }
    }
  }
  #filter the clones if requested
  if (minimumClonesPerSubject > 1) {
    #filter out unique/rare clones
    metadata <- metadata |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c("SubjectId", gene_segments_and_chains)))) |>
      dplyr::reframe(count = dplyr::n()) |>
      unique.data.frame()
    metadata <- metadata |>
      dplyr::filter(count >= minimumClonesPerSubject)
  }
  #isolate necessary columns and de-duplicate the chains
  metadata <- metadata |>
    dplyr::select(dplyr::all_of(gene_segments_and_chains)) |>
    dplyr::as_tibble() |>
    unique.data.frame()

  if (!is.null(spikeInDataframe)) {
    #check if the spikeInDataframe is a dataframe
    if (!is.data.frame(spikeInDataframe)) {
      stop("The spikeInDataframe must be a dataframe.")
    }
    #check if the spikeInDataframe has the same columns as the metadata
    if (!all(colnames(spikeInDataframe) %in% colnames(metadata))) {
      stop(paste0("The spikeInDataframe must have the same columns as the metadata. These columns are determined by the requested chains: ",
                  chains, ". The columns in the spikeInDataframe are: ", colnames(spikeInDataframe), "."))
    }
    metadata <- rbind(metadata, unsequencedTCRs)
  }
  #normalize the path for outputCsv
  outputCsv <- R.utils::getAbsolutePath(outputCsv)

  #format for conga.
  #conga requires a list of tuples to be literally interpreted by python, so that's how we must write the csv
  #the individual chains are written differently than joint chain couplets, so there's a case for each version.
  if (all(chains == "TRA")) {
    utils::write.csv(file = outputCsv,
                     paste0("[",
                            paste0("(", "('", metadata$TRA_V, "*01','", metadata$TRA_J, "*01','", metadata$TRA, "'))" , collapse = ","),
                            "]"), row.names = FALSE, quote = FALSE)

  } else if (all(chains == "TRB")) {
    utils::write.csv(file = outputCsv,
                     paste0("[",
                            paste0("(", "('", metadata$TRB_V, "*01','", metadata$TRB_J, "*01','", metadata$TRB, "'))" , collapse = ","),
                            "]"), row.names = FALSE, quote = FALSE)

  } else if (all(chains == "TRG")) {
    utils::write.csv(file = outputCsv,
                     paste0("[",
                            paste0("(", "('", metadata$TRG_V, "*01','", metadata$TRG_J, "*01','", metadata$TRG, "'))" , collapse = ","),
                            "]"), row.names = FALSE, quote = FALSE)
  } else if (all(chains == "TRD")) {
    utils::write.csv(file = outputCsv,
                     paste0("[",
                            paste0("(", "('", metadata$TRD_V, "*01','", metadata$TRD_J, "*01','", metadata$TRD, "'))" , collapse = ","),
                            "]"), row.names = FALSE, quote = FALSE)
  } else {
    stop(paste0("Chain ", chains, " is not supported."))
  }
}

.PullCongaDb <- function(organism = 'human',
                         outputFilePath = './conga_gene_segments.txt',
                         pythonExecutable = NULL) {
  if (is.null(pythonExecutable)) {
    pythonExecutable <- reticulate::py_exe()
  }
  outputFilePath <- R.utils::getAbsolutePath(outputFilePath)
  template <- readr::read_file(system.file("scripts/PullCongaDb.py", package = "tcrClustR"))
  script <- tempfile()
  readr::write_file(template, script)
  #format and write the python function to the end of the script
  command <- paste0("PullCongaDb(organism = '", organism,
                    "', outputFilePath = '", outputFilePath,
                    "')")
  readr::write_file(command, script, append = TRUE)
  #execute
  system2(pythonExecutable, script)
}
