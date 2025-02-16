
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
#' @return NULL
#' @export

FormatMetadataForTcrDist3 <- function(metadata,
                                      outputCsv = './tcrDist3Input.csv',
                                      chains = c("TRA", "TRB"),
                                      organism = 'human',
                                      cleanMetadata = T,
                                      summarizeClones = T,
                                      imputeCloneNames = T,
                                      minimumClonesPerSubject = 100,
                                      writeUnannotatedGeneSegmentsToFile = T
) {
  if (cleanMetadata) {
    for (chain in chains) {
      #filter rows with NA values in the requested chains
      metadata <- metadata[!is.na(metadata[[chain]]), ]
      #filter rows with commas (multiple segments detected in a cell) in the requested chains
      metadata <- metadata[!grepl(",", metadata[[chain]]), ]

      .PullTcrdist3Db(organism = organism,
                      outputFilePath = './tcrdist3_gene_segments.txt')
      gene_segments_in_db <- readr::read_csv('./tcrdist3_gene_segments.txt', show_col_types = FALSE) |>
        dplyr::mutate(`gene_segments` = gsub("\\*[0-9]+$", "", `gene_segments`)) |>
        unlist() |>
        unique()
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
            print(paste0("Writing TRA segments present in the data, but missing in tcrdist3 database to file: ", R.utils::getAbsolutePath('./filtered_TRA_gene_segments.csv')))
            utils::write.csv(filtered_genes, file = './filtered_TRA_gene_segments.csv', row.names = FALSE)
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRA_V %in% gene_segments_in_db) |>
          dplyr::filter(TRA_J %in% gene_segments_in_db)
      } else if (chain == "TRB") {
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
            print(paste0("Writing TRB segments present in the data, but missing in tcrdist3 database to file: ", R.utils::getAbsolutePath('./filtered_TRB_gene_segments.csv')))
            utils::write.csv(filtered_genes, file = './filtered_TRB_gene_segments.csv', row.names = FALSE)
          }
        }
        metadata <- metadata |>
          dplyr::filter(TRB_V %in% gene_segments_in_db) |>
          dplyr::filter(TRB_J %in% gene_segments_in_db)

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
            print(paste0("Writing TRG segments present in the data, but missing in tcrdist3 database to file: ", R.utils::getAbsolutePath('./filtered_TRG_gene_segments.csv')))
            utils::write.csv(filtered_genes, file = './filtered_TRG_gene_segments.csv', row.names = FALSE)
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
            print(paste0("Writing TRD segments present in the data, but missing in tcrdist3 database to file: ", R.utils::getAbsolutePath('./filtered_TRD_gene_segments.csv')))
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
  #impute clone names if asked
  if (imputeCloneNames) {
    if (!"CloneNames" %in% colnames(metadata)) {
      #initialize the CloneNames metadata column
      metadata$CloneNames <- NA
    }
    #assume that clone names are set by Rdiscvr, but if they're NA (like for the tests, we need to impute them)
    metadata$CloneNames <- ifelse(is.na(metadata$CloneNames),
                                  yes = paste0(metadata$SubjectId, "_", seq_len(nrow(metadata))),
                                  no = metadata$CloneNames)
  }
  if (summarizeClones) {
    #TODO: figure out if we need to index clones jointly (across both chains)
    #or singly (TRAs would have a clone ID and TRBs would have their own clone ID)
    metadata <- metadata |>
      dplyr::group_by(SubjectId, TRA, TRB, TRA_V, TRA_J, TRB_V, TRB_J) |>
      dplyr::reframe(count = dplyr::n(), CloneNames) |>
      unique.data.frame()
    #filter out unique/rare clones
    #TODO: check for parameter nesting between summarizeClones and imputeCloneNames appropriately
    if (minimumClonesPerSubject > 1) {
      metadata <- metadata |>
        dplyr::filter(count >= minimumClonesPerSubject)
    }
  }


  #reformat data
  formatted_data <- data.frame(
    subject = metadata$SubjectId,
    epitope = rep("dummy_epitope", nrow(metadata)),  #epitope information is not available in seuratMetadata.csv
    count = metadata$count,    #TODO: I could tabulate count info if necessary
    v_a_gene = paste0(metadata$TRA_V, "*01"),
    j_a_gene = paste0(metadata$TRA_J, "*01"),
    cdr3_a_aa = metadata$TRA,
    cdr3_a_nucseq = sapply(metadata$TRA, .reverse_translate_cdr3),  #TODO: nucleotides aren't given, but could be reverse-translated
    v_b_gene = paste0(metadata$TRB_V, "*01"),
    j_b_gene = paste0(metadata$TRB_J, "*01"),
    cdr3_b_aa = metadata$TRB,
    cdr3_b_nucseq = sapply(metadata$TRB, .reverse_translate_cdr3),  #TODO: nucleotides aren't given, but could be reverse-translated
    clone_id = metadata$CloneNames
  )
  # Write the formatted data to the output CSV file
  utils::write.csv(formatted_data, outputCsv, row.names = FALSE)
}

.reverse_translate_cdr3 <- function(cdr3_aa_seq) {
  # Codon table
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
                            pythonExecutable = NULL) {
  if (is.null(pythonExecutable)) {
    pythonExecutable <- reticulate::py_exe()
  }
  outputFilePath <- R.utils::getAbsolutePath(outputFilePath)
  template <- readr::read_file(system.file("scripts/PullTcrdist3Db.py", package = "tcrClustR"))
  script <- tempfile()
  readr::write_file(template, script)
  #format and write the python function to the end of the script
  command <- paste0("PullTcrdist3Db(organism = '", organism,
                    "', outputFilePath = '", outputFilePath,
                    "')")
  readr::write_file(command, script, append = TRUE)
  #execute
  system2(pythonExecutable, script)
}
