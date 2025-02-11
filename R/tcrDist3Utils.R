
utils::globalVariables(
  names = c('SubjectId', 'TRA_V', 'TRA_J', 'TRB_V', 'TRB_J', 'CloneNames', 'count'),
  package = 'tcrClustR',
  add = TRUE
)

#' Format metadata for tcrDist3
#' @description This function formats a seurat object's metadata (with TCR information appended) for tcrDist3 distance caluclations.
#'
#' @param metadata Data frame containing metadata.
#' @param chains TCR chains to include in the analysis. TRA/TRB supported and tested, but others likely work
#' @param cleanMetadata Boolean controlling whether to clean the metadata by removing rows with NA values or commas in the specified chains.
#' @param summarizeClones Boolean controlling whether to summarize clones by SubjectId, TRA, TRB, TRA_V, TRA_J, TRB_V, and TRB_J.
#' @param imputeCloneNames Boolean controlling whether to impute clone names if they are missing.
#' @param outputCsv Path to the output CSV file.
#' @return NULL
#' @export

FormatMetadataForTcrDist3 <- function(metadata,
                                      outputCsv = './tcrDist3Input.csv',
                                      chains = c("TRA", "TRB"),
                                      cleanMetadata = T,
                                      summarizeClones = T,
                                      imputeCloneNames = T,
                                      minimumClonesPerSubject = 100
                                      ) {
  if (cleanMetadata) {
    for (chain in chains) {
      #filter rows with NA values in the requested chains
      metadata <- metadata[!is.na(metadata[[chain]]), ]
      #filter rows with commas (multiple chains detected in a cell) in the requested chains
      metadata <- metadata[!grepl(",", metadata[[chain]]), ]
    }
  }

  if (summarizeClones) {
    #TODO: figure out if we need to index clones jointly (across both chains)
    #or singly (TRAs would have a clone ID and TRBs would have their own clone ID)
    metadata <- metadata |>
      dplyr::group_by(SubjectId, TRA, TRB, TRA_V, TRA_J, TRB_V, TRB_J) |>
      dplyr::reframe(count = dplyr::n(), CloneNames) |>
      unique.data.frame()
    #filter out unique/rare clones
    #TODO: check for parameter nesting appropriately
    if (minimumClonesPerSubject > 1) {
      metadata <- metadata |>
        dplyr::filter(count >= minimumClonesPerSubject)
    }
  }

  if (imputeCloneNames) {
    #assume that clone names are set by Rdiscvr, but if they're NA (like for the tests, we need to impute them)
    metadata$CloneNames <- ifelse(is.na(metadata$CloneNames),
                                  yes = paste0(metadata$SubjectId, "_", seq_len(nrow(metadata))),
                                  no = metadata$CloneNames)
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
  write.csv(formatted_data, outputCsv, row.names = FALSE)
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
