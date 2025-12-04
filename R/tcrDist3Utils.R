utils::globalVariables(
  names = c('SubjectId', 'TRA', 'TRA_V', 'TRA_J', 'TRB', 'TRB_V', 'TRB_J', 'CloneNames', 'count', 'Cluster', 'DistanceSum'),
  package = 'tcrClustR',
  add = TRUE
)

#' Format metadata for tcrDist3
#' @description This function formats a seurat object's metadata (with TCR information appended) for tcrDist3 distance caluclations.
#'
#' @param metadata Data frame containing metadata.
#' @param chains TCR chains to include in the analysis. TRA/TRB supported and tested, but others likely work.
#' @param organism Organism to use for tcrDist3. Default is 'human'.
#' @param sampleGroupingColumns The set of columns on which to group the data
#' @param calculateChainPairs If true, this will prepare the columns needed for A/B and/or  G/D (depending on values of chains
#' @param spikeInDataframe Data frame containing spike-in data. Default is NULL. See examples for formatting requirements.
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

FormatMetadataForTcrDist3 <- function(metadata,
                                      chains = c("TRA", "TRB"),
                                      organism = 'human',
                                      sampleGroupingColumns = c('SubjectId'),
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

  if (any(!sampleGroupingColumns %in% colnames(metadata))) {
    missing <- sampleGroupingColumns[!sampleGroupingColumns %in% colnames(metadata)]
    stop(paste0('Metadata missing columns: ', paste0(missing, collapse = ',')))
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
    for (chainName in chains) {
      if (! .has_chain_columns(spikeInDataframe, chain)) {
        stop(paste0("The spikeInDataframe must have the columns '", chain, "' (the CDR3), '", paste0(chain, '_V'), "', and '", paste0(chain, '_J'), "'"))
      }
    }

    # Check that the spikeInDataframe has the columns 'CloneNames'
    if (!"CloneNames" %in% colnames(spikeInDataframe)) {
      stop("The spikeInDataframe must have the column 'CloneNames'")
    }

    for (colName in sampleGroupingColumns) {
      # Ensure the source is a string for rbind
      metadata[[colName]] <- as.character(metadata[[colName]])

      if (!colName %in% names(spikeInDataframe)) {
        spikeInDataframe[[colName]] <- paste0("SpikeIn_", seq_len(nrow(spikeInDataframe)))
      }
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

  # Group and add unique columns for later joins:
  metadata$GlobalRowIdx <- seq_len(nrow(metadata))
  metadata <- metadata |>
    tibble::rownames_to_column('RowNames_') |>
    dplyr::group_by(dplyr::across(dplyr::all_of(sampleGroupingColumns))) |>
    dplyr::mutate(SampleIdx = paste0('Sample_', dplyr::cur_group_id())) |>
    dplyr::ungroup() %>%
    tibble::column_to_rownames('RowNames_')

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
      dplyr::mutate(`_CloneIdx_` = paste0(chain, '-', 'Clone-', dplyr::cur_group_id())) |>
      dplyr::ungroup() %>%
      tibble::column_to_rownames('RowNames_')

    metadata[['_CloneIdx_']][!metadata$IsValidForChain] <- NA

    # There might be a more elegant solution to this:
    names(metadata)[names(metadata) == '_CloneIdx_'] <- paste0(chain, '-CloneIdx')
    metadata$IsValidForChain <- NULL
  }

  # Repeat for A/B and G/D:
  if (calculateChainPairs) {
    toTest <- c()
    if ('TRA' %in% chains && 'TRB' %in% chains) {
      toTest <- c(toTest, 'TRA_TRB')
    }

    if ('TRG' %in% chains && 'TRD' %in% chains) {
      toTest <- c(toTest, 'TRG_TRD')
    }

    for (chainId in toTest) {
      chains <- unlist(strsplit(chainId, split = '_'))
      metadata$IsValidForChain <- metadata[[paste0(chains[1], '_ValidForClustering')]] & metadata[[paste0(chains[2], '_ValidForClustering')]]
      tcr_grouping_columns <- c(chains[1], paste0(chains[1], "_V"), paste0(chains[1], "_J"), chains[2], paste0(chains[2], "_V"), paste0(chains[2], "_J"))
      metadata <- metadata |>
        tibble::rownames_to_column('RowNames_') |>
        dplyr::group_by(dplyr::across(dplyr::all_of(tcr_grouping_columns))) |>
        dplyr::mutate(`_CloneIdx_` = paste0(chainId, '-', 'Clone-', dplyr::cur_group_id())) |>
        dplyr::ungroup() %>%
        tibble::column_to_rownames('RowNames_')

      metadata[['_CloneIdx_']][!metadata$IsValidForChain] <- NA

      # There might be a more elegant solution to this:
      names(metadata)[names(metadata) == '_CloneIdx_'] <- paste0(chainId, '-CloneIdx')
      metadata$IsValidForChain <- NULL
    }
  }

  if (any(rownames(metadata) != metadata$RowName)) {
    print(head(rownames(metadata)))
    stop('Row names did not match after FormatMetadataForTcrDist3')
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

  if (pythonExecutable == "" || !file.exists(pythonExecutable)) {
    stop("No valid Python executable found. Please install Python 3.8+ or run SetupPythonEnvironment() to configure.")
  }

  #validate Python modules are available before proceeding
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
         "\n\nFor tcrdist3, use: ", pythonExecutable, " -m pip install git+https://github.com/kmayerb/tcrdist3.git@0.2.2")
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
#' @param verbose Boolean controlling whether to display processing steps. Default is FALSE.
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
#'   annotate_clusters = TRUE,
#'   verbose = TRUE
#' )
#' }
.TCRDistanceHeatmap <- function(
  seuratObj_TCR = NULL,
  assay = NULL,
  cluster_info = NULL,
  cluster_colors = NULL,
  annotate_clusters = TRUE,
  verbose = FALSE
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

  if (verbose) message("Creating heatmap for assay: ", assay)

  m <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))

  if (verbose) {
    message("Distance matrix dimensions: ", nrow(m), " x ", ncol(m))
    message("Number of cluster levels: ", length(levels(factor(cluster_info))))
    message("Cluster info length: ", length(cluster_info))
  }

  # Validate dimensions
  if (nrow(m) == 0 || ncol(m) == 0) {
    stop("Distance matrix has zero dimensions for assay ", assay, ": ", nrow(m), " x ", ncol(m))
  }

  if (length(cluster_info) != ncol(m)) {
    stop("Cluster info length (", length(cluster_info), ") does not match matrix columns (", ncol(m), ") for assay ", assay)
  }

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

    heatmap <- ComplexHeatmap::Heatmap(
      m,
      name               = assay,
      column_title       = assay,
      border_gp          = grid::gpar(col = "black", lty = 2),
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
    heatmap <- ComplexHeatmap::Heatmap(
      m,
      name               = assay,
      column_title       = assay,
      border_gp          = grid::gpar(col = "black", lty = 2),
      use_raster         = TRUE,
      cluster_columns    = FALSE,
      cluster_rows       = FALSE,
      show_heatmap_legend = TRUE,
      show_column_names    = FALSE,
      show_row_names       = FALSE
    )
  }
  return(heatmap)
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
#' @param verbose Boolean controlling whether to display processing steps. Default is FALSE.
#'
#' @return Returns the patchwork object containing ComplexHeatmaps of the TCR distance data.
#' @export
#' @examples
#' \dontrun{
#' TCRDistanceHeatmaps(
#'   seuratObj_TCR = seuratObj,
#'   assayList = NULL,
#'   resolution = 0.1,
#'   annotate_clusters = TRUE,
#'   verbose = TRUE
#' )
#' }
TCRDistanceHeatmaps <- function(
  seuratObj_TCR = NULL,
  assayList = NULL,
  resolution = 0.1,
  annotate_clusters = TRUE,
  verbose = FALSE
) {
  if (is.null(seuratObj_TCR)) {
    stop("Please provide a Seurat Object with TCR distance assays.")
  }

  if (verbose) {
    message("Starting TCRDistanceHeatmaps with resolution: ", resolution)
    message("Available assays in Seurat object: ", paste(SeuratObject::Assays(seuratObj_TCR), collapse = ", "))
    message("Available metadata columns: ", paste(colnames(seuratObj_TCR@meta.data), collapse = ", "))
  }

  assays_to_use <- if (is.null(assayList)) {
    SeuratObject::Assays(seuratObj_TCR)
  } else {
    assayList
  }

  if (verbose) message("Assays to process: ", paste(assays_to_use, collapse = ", "))

  heatmaps <- list()

  for (assay in assays_to_use) {
    if (verbose) message("Processing assay: ", assay)

    distance_matrix <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))

    if (verbose) message("Distance matrix dimensions for ", assay, ": ", nrow(distance_matrix), " x ", ncol(distance_matrix))

    # Find metadata column from clustering pipeline
    cluster_column <- paste0("TcrClustR_", assay, "_", resolution)
    if (verbose) message("Looking for cluster column: ", cluster_column)

    if (!(cluster_column %in% colnames(seuratObj_TCR@meta.data))) {
      if (verbose) message("Skipping ", assay, " - cluster column ", cluster_column, " not found")
      message(paste("Skipping", assay, "- no", cluster_column))
      next
    }

    # Cluster metadata is an array of assay_size x n_assay, so iterate and slice
    # the correct portion of that array
    full_cluster_info <- seuratObj_TCR@meta.data[[cluster_column]]
    n_cells_in_assay <- ncol(distance_matrix)

    if (verbose) {
      message("Full cluster info length: ", length(full_cluster_info))
      message("Cells in current assay (", assay, "): ", n_cells_in_assay)
    }

    assay_start_index <- 1
    for (a in SeuratObject::Assays(seuratObj_TCR)) {
      if (a == assay) break
      assay_cells <- ncol(Seurat::GetAssayData(seuratObj_TCR, assay = a, layer = "counts"))
      if (verbose) message("Assay ", a, " has ", assay_cells, " cells, start index now: ", assay_start_index + assay_cells)
      assay_start_index <- assay_start_index + assay_cells
    }
    assay_end_index <- assay_start_index + n_cells_in_assay - 1

    if (verbose) message("Slicing cluster info from index ", assay_start_index, " to ", assay_end_index)

    if (assay_end_index > length(full_cluster_info)) {
      if (verbose) message("ERROR: End index (", assay_end_index, ") exceeds cluster info length (", length(full_cluster_info), ")")
      message("Skipping ", assay, " - cluster indexing error")
      next
    }

    cluster_info <- full_cluster_info[assay_start_index:assay_end_index]
    cluster_info <- as.factor(cluster_info)
    cluster_levels <- levels(cluster_info)

    if (verbose) {
      message("Cluster levels for ", assay, ": ", paste(cluster_levels, collapse = ", "))
      message("Number of cluster levels: ", length(cluster_levels))
    }
    cluster_colors <- stats::setNames(
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
    if (verbose) message("Creating ComplexHeatmap for ", assay)

    heatmap_obj <- .TCRDistanceHeatmap(seuratObj_TCR,
                                       assay,
                                       cluster_info,
                                       cluster_colors,
                                       annotate_clusters,
                                       verbose = verbose)
    drawn_heatmap <- draw(heatmap_obj,
                          merge_legend = FALSE,
                          heatmap_legend_side = "right",
                          annotation_legend_side = "right",
                          newpage = FALSE)

    if (!is.null(drawn_heatmap)) {
      if (verbose) message("Successfully created heatmap for ", assay)
      heatmaps[[assay]] <- drawn_heatmap
    } else {
      if (verbose) message("Failed to create heatmap for ", assay)
    }
  }

  if (verbose) message("Created ", length(heatmaps), " heatmaps total")

  # Composite with patchwork
  combined_heatmaps <- patchwork::wrap_plots(lapply(heatmaps, function(hm) grid::grid.grabExpr(draw(hm))), ncol = 1)

  final_plot <- combined_heatmaps +
    patchwork::plot_annotation(title = "TCR Similarity",
                               theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5)))

  #print(final_plot)

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
#' @importFrom stats setNames
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
  cell_counts <- stats::setNames(
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
    cl_cols <- stats::setNames(pal, levels(cluster_info))

    df <- data.frame(
      DistanceSum = rowSums(dist_mat),
      Cluster     = cluster_info
    )

    p <- ggplot2::ggplot(df, ggplot2::aes(x = DistanceSum, fill = Cluster)) +
      ggplot2::geom_histogram(bins = 50, color = "black") +
      ggplot2::scale_fill_manual(values = cl_cols) +
      ggplot2::facet_wrap(~ Cluster, nrow = 1, scales = "free_y") +
      ggplot2::labs(
        title = assay,
        x     = "Summed Distances",
        y     = "# of cells"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title      = ggplot2::element_text(size = 14, face = "bold", hjust = 0),
        axis.title      = ggplot2::element_text(size = 12),
        axis.text       = ggplot2::element_text(size = 10),
        strip.text      = ggplot2::element_text(face = "bold"),
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

#' @title Get Example Markdown
#'
#' @description Save a template R markdown file, showing usage of this package
#' @param dest The destination filepath, where the file will be saved
#' @export
GetExampleMarkdown <- function(dest) {
  source <- system.file("rmd/tcrClustR.rmd", package = "tcrClustR")
  if (!file.exists(source)) {
    stop(paste0('Unable to find file: ', source))
  }

  dest <- normalizePath(dest, mustWork = F)
  success <- file.copy(source, dest, overwrite = TRUE)

  if (!success) {
    stop(paste0('Unable to copy file to: ', dest))
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
      stop(paste0('Missing columns: ', paste0(missingCols, collapse = ','), ', available columns: ', paste(colnames(metadata), collapse = ", ")))
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
      warning('The following ', length(unk), ' ', vCol, ' values were not found in the DB: ', paste0(unk, collapse = ','))
      warning(paste0("Run tcrClustR:::.PullTcrdist3Db(organism = '", organism, "', outputFilePath = '...') to obtain the list of known segments."))
      metadata[[validCol]][unknown_v_segments] <- FALSE
    }

    toTest <- as.character(metadata[[jCol]])
    unknown_j_segments <- metadata[[validCol]] & !is.na(toTest) & !(toTest %in% gene_segments_in_db)
    if (sum(unknown_j_segments) > 0) {
      unk <- sort(unique(toTest[unknown_j_segments]))
      warning('The following ', length(unk), ' ', jCol, ' values were not found in the DB: ', paste0(unk, collapse = ','))
      warning(paste0("Run tcrClustR:::.PullTcrdist3Db(organism = '", organism, "', outputFilePath = '...') to obtain the list of known segments."))
      metadata[[validCol]][unknown_j_segments] <- FALSE
    }

    if (verbose) {
      message(paste0("Finished metadata cleaning for ", chain, ". Initial rows: ", nrow(metadata), ", remaining after filters: ", sum(metadata[[validCol]])))
    }
  }

  return(metadata)
}