#' @title JoinClusteringResults
#' @description Join TCR clustering results from parquet files back to a Seurat object.
#' This function reads parquet files produced by RunTcrClustering() and joins cluster
#' assignments to the Seurat object metadata by matching V gene, J gene, and CDR3 amino acid
#' columns. Handles both single-chain and paired-chain parquet schemas.
#'
#' @param seuratObj Seurat object to add clustering results to (typically the original object with cell-level data, but not required to be).
#' @param parquetFiles Character vector of paths to parquet files containing clustering results.
#'   If NULL, will auto-detect files in `parquetDir` matching the `filePattern`.
#' @param parquetDir Directory containing parquet files. Default = "./tcrclustering_output/".
#'   Only used if `parquetFiles` is NULL.
#' @param filePattern Regular expression pattern to match parquet files. Default = "\\.parquet$".
#'   Only used if `parquetFiles` is NULL.
#' @param metadataColumnPrefix Prefix for new metadata columns added to the Seurat object.
#'   Default = "TcrFamily". Final column names will be like "TcrFamily_TRA", "TcrFamily_TRB", etc.
#' @param stripAlleles Logical. If TRUE (default), strip allele suffixes (e.g., "*01") from V/J gene names
#'   in the Seurat metadata before matching to parquet data (which should already be stripped by RunTcrClustering).
#' @param overwriteExisting Logical. If TRUE, overwrite existing cluster columns. If FALSE (default), skip
#'   files whose cluster column already exists in metadata.
#' @param verbose Logical. Print progress messages. Default = TRUE.
#'
#' @return Updated Seurat object with new metadata columns containing cluster assignments.
#'   Column names follow the pattern: `{metadataColumnPrefix}_{chain}` (e.g., "TcrFamily_TRA", "TcrFamily_TRB").
#'   Cells/clones without cluster assignments will have NA values.
#'
#' @details
#' **Matching Logic:**
#' - **Single-chain files** (schema: chain, v_gene, j_gene, CDR3, Cluster, ...):
#'   Matches on V gene + J gene + CDR3 for the specified chain (e.g., TRA_V + TRA_J + TRA).
#'
#' - **Paired-chain files** (schema: chain_1, v_gene_1, j_gene_1, CDR3_1, chain_2, v_gene_2, j_gene_2, CDR3_2, Cluster, ...):
#'   Requires a simultaneous match on BOTH chains (AND join). A cell is assigned the cluster only when
#'   the combination of Chain 1 (v_gene_1 + j_gene_1 + CDR3_1) AND Chain 2 (v_gene_2 + j_gene_2 + CDR3_2)
#'   matches a row in the parquet data. Cells with only one chain detected will not receive a paired-chain
#'   assignment (remain NA in the paired column). For single-chain matching, use the single-chain parquet outputs.
#'
#' **Allele Handling:**
#' If `stripAlleles = TRUE`, the function removes allele suffixes (e.g., "*01") from V/J gene names
#' in the Seurat metadata to match the stripped values in parquet files (assuming RunTcrClustering
#' was run with stripAlleles = TRUE, which is the default).
#'
#' **Column Naming:**
#' New columns are named based on the detected chain(s) in the parquet file:
#' - Single-chain: `{metadataColumnPrefix}_TRA`, `{metadataColumnPrefix}_TRB`, etc.
#' - Paired-chain: `{metadataColumnPrefix}_TRA_TRB`, `{metadataColumnPrefix}_TRACDR3_TRB`, etc.
#'
#' @examples
#' \dontrun{
#'   # Basic usage: join clustering results after running RunTcrClustering()
#'   clustering_output <- RunTcrClustering(seuratObj_TCR, outputPrefix = "my_clusters")
#'   seuratObj <- JoinClusteringResults(seuratObj,
#'                                       parquetFiles = clustering_output$parquet_files)
#'
#'   # Auto-detect parquet files in a directory
#'   seuratObj <- JoinClusteringResults(seuratObj,
#'                                       parquetDir = "./tcrclustering_output/",
#'                                       metadataColumnPrefix = "TcrFamily")
#'
#'   # Custom column prefix
#'   seuratObj <- JoinClusteringResults(seuratObj,
#'                                       parquetFiles = clustering_output$parquet_files,
#'                                       metadataColumnPrefix = "MyCluster")
#' }
#'
#' @export
JoinClusteringResults <- function(seuratObj,
                                   parquetFiles = NULL,
                                   parquetDir = "./tcrclustering_output/",
                                   filePattern = "\\.parquet$",
                                   metadataColumnPrefix = "TcrFamily",
                                   stripAlleles = TRUE,
                                   overwriteExisting = FALSE,
                                   verbose = TRUE) {

  # Validate inputs
  if (!inherits(seuratObj, "Seurat")) {
    stop("seuratObj must be a Seurat object")
  }

  # Auto-detect parquet files if not provided
  if (is.null(parquetFiles)) {
    if (!dir.exists(parquetDir)) {
      stop("parquetDir '", parquetDir, "' does not exist. Please provide valid parquetFiles or parquetDir.")
    }
    parquetFiles <- list.files(parquetDir, pattern = filePattern, full.names = TRUE)
    if (length(parquetFiles) == 0) {
      stop("No parquet files found in '", parquetDir, "' matching pattern '", filePattern, "'")
    }
    if (verbose) message("Auto-detected ", length(parquetFiles), " parquet file(s) in ", parquetDir)
  }

  # Validate parquet files exist
  missing_files <- parquetFiles[!file.exists(parquetFiles)]
  if (length(missing_files) > 0) {
    stop("The following parquet files do not exist:\n  ", paste(missing_files, collapse = "\n  "))
  }

  metadata <- seuratObj@meta.data

  # Strip alleles from metadata if requested (to match parquet data)
  if (stripAlleles) {
    if (verbose) message("Stripping allele suffixes from Seurat metadata V/J genes...")
    for (col in colnames(metadata)) {
      # Strip alleles from V and J gene columns (TRA_V, TRA_J, TRB_V, TRB_J, etc.)
      if (grepl("_(V|J)$", col)) {
        metadata[[col]] <- sapply(metadata[[col]], .strip_allele_suffix)
      }
    }
  }

  # Process each parquet file
  for (parquet_file in parquetFiles) {
    if (verbose) message("\nProcessing: ", basename(parquet_file))

    # Read parquet file
    cluster_data <- tryCatch({
      arrow::read_parquet(parquet_file)
    }, error = function(e) {
      warning("Failed to read parquet file '", parquet_file, "': ", e$message)
      return(NULL)
    })

    if (is.null(cluster_data)) next

    # Detect schema type (single-chain vs paired-chain)
    is_paired <- all(c("chain_1", "v_gene_1", "j_gene_1", "CDR3_1") %in% colnames(cluster_data))
    is_single <- all(c("chain", "v_gene", "j_gene", "CDR3") %in% colnames(cluster_data))

    if (!is_paired && !is_single) {
      warning("Parquet file '", basename(parquet_file), "' has unexpected schema. Skipping.")
      next
    }

    # Determine column name for cluster assignments
    if (is_single) {
      # Extract chain name from parquet data
      chain_name <- unique(cluster_data$chain)[1]  # Should be uniform (e.g., "TRA", "TRB")
      if (is.na(chain_name) || chain_name == "unknown") {
        warning("Could not determine chain name from parquet file '", basename(parquet_file), "'. Skipping.")
        next
      }
      new_col_name <- paste0(metadataColumnPrefix, "_", chain_name)

    } else {  # is_paired
      # Extract both chain names
      chain_1 <- unique(cluster_data$chain_1)[1]
      chain_2 <- unique(cluster_data$chain_2)[1]

      # Handle NA chains (only one chain present in some rows)
      chain_1 <- if (is.na(chain_1)) "NA" else chain_1
      chain_2 <- if (is.na(chain_2)) "NA" else chain_2

      # Create combined name (e.g., "TRA_TRB")
      chain_combo <- paste(sort(c(chain_1, chain_2)), collapse = "_")
      chain_combo <- gsub("_NA|NA_", "", chain_combo)  # Remove NA placeholders

      new_col_name <- paste0(metadataColumnPrefix, "_", chain_combo)
    }

    # Check if column already exists
    if (new_col_name %in% colnames(metadata) && !overwriteExisting) {
      if (verbose) message("  Column '", new_col_name, "' already exists. Skipping (set overwriteExisting=TRUE to replace).")
      next
    }

    # Join cluster assignments to metadata
    if (is_single) {
      matches <- .join_single_chain(metadata, cluster_data, chain_name, verbose)
    } else {
      matches <- .join_paired_chain(metadata, cluster_data, verbose)
    }

    # Add cluster column to metadata
    metadata[[new_col_name]] <- matches
    n_matched <- sum(!is.na(matches))
    if (verbose) message("  Joined ", n_matched, " cluster assignments to column '", new_col_name, "'")
  }

  # Update Seurat object metadata
  seuratObj@meta.data <- metadata

  if (verbose) message("\nSuccessfully joined clustering results to Seurat object.")
  return(seuratObj)
}


# Internal helper: Join single-chain parquet data to metadata
# Returns a vector of cluster assignments (same length as nrow(metadata))
.join_single_chain <- function(metadata, cluster_data, chain_name, verbose = FALSE) {

  # Identify metadata columns for this chain
  v_col <- paste0(chain_name, "_V")
  j_col <- paste0(chain_name, "_J")
  cdr3_col <- chain_name

  # Check if required columns exist in metadata
  missing_cols <- c(v_col, j_col, cdr3_col)[!c(v_col, j_col, cdr3_col) %in% colnames(metadata)]
  if (length(missing_cols) > 0) {
    warning("  Missing metadata columns for chain ", chain_name, ": ", paste(missing_cols, collapse = ", "), ". Returning NA.")
    return(rep(NA_character_, nrow(metadata)))
  }

  # Create join keys
  # Metadata join key: paste V_gene, J_gene, CDR3
  metadata_key <- paste(metadata[[v_col]], metadata[[j_col]], metadata[[cdr3_col]], sep = "|||")

  # Parquet join key: paste v_gene, j_gene, CDR3
  parquet_key <- paste(cluster_data$v_gene, cluster_data$j_gene, cluster_data$CDR3, sep = "|||")

  # Create lookup table (parquet_key -> Cluster)
  lookup <- setNames(cluster_data$Cluster, parquet_key)

  # Match metadata keys to parquet keys
  cluster_assignments <- lookup[metadata_key]
  names(cluster_assignments) <- NULL

  return(cluster_assignments)
}


# Internal helper: Join paired-chain parquet data to metadata
# Requires BOTH chains to match simultaneously (AND join)
# Returns a vector of cluster assignments (same length as nrow(metadata))
.join_paired_chain <- function(metadata, cluster_data, verbose = FALSE) {

  # Extract chain names from parquet data (should be uniform per file)
  chain_1 <- unique(cluster_data$chain_1)[1]
  chain_2 <- unique(cluster_data$chain_2)[1]

  # Validate chain names
  if (is.na(chain_1) || is.na(chain_2) || chain_1 == "" || chain_2 == "") {
    if (verbose) message("  Paired parquet has missing chain names; cannot perform AND join. Returning NA.")
    return(rep(NA_character_, nrow(metadata)))
  }

  # Identify required metadata columns for both chains
  v_col_1 <- paste0(chain_1, "_V"); j_col_1 <- paste0(chain_1, "_J"); cdr3_col_1 <- chain_1
  v_col_2 <- paste0(chain_2, "_V"); j_col_2 <- paste0(chain_2, "_J"); cdr3_col_2 <- chain_2

  required_cols <- c(v_col_1, j_col_1, cdr3_col_1, v_col_2, j_col_2, cdr3_col_2)
  missing_cols <- required_cols[!required_cols %in% colnames(metadata)]
  if (length(missing_cols) > 0) {
    warning("  Missing metadata columns for paired join (", chain_1, "+", chain_2, "): ", paste(missing_cols, collapse = ", "), ". Returning NA.")
    return(rep(NA_character_, nrow(metadata)))
  }

  # Filter parquet rows to those with complete info for BOTH chains
  complete_rows <- !is.na(cluster_data$v_gene_1) & !is.na(cluster_data$j_gene_1) & !is.na(cluster_data$CDR3_1) &
                   !is.na(cluster_data$v_gene_2) & !is.na(cluster_data$j_gene_2) & !is.na(cluster_data$CDR3_2)
  if (!any(complete_rows)) {
    if (verbose) message("  No complete paired entries in parquet for ", chain_1, "+", chain_2, ". Returning NA.")
    return(rep(NA_character_, nrow(metadata)))
  }
  cluster_data_complete <- cluster_data[complete_rows, , drop = FALSE]

  # Build composite join keys: [v1|j1|cdr3_1||v2|j2|cdr3_2]
  # Use a delimiter unlikely to appear in gene names/seqs
  make_key <- function(v, j, c) paste(v, j, c, sep = "|||")

  metadata_key <- paste(
    make_key(metadata[[v_col_1]], metadata[[j_col_1]], metadata[[cdr3_col_1]]),
    make_key(metadata[[v_col_2]], metadata[[j_col_2]], metadata[[cdr3_col_2]]),
    sep = "||||"
  )

  parquet_key <- paste(
    make_key(cluster_data_complete$v_gene_1, cluster_data_complete$j_gene_1, cluster_data_complete$CDR3_1),
    make_key(cluster_data_complete$v_gene_2, cluster_data_complete$j_gene_2, cluster_data_complete$CDR3_2),
    sep = "||||"
  )

  # Create lookup: composite key -> Cluster
  # If duplicates exist with conflicting clusters, keep the first and warn
  if (any(duplicated(parquet_key))) {
    dup_keys <- unique(parquet_key[duplicated(parquet_key)])
    # Check for conflicting cluster labels among duplicates
    for (dk in dup_keys) {
      labs <- unique(as.character(cluster_data_complete$Cluster[parquet_key == dk]))
      if (length(labs) > 1 && verbose) {
        warning("  Paired parquet has conflicting clusters for the same composite key; using first encountered.")
        break
      }
    }
  }

  lookup <- stats::setNames(as.character(cluster_data_complete$Cluster), parquet_key)
  matches <- lookup[metadata_key]
  names(matches) <- NULL
  return(matches)
}


# Helper function to strip allele suffixes from gene names (reused from RunTcrClustering.R)
# This is duplicated here to keep JoinClusteringResults.R self-contained
.strip_allele_suffix <- function(gene_name) {
  if (is.na(gene_name) || length(gene_name) == 0) return(gene_name)
  gsub("\\*[0-9]+$", "", as.character(gene_name))
}
