#' @title RunTcrClustering
#' @description Cluster TCR distance assays in a single Seurat/tcr object and write compact tables to disk.
#' This function expects the user to have already subsetted their data (no grouping column).
#' It supports single-chain assays and multi-chain (paired) assays. For single-chain assays the
#' output table has columns: chain, v_gene, j_gene, CDR3, Cluster, Cluster_Size_Threshold, Clustering_Method, ...
#' For paired assays (two chains) the output table is written with columns: chain_1, v_gene_1, j_gene_1, CDR3_1, chain_2, v_gene_2, j_gene_2, CDR3_2, Cluster, ...
#'
#' @param seuratObj_TCR Seurat object produced by RunTcrdist3 (contains distance assays and metadata). If NULL, provide `metadata` and `distance_matrices` instead.
#' @param metadata Optional metadata data.frame (row.names must match clone IDs). If NULL the function will use `seuratObj_TCR@meta.data`.
#' @param assaysToCluster Character vector of assay names to cluster. If NULL, all assays in `seuratObj_TCR` will be processed.
#' @param clusteringMethod "DIANA" or "Leiden". Default = "DIANA". If "Leiden", ClusterTcrs() is used.
#' @param dianaHeight DIANA cut height (used when clusteringMethod == "DIANA"). Default = 20.
#' @param leidenResolution Resolution for Leiden (passed to ClusterTcrs if used). Default = 0.3.
#' @param clusterSizeThreshold Minimum cluster size to keep (clusters below set to 0). Default = 2.
#' @param outputDir Directory to write parquet files. Default = "./tcrclustering_output/".
#' @param outputPrefix Optional prefix for output filenames. Default = "clustering".
#' @param maxChainsToRecord Max number of chains to record for paired output. Default = 2.
#' @param stripAlleles Logical. If TRUE (default), strip allele suffixes (e.g., "*01") from V/J gene names.
#' @param filterInvalidClones Logical. If TRUE (default), filter clones with concatenated values (commas/semicolons) or NA values in chain-specific columns before clustering.
#' @param verbose Logical. Print progress. Default = TRUE.
#' @param ... Additional arguments forwarded to ClusterTcrs when Leiden is used.
#' @return A list with `parquet_files` (vector of files written) and `summary` (data.frame summary per assay)
#' @export
RunTcrClustering <- function(seuratObj_TCR = NULL,
                             metadata = NULL,
                             assaysToCluster = NULL,
                             clusteringMethod = "DIANA",
                             dianaHeight = 20,
                             leidenResolution = 0.3,
                             clusterSizeThreshold = 2,
                             outputDir = "./tcrclustering_output/",
                             outputPrefix = "clustering",
                             maxChainsToRecord = 2,
                             stripAlleles = TRUE,
                             filterInvalidClones = TRUE,
                             verbose = TRUE,
                             ...) {
  if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)

  if (is.null(metadata)) {
    if (is.null(seuratObj_TCR)) stop("Please provide either seuratObj_TCR or metadata")
    metadata <- seuratObj_TCR@meta.data
  }

  # determine assays to process
  if (!is.null(seuratObj_TCR)) {
    available_assays <- SeuratObject::Assays(seuratObj_TCR)
  } else {
    available_assays <- character(0)
  }

  if (is.null(assaysToCluster)) {
    assays <- available_assays
  } else {
    assays <- intersect(assaysToCluster, available_assays)
    if (length(assays) == 0) stop("No assays found to cluster. Provide valid assay names or a seuratObj_TCR with assays.")
  }

  parquet_files <- character()
  summary_rows <- list()

  for (assay in assays) {
    if (verbose) message("Processing assay: ", assay)

    # make sure the assay data layer is accessible
    if (methods::is(seuratObj_TCR[[assay]], "Assay5")) {
      seuratObj_TCR[[assay]] <- SeuratObject::JoinLayers(seuratObj_TCR[[assay]])
    }

    dist_matrix <- tryCatch({
      # For Assay5, try counts layer first (where tcrdist3 stores distance matrices)
      if (methods::is(seuratObj_TCR[[assay]], "Assay5")) {
        as.matrix(SeuratObject::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))
      } else {
        as.matrix(SeuratObject::GetAssayData(seuratObj_TCR, assay = assay, layer = "data"))
      }
    }, error = function(e) {
      warning("Could not extract distance matrix for assay ", assay, ": ", e$message); return(NULL)
    })
    if (is.null(dist_matrix)) next
    
    # check if matrix is too small or empty
    if (nrow(dist_matrix) < 2 || ncol(dist_matrix) < 2) {
      if (verbose) message("Skipping assay ", assay, ": distance matrix too small (", nrow(dist_matrix), "x", ncol(dist_matrix), ")")
      next
    }

    # clustering
    if (clusteringMethod == "DIANA") {
      clustering_result <- .DianaClustering(dist_matrix, cutHeight = dianaHeight, verbose = verbose)
      cluster_assignments <- clustering_result$clustering
      clustering_params <- list(cutHeight = dianaHeight)
    } else if (clusteringMethod == "Leiden") {
      if (verbose) message("Running Leiden via ClusterTcrs on assay: ", assay)
      clustered_obj <- ClusterTcrs(seuratObj_TCR = seuratObj_TCR,
                                   resolutionParameters = leidenResolution,
                                   computeMultiChain = TRUE,
                                   verbose = verbose,
                                   ...)
      # try to find column
      cluster_col <- paste0("TcrClustR_", ifelse(grepl("_cdr3", assay), gsub("_cdr3", "CDR3", assay), assay), "_", leidenResolution)
      if (cluster_col %in% colnames(clustered_obj[[assay]]@meta.data)) {
        cluster_assignments <- clustered_obj[[assay]]@meta.data[[cluster_col]]
        names(cluster_assignments) <- rownames(clustered_obj[[assay]]@meta.data)
      } else {
        warning("Leiden clustering column not found for assay ", assay, ". Skipping.")
        next
      }
      clustering_params <- list(resolution = leidenResolution)
    } else {
      stop("Unsupported clusteringMethod: ", clusteringMethod)
    }

    # apply cluster size threshold
    filtered_assignments <- .ThresholdClustersBySize(cluster_assignments, minClusterSize = clusterSizeThreshold, verbose = FALSE)

    # detect chains present in assay name (reuse pattern from Clustering.R)
    assay_normalized <- gsub("CDR3", "_cdr3", assay)
    assay_parts <- strsplit(assay_normalized, "_")[[1]]
    all_possible_chains <- c("TRA", "TRB", "TRG", "TRD")
    chains_in_assay <- all_possible_chains[sapply(all_possible_chains, function(chain) any(grepl(chain, assay_parts, fixed = TRUE)))]
    
    # skip if no chains detected
    if (length(chains_in_assay) == 0) {
      if (verbose) message("Skipping assay ", assay, ": could not detect chain type")
      next
    }

    # Filter clones if requested
    if (filterInvalidClones) {
      valid_clones <- .filter_clones_for_assay(metadata, chains_in_assay, verbose = verbose)
      
      # Filter distance matrix and cluster assignments
      if (sum(valid_clones) < nrow(dist_matrix)) {
        # Get IDs of valid clones (handle hyphen/underscore conversion)
        valid_ids <- rownames(metadata)[valid_clones]
        valid_ids_hyphen <- gsub("_", "-", valid_ids, fixed = TRUE)
        
        # Find matching rows in distance matrix
        valid_idx <- which(rownames(dist_matrix) %in% valid_ids_hyphen | rownames(dist_matrix) %in% valid_ids)
        
        if (length(valid_idx) < 2) {
          if (verbose) message("Skipping assay ", assay, ": too few valid clones after filtering (", length(valid_idx), ")")
          next
        }
        
        # Subset distance matrix
        dist_matrix <- dist_matrix[valid_idx, valid_idx, drop = FALSE]
        
        # Re-run clustering on filtered matrix
        if (clusteringMethod == "DIANA") {
          clustering_result <- .DianaClustering(dist_matrix, cutHeight = dianaHeight, verbose = verbose)
          filtered_assignments <- .ThresholdClustersBySize(clustering_result$clustering, minClusterSize = clusterSizeThreshold, verbose = FALSE)
        } else {
          # For Leiden, need to re-run on subset - more complex, skip for now or filter assignments
          # Filter existing assignments
          filtered_assignments <- filtered_assignments[names(filtered_assignments) %in% rownames(dist_matrix)]
        }
        
        if (verbose) message("Filtered from ", nrow(metadata), " to ", length(filtered_assignments), " clones")
      }
    }

    # prepare table depending on how many chains
    n_rows <- length(filtered_assignments)
    ids <- names(filtered_assignments)
    if (is.null(ids) || any(ids == "")) {
      # try to fallback to rownames of distance matrix
      if (!is.null(rownames(dist_matrix))) ids <- rownames(dist_matrix)
      names(filtered_assignments) <- ids
    }

    if (length(chains_in_assay) <= 1) {
      # single-chain schema: chain, v_gene, j_gene, CDR3
      chain_name <- ifelse(length(chains_in_assay) == 1, chains_in_assay, NA_character_)
      
      # Map to tcrdist3 column names (v_a_gene, v_b_gene, etc.)
      chain_letter <- tolower(substr(chain_name, 3, 3))  # "TRB" -> "b", "TRA" -> "a", etc.
      v_col_tcrdist <- paste0("v_", chain_letter, "_gene")
      j_col_tcrdist <- paste0("j_", chain_letter, "_gene")  
      cdr3_col_tcrdist <- paste0("cdr3_", chain_letter, "_aa")
      
      # Try tcrdist3 column names first, fall back to original names
      v_col <- if (v_col_tcrdist %in% colnames(metadata)) v_col_tcrdist else if (!is.na(chain_name)) paste0(chain_name, "_V") else NA_character_
      j_col <- if (j_col_tcrdist %in% colnames(metadata)) j_col_tcrdist else if (!is.na(chain_name)) paste0(chain_name, "_J") else NA_character_
      cdr3_col <- if (cdr3_col_tcrdist %in% colnames(metadata)) cdr3_col_tcrdist else if (!is.na(chain_name)) chain_name else NA_character_

      # Fix rowname matching - distance matrix has hyphens, metadata has underscores
      ids_fixed <- gsub("-", "_", ids, fixed = TRUE)

      # Extract metadata values
      v_gene_vals <- if (!is.na(v_col) && v_col %in% colnames(metadata)) metadata[[v_col]][match(ids_fixed, rownames(metadata))] else rep(NA_character_, n_rows)
      j_gene_vals <- if (!is.na(j_col) && j_col %in% colnames(metadata)) metadata[[j_col]][match(ids_fixed, rownames(metadata))] else rep(NA_character_, n_rows)
      
      # Strip alleles if requested
      if (stripAlleles) {
        v_gene_vals <- sapply(v_gene_vals, .strip_allele_suffix)
        j_gene_vals <- sapply(j_gene_vals, .strip_allele_suffix)
      }

      out_df <- data.frame(
        chain = rep(ifelse(is.na(chain_name), "unknown", chain_name), n_rows),
        v_gene = v_gene_vals,
        j_gene = j_gene_vals,
        CDR3 = if (!is.na(cdr3_col) && cdr3_col %in% colnames(metadata)) metadata[[cdr3_col]][match(ids_fixed, rownames(metadata))] else rep(NA_character_, n_rows),
        Cluster = as.character(filtered_assignments),
        Cluster_Size_Threshold = rep(clusterSizeThreshold, n_rows),
        Clustering_Method = rep(clusteringMethod, n_rows),
        stringsAsFactors = FALSE
      )

      # add clustering params columns if any
      if (length(clustering_params) > 0) {
        for (p in names(clustering_params)) {
          out_df[[paste0("param_", p)]] <- rep(as.character(clustering_params[[p]]), n_rows)
        }
      }

      # write parquet
      out_file <- file.path(outputDir, paste0(outputPrefix, "_", gsub("[^[:alnum:]_]", "_", assay), "_single.parquet"))
      arrow::write_parquet(out_df, out_file)
      parquet_files <- c(parquet_files, out_file)
      summary_rows[[assay]] <- data.frame(assay = assay, chains_used = paste(chains_in_assay, collapse = "+"), n_clones = n_rows, n_clusters = length(unique(out_df$Cluster[out_df$Cluster != "0"])), stringsAsFactors = FALSE)

    } else {
      # multi-chain: create chain_1, v_gene_1, j_gene_1, CDR3_1, chain_2, ... up to maxChainsToRecord
      n_chains_to_record <- min(length(chains_in_assay), maxChainsToRecord)
      selected_chains <- sort(chains_in_assay)[seq_len(n_chains_to_record)]
      
      # Fix rowname matching
      ids_fixed <- gsub("-", "_", ids, fixed = TRUE)
      
      out_df <- data.frame(stringsAsFactors = FALSE)
      for (i in seq_along(selected_chains)) {
        ch <- selected_chains[i]
        chain_letter <- tolower(substr(ch, 3, 3))
        
        # Try tcrdist3 column names first, fall back to original
        v_col_tcrdist <- paste0("v_", chain_letter, "_gene")
        j_col_tcrdist <- paste0("j_", chain_letter, "_gene")
        cdr3_col_tcrdist <- paste0("cdr3_", chain_letter, "_aa")
        
        v_col <- if (v_col_tcrdist %in% colnames(metadata)) v_col_tcrdist else paste0(ch, "_V")
        j_col <- if (j_col_tcrdist %in% colnames(metadata)) j_col_tcrdist else paste0(ch, "_J")
        cdr3_col <- if (cdr3_col_tcrdist %in% colnames(metadata)) cdr3_col_tcrdist else ch
        
        # Extract values
        v_gene_vals <- if (v_col %in% colnames(metadata)) metadata[[v_col]][match(ids_fixed, rownames(metadata))] else rep(NA_character_, n_rows)
        j_gene_vals <- if (j_col %in% colnames(metadata)) metadata[[j_col]][match(ids_fixed, rownames(metadata))] else rep(NA_character_, n_rows)
        
        # Strip alleles if requested
        if (stripAlleles) {
          v_gene_vals <- sapply(v_gene_vals, .strip_allele_suffix)
          j_gene_vals <- sapply(j_gene_vals, .strip_allele_suffix)
        }
        
        out_df[[paste0("chain_", i)]] <- rep(ch, n_rows)
        out_df[[paste0("v_gene_", i)]] <- v_gene_vals
        out_df[[paste0("j_gene_", i)]] <- j_gene_vals
        out_df[[paste0("CDR3_", i)]] <- if (cdr3_col %in% colnames(metadata)) metadata[[cdr3_col]][match(ids_fixed, rownames(metadata))] else rep(NA_character_, n_rows)
      }

      # if fewer than maxChainsToRecord, pad with NA columns
      if (length(selected_chains) < maxChainsToRecord) {
        for (i in (length(selected_chains)+1):maxChainsToRecord) {
          out_df[[paste0("chain_", i)]] <- rep(NA_character_, n_rows)
          out_df[[paste0("v_gene_", i)]] <- rep(NA_character_, n_rows)
          out_df[[paste0("j_gene_", i)]] <- rep(NA_character_, n_rows)
          out_df[[paste0("CDR3_", i)]] <- rep(NA_character_, n_rows)
        }
      }

      out_df$Cluster <- as.character(filtered_assignments)
      out_df$Cluster_Size_Threshold <- rep(clusterSizeThreshold, n_rows)
      out_df$Clustering_Method <- rep(clusteringMethod, n_rows)

      if (length(clustering_params) > 0) {
        for (p in names(clustering_params)) {
          out_df[[paste0("param_", p)]] <- rep(as.character(clustering_params[[p]]), n_rows)
        }
      }

      out_file <- file.path(outputDir, paste0(outputPrefix, "_", gsub("[^[:alnum:]_]", "_", assay), "_paired.parquet"))
      arrow::write_parquet(out_df, out_file)
      parquet_files <- c(parquet_files, out_file)
      summary_rows[[assay]] <- data.frame(assay = assay, chains_used = paste(chains_in_assay, collapse = "+"), n_clones = n_rows, n_clusters = length(unique(out_df$Cluster[out_df$Cluster != "0"])), stringsAsFactors = FALSE)
    }

    if (verbose) message("Wrote output for assay: ", assay)
  }

  summary_df <- do.call(rbind, summary_rows)
  return(list(parquet_files = parquet_files, summary = summary_df))
}

# Helper function to check if a value is concatenated (contains comma or semicolon)
.is_concatenated <- function(x) {
  if (is.na(x) || length(x) == 0) return(FALSE)
  grepl("[,;]", as.character(x))
}

# Helper function to filter clones based on chain requirements
# For single-chain: filter if chain columns have NA or concatenated values
# For multi-chain: filter if ALL chains in the assay have NA or concatenated values
#   (i.e., a clone is valid if at least one chain is complete)
.filter_clones_for_assay <- function(metadata, chains_in_assay, verbose = FALSE) {
  if (length(chains_in_assay) == 0) {
    if (verbose) message("No chains detected, returning all clones")
    return(rep(TRUE, nrow(metadata)))
  }
  
  # Build column names to check (both original and tcrdist3 formats)
  cols_to_check <- character(0)
  for (chain in chains_in_assay) {
    # Original format: TRA_V, TRA_J, TRA
    cols_to_check <- c(cols_to_check, paste0(chain, "_V"), paste0(chain, "_J"), chain)
    
    # tcrdist3 format: v_a_gene, j_a_gene, cdr3_a_aa
    chain_letter <- tolower(substr(chain, 3, 3))
    cols_to_check <- c(cols_to_check, 
                      paste0("v_", chain_letter, "_gene"),
                      paste0("j_", chain_letter, "_gene"),
                      paste0("cdr3_", chain_letter, "_aa"))
  }
  
  # Keep only columns that exist in metadata
  cols_to_check <- intersect(cols_to_check, colnames(metadata))
  
  if (length(cols_to_check) == 0) {
    if (verbose) message("No relevant chain columns found in metadata")
    return(rep(TRUE, nrow(metadata)))
  }
  
  if (length(chains_in_assay) == 1) {
    # SINGLE-CHAIN: Clone is valid only if ALL chain-specific columns are non-NA and non-concatenated
    keep <- rep(TRUE, nrow(metadata))
    for (col in cols_to_check) {
      keep <- keep & !is.na(metadata[[col]]) & !sapply(metadata[[col]], .is_concatenated)
    }
    
    if (verbose) {
      n_filtered <- sum(!keep)
      if (n_filtered > 0) {
        message("Single-chain filtering: removed ", n_filtered, " clones with NA or concatenated values in ", paste(chains_in_assay, collapse = "+"))
      }
    }
    
  } else {
    # MULTI-CHAIN: Clone is valid if AT LEAST ONE chain is complete (all its columns non-NA and non-concatenated)
    # Group columns by chain
    chain_valid <- matrix(TRUE, nrow = nrow(metadata), ncol = length(chains_in_assay))
    
    for (i in seq_along(chains_in_assay)) {
      chain <- chains_in_assay[i]
      chain_letter <- tolower(substr(chain, 3, 3))
      
      # Find this chain's columns
      chain_cols <- c(
        paste0(chain, "_V"), paste0(chain, "_J"), chain,  # Original
        paste0("v_", chain_letter, "_gene"),              # tcrdist3
        paste0("j_", chain_letter, "_gene"),
        paste0("cdr3_", chain_letter, "_aa")
      )
      chain_cols <- intersect(chain_cols, colnames(metadata))
      
      if (length(chain_cols) > 0) {
        # Chain is valid if all its columns are non-NA and non-concatenated
        for (col in chain_cols) {
          chain_valid[, i] <- chain_valid[, i] & !is.na(metadata[[col]]) & !sapply(metadata[[col]], .is_concatenated)
        }
      }
    }
    
    # Keep clone if ANY chain is valid
    keep <- rowSums(chain_valid) > 0
    
    if (verbose) {
      n_filtered <- sum(!keep)
      if (n_filtered > 0) {
        message("Multi-chain filtering: removed ", n_filtered, " clones where ALL chains had NA or concatenated values")
      }
    }
  }
  
  return(keep)
}

# Helper function to strip allele suffixes from gene names
.strip_allele_suffix <- function(gene_name) {
  if (is.na(gene_name) || length(gene_name) == 0) return(gene_name)
  gsub("\\*[0-9]+$", "", as.character(gene_name))
}
