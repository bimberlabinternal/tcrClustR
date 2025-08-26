
utils::globalVariables(
  names = c('matrix_rowname.x', 'matrix_rowname.y', '.data', '.', 'distanceMatrix', 'combined_matrix',
            'key', 'seed', 'seuratObj'),
  package = 'tcrClustR',
  add = TRUE
)

#TODO: summarize data lossy-ness vignette.

#' @title Cluster TCRs using TcrClustR
#' @description This function clusters TCRs in a Seurat object using TcrClustR.
#' It performs PCA or kernel PCA, computes distance matrices, and applies clustering algorithms.
#' The clustering results are stored in the metadata of the returned Seurat objects.
#' Users can join clustering results back to their original object using clonotypic metadata.
#' @param seuratObj_TCR Seurat object with TCR distance matrices.
#' @param resolutionParameter Resolution parameter for clustering. Default is 0.1.
#' @param pcaComponents Number of components for PCA or kernel PCA. Default is 50.
#' @param kpcaKernel Kernel type for kernel PCA. Default is "rbfdot". Ignored if usePCA is TRUE.
#' @param usePCA Boolean indicating whether to use standard PCA instead of kernel PCA. Default is TRUE.
#' @param proportionOfGraphAsNeighbors Proportion of the graph to consider as neighbors. Default is 0.1.
#' @param jaccardIndexThreshold Jaccard index threshold for pruning edges. Default is 0.1.
#' @param seed Random seed for reproducibility. Default is 1234.
#' @param computeMultiChain Boolean indicating whether to compute multi-chain graphs. Default is TRUE.
#' @param verbose Boolean indicating whether to print verbose debugging output. Default is FALSE.
#' @return A list containing single-chain and multi-chain Seurat objects with clustering results in their metadata.
#'   Users should perform clonotypic joins to transfer clustering information to their original Seurat object.
#' @export


ClusterTcrs <- function(seuratObj_TCR = NULL,
                        resolutionParameter = 0.1,
                        pcaComponents = 50,
                        kpcaKernel = "rbfdot",
                        usePCA = TRUE,
                        proportionOfGraphAsNeighbors = 0.1,
                        jaccardIndexThreshold = 0.1,
                        seed = 1234,
                        computeMultiChain = T,
                        verbose = FALSE) {

  #perform leiden clustering on the single chain distance matrices, create the multichain distance matrices, and cluster them.
  clusteredSeuratObjects <- .DistanceMatrixToClusteredGraphs(seuratObj_TCR = seuratObj_TCR,
                                                             pcaComponents = pcaComponents,
                                                             kpcaKernel = kpcaKernel,
                                                             usePCA = usePCA,
                                                             proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                                             jaccardIndexThreshold = jaccardIndexThreshold,
                                                             seed = seed,
                                                             computeMultiChain = computeMultiChain,
                                                             verbose = verbose)

  #filter out NULL objects to prevent assay access errors when iterating
  clusteredSeuratObjects <- clusteredSeuratObjects[!sapply(clusteredSeuratObjects, is.null)]

  #return the clustered TCR objects
  #users can join clustering results to their original object using clonotypic metadata
  #see the vignette for examples of how to do this join
  return(clusteredSeuratObjects)
}



#' @title .DistanceMatrixToClusteredGraphs
#' @description
#' This function takes a Seurat object with TCR distance matrices and computes clustered graphs for each chain.
#' @param seuratObj_TCR Seurat object containing TCR distance matrices.
#' @param pcaComponents Number of components for PCA or kernel PCA. Default is 50.
#' @param kpcaKernel Kernel type for kernel PCA. Default is "rbfdot". Ignored if usePCA is TRUE.
#' @param usePCA Boolean indicating whether to use standard PCA instead of kernel PCA. Default is TRUE.
#' @param partitionType Type of partitioning algorithm to use. Default is "CPMVertexPartition".
#' @param proportionOfGraphAsNeighbors Proportion of the graph to consider as neighbors. Default is 0.1.
#' @param jaccardIndexThreshold Jaccard index threshold for pruning edges. Default is 0.1.
#' @param resolutions Vector of resolution parameters for clustering. Default is c(0.1, 0.2, 0.3).
#' @param seed Random seed for reproducibility. Default is 1234.
#' @param computeMultiChain Boolean indicating whether to compute multi-chain graphs. Default is TRUE.
#' @param verbose Boolean indicating whether to print verbose debugging output. Default is FALSE.
#' @return Single Chain and multi-chain Seurat objects

.DistanceMatrixToClusteredGraphs <- function(seuratObj_TCR = NULL,
                                             pcaComponents = 50,
                                             kpcaKernel = "rbfdot",
                                             usePCA = TRUE,
                                             partitionType = "CPMVertexPartition",
                                             proportionOfGraphAsNeighbors = 0.1,
                                             jaccardIndexThreshold = 0.1,
                                             resolutions = c(0.1, 0.2, 0.3),
                                             seed = 1234,
                                             computeMultiChain = T,
                                             verbose = FALSE) {

  #check the Assays in the Seurat Object and compute graphs
  assays <- Seurat::Assays(seuratObj_TCR)

  #initialize composite object for multi-chain analysis
  seuratObj_TCR_composite <- NULL

  single_chain_graphs <- list()

  #calculate the single chain graphs
  for (assay in assays){
    if (verbose) print(assay)
    #get the distance matrix
    #handle Seurat v5 layer system - join layers if necessary
    tryCatch({
      #try to join layers to avoid v5 compatibility issues
      if (methods::is(seuratObj_TCR[[assay]], "Assay5")) {
        seuratObj_TCR <- SeuratObject::JoinLayers(seuratObj_TCR, assay = assay)
      }
    }, error = function(e) {
      #if JoinLayers fails or doesn't exist, continue without it
      warning(paste("Could not join layers for assay", assay, ":", e$message))
    })

    #try to get data without specifying layer first, then with layer if needed
    distanceMatrix <- tryCatch({
      as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay))
    }, error = function(e) {
      #if that fails, try with layer specification
      tryCatch({
        as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))
      }, error = function(e2) {
        #if both fail, try data slot
        as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, slot = "data"))
      })
    })
    graph_and_pca_results <- .PcaAndClustering(distanceMatrix = distanceMatrix,
                                               pcaComponents = pcaComponents,
                                               kpcaKernel = kpcaKernel,
                                               usePCA = usePCA,
                                               proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                               jaccardIndexThreshold = jaccardIndexThreshold)
    pruned_graph <- graph_and_pca_results$graph
    if (verbose) print(assay)
    single_chain_graphs[[assay]] <- pruned_graph
    for (resolution in resolutions) {
      partition <- leidenbase::leiden_find_partition(pruned_graph,
                                                     partition_type = partitionType,
                                                     initial_membership = NULL,
                                                     edge_weights = NULL,
                                                     node_sizes = NULL,
                                                     seed = seed,
                                                     resolution_parameter = resolution,
                                                     num_iter = 2,
                                                     verbose = verbose)
      #add the partition to the seurat object's metadata
      partition_metadata <- data.frame(partition$membership)
      rownames(partition_metadata) <- colnames(seuratObj_TCR[[assay]])

      seuratObj_TCR <- Seurat::AddMetaData(seuratObj_TCR,
                                           partition_metadata,
                                           col.name = paste0("TcrClustR_", assay, "_", resolution))
    }
    #add single-chain PCA/KPCA reductions and UMAPs
    reductionName <- paste0("TcrClustR_pca.", gsub("_",".", assay))
    pca_result <- graph_and_pca_results$pca_result
    seuratObj_TCR <- .AddDimensionalityReductions(seuratObj_TCR,
                                                  pca_result,
                                                  reductionName = reductionName,
                                                  assayName = assay,
                                                  distanceMatrix = distanceMatrix
    )
  }

  #TODO: this seems stylistically poor. A more upfront nested elif is probably better.
  #bail out of the multi-chain computation if it's not requested
  if (!computeMultiChain) {
    return(list(singleChainSeuratObject = seuratObj_TCR, multiChainSeuratObject = NULL))
  }

  #check if multichain assays actually exist before proceeding
  has_multichain_assays <- any(sapply(assays, function(assay) {
    #count the number of TCR chain identifiers in the assay name
    parts <- strsplit(assay, "_")[[1]]
    chain_count <- sum(parts %in% c("TRA", "TRB", "TRG", "TRD"))
    return(chain_count >= 2)
  }))

  if (!has_multichain_assays) {
    if (verbose) {
      print("No multichain assays found. Skipping multichain computation.")
      print("To compute multichain clustering, run RunTcrdist3 with multichain = TRUE.")
    }
    return(list(singleChainSeuratObject = seuratObj_TCR, multiChainSeuratObject = NULL))
  }

  #create a list of all possible combinations of chains
  chain_combinations <- utils::combn(names(single_chain_graphs), 2, simplify = FALSE)
  chain_combinations <- lapply(chain_combinations, function(x) paste(sort(x), collapse = "_"))
  #filter self-chain combinations
  chain_combinations <- chain_combinations[unlist(lapply(chain_combinations, FUN = function(x) {
    split <- strsplit(x, split = "_")[[1]]
    #check that the first three characters are not identical
    if (length(split) > 1 && substr(split[[1]], 1, 3) != substr(split[[2]], 1, 3)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }))]

  #calculate the combined chain graphs
  multi_chain_graphs <- list()
  for (joint_graph in chain_combinations) {
    #initialize the vectors to store group_by_variables and assays_to_access
    group_by_variables <- c()
    assays_to_access <- c()

    joint_graph <- gsub("_cdr3", "CDR3", joint_graph)

    #the graphs are named with the chains (TRA if V+J+CDR3, or TRACDR3 if only CDR3) separated by underscores
    #(e.g. TRACDR3_TRB for only TRA's CDR3 and V+J+CDR3 for the beta chain)
    #parse these strings and assign them to how they should inform the dplyr grouping
    #to determine unique combinations of TRA+TRB observed in the data.
    if (length(strsplit(joint_graph, "_")[[1]]) > 1 ) {
      chains <- strsplit(joint_graph, "_")[[1]]
    } else {
      chains <- joint_graph
    }
    if (verbose) print(chains)

    #get metadata once outside the loop
    metadata <- seuratObj_TCR@meta.data

    for (chain in chains) {
      #determine if the chain is CDR3-only and extract type
      is_cdr3_only <- grepl("CDR3$", chain)
      type <- sub("CDR3$", "", chain)

      #validate chain type
      valid_types <- c("TRA", "TRB", "TRG", "TRD")
      if (!type %in% valid_types) {
        stop(paste("Invalid chain type:", type))
      }

      if (is_cdr3_only) {
        # CDR3-only chain: group by CDR3 column, access CDR3 assay
        group_by_variables <- c(group_by_variables, type)
        names(group_by_variables)[length(group_by_variables)] <- paste0(type, "_cdr3")
        assays_to_access <- c(assays_to_access, paste0(type, "_cdr3"))
      } else {
        # Full chain: group by V/J/CDR3 columns, access main assay
        vj_columns <- c(paste0(type, "_V"), paste0(type, "_J"), type)
        group_by_variables <- c(group_by_variables, vj_columns)
        names(group_by_variables)[length(group_by_variables) - 2:0] <- c(paste0(type, "_V"), paste0(type, "_J"), type)
        assays_to_access <- c(assays_to_access, type)
      }
      if (verbose) print(paste0("group_variables:", group_by_variables))
    }

    #process multichain combinations after parsing all chains
    if (length(assays_to_access) > 1) {
      #translate group_by_variables to tcrdist3 column names for metadata access
      translated_group_by_variables <- .TranslateGroupByVariablesToTcrdist3(group_by_variables)
      names(translated_group_by_variables) <- names(group_by_variables)

      combined_matrix <- .ComputeMultiTCRDistanceMatrix(seuratObj_TCR = seuratObj_TCR,
                                                        group_by_variables = translated_group_by_variables,
                                                        assays_to_access = assays_to_access,
                                                        metadata = metadata,
                                                        verbose = verbose)

        graph_and_pca_results <- .PcaAndClustering(distanceMatrix = as.matrix(combined_matrix),
                                                   pcaComponents = pcaComponents,
                                                   kpcaKernel = kpcaKernel,
                                                   usePCA = usePCA,
                                                   proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                                   jaccardIndexThreshold = jaccardIndexThreshold)
        combined_graph <- graph_and_pca_results$graph

        multi_chain_graphs[[joint_graph]] <- combined_graph

        #create the composite TCR seuratObj distance object/append to the existing one
        if (is.null(seuratObj_TCR_composite)) {
          seuratObj_TCR_composite <- Seurat::CreateSeuratObject(counts = combined_matrix,
                                                                assay = joint_graph)
          seuratObj_TCR_composite$orig.ident <- joint_graph
          seuratObj_TCR_composite <- Seurat::AddMetaData(seuratObj_TCR_composite, rownames(combined_matrix), col.name = "composite_id")
        } else {
          #create a new assay in the existing composite object instead of merging separate objects
          seuratObj_TCR_composite[[joint_graph]] <- Seurat::CreateAssayObject(counts = combined_matrix)
        }

        #add multi-chain PCA/KPCA reductions and UMAPs
        reductionName <- paste0("TcrClustR_pca.", gsub("_",".", joint_graph))
        pca_result <- graph_and_pca_results$pca_result
        seuratObj_TCR_composite <- .AddDimensionalityReductions(seuratObj_TCR_composite,
                                                                pca_result,
                                                                reductionName,
                                                                assayName = joint_graph,
                                                                distanceMatrix = combined_matrix
        )

        for (resolution in resolutions) {
          partition <- leidenbase::leiden_find_partition(combined_graph,
                                                         partition_type = c("CPMVertexPartition", "ModularityVertexPartition", "RBConfigurationVertexPartition", "RBERVertexPartition", "SignificanceVertexPartition", "SurpriseVertexPartition"),
                                                         initial_membership = NULL,
                                                         edge_weights = NULL,
                                                         node_sizes = NULL,
                                                         seed = seed,
                                                         resolution_parameter = resolution,
                                                         num_iter = 2,
                                                         verbose = verbose)
          #add the partition to the seurat object's metadata
          partition_metadata <- data.frame(partition$membership)
          rownames(partition_metadata) <- colnames(seuratObj_TCR_composite[[joint_graph]])

          if (any(grepl("_cdr3",  joint_graph))) {
            joint_graph_name <- gsub("_cdr3", "CDR3", joint_graph)
          } else {
            joint_graph_name <- joint_graph
          }
          if (verbose) print(joint_graph_name)
          seuratObj_TCR_composite <- Seurat::AddMetaData(seuratObj_TCR_composite, partition_metadata, col.name = paste0("TcrClustR_", joint_graph_name, "_", resolution))
        }
      }
    }
  #rename the seuratObj_TCR object's assays to have parity with the multi-chain seuratObj_TCR_composite
  for (assay in Seurat::Assays(seuratObj_TCR)) {
    if (endsWith(assay, "_cdr3")) {
      new_assay_name <- gsub("_cdr3", "CDR3", assay)
      seuratObj_TCR <- SeuratObject::RenameAssays(seuratObj_TCR, assay = assay, new.assay.name = new_assay_name)
    }
  }
  return(list(singleChainSeuratObject = seuratObj_TCR, multiChainSeuratObject = seuratObj_TCR_composite))
}

.PcaAndClustering <- function(distanceMatrix = distanceMatrix,
                              pcaComponents = pcaComponents,
                              kpcaKernel = kpcaKernel,
                              usePCA = TRUE,
                              proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                              jaccardIndexThreshold = jaccardIndexThreshold){

  #validate distance matrix dimensions
  if (any(dim(distanceMatrix) == 0)) {
    stop(paste("Distance matrix has zero dimensions:", paste(dim(distanceMatrix), collapse = "x")))
  }

  #validate number of components vs distance matrix dimensions
  if (pcaComponents > ncol(distanceMatrix)) {
   warning("Number of requested components exceeds available dimensions. Setting pcaComponents to the number of columns in the distance matrix.")
   pcaComponents <- ncol(distanceMatrix)
  }

  if (usePCA) {
    #use standard PCA
    pca_result <- stats::prcomp(distanceMatrix, center = TRUE, scale. = TRUE)

    #create a compatible object structure similar to kernlab::kpca, so downstream parsing can be identical
    pca_result_obj <- list(
      rotated = pca_result$x,
      sdev = pca_result$sdev
    )
    class(pca_result_obj) <- "pca_result"

    #add methods for compatibility
    attr(pca_result_obj, "rotated") <- pca_result$x

  } else {
    #otherwise, kernel PCA.
    pca_result_obj <- kernlab::kpca(x = distanceMatrix, kernel = kpcaKernel)
  }

  #get the rotated data from the kernlab compatible object
  if (usePCA) {
    rotated_data <- pca_result_obj$rotated
  } else {
    rotated_data <- kernlab::rotated(pca_result_obj)
  }

  #reduce the data to the first n_components
  n_components <- min(c(pcaComponents, nrow(distanceMatrix), ncol(rotated_data)))
  reduced_data <- rotated_data[, 1:n_components, drop = FALSE]

  #take 10% of the graph as nearest neighbors, in the style of conga by default
  k <-  round(proportionOfGraphAsNeighbors * ncol(distanceMatrix))
  #validate the number of neighbors
  if (k < 1) k <- 1
  if (k >= nrow(distanceMatrix)) k <- max(1, nrow(distanceMatrix) - 1)

  knn_result <- FNN::get.knn(reduced_data, k = k)

  #make graph
  edges <- cbind(rep(seq_len(nrow(distanceMatrix)), each = k), c(knn_result$nn.index))
  g <- igraph::graph_from_edgelist(edges, directed = FALSE)
  #remove loops
  g <- igraph::simplify(g)
  #prune graph using jaccard index thresholding on the snn graph, similar to Seurat.
  jaccard_index <- igraph::similarity(g, vids = igraph::V(g), mode = 'all')
  adj_matrix <- igraph::as_adjacency_matrix(g)
  #prune edges below the jaccard index threshold
  adj_matrix[jaccard_index < jaccardIndexThreshold] <- 0

  #TODO: evaluate some of these parameters in different contexts ('barnyard' experiment on TCRs?)
  pruned_graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = 'undirected')
  return(list(graph = pruned_graph, pca_result = pca_result_obj))
}

.ComputeMultiTCRDistanceMatrix <- function(seuratObj_TCR = NULL,
                                           group_by_variables = NULL,
                                           assays_to_access = NULL,
                                           metadata = NULL,
                                           verbose = FALSE) {

  if (verbose) {
    print(paste("Debug: group_by_variables =", paste(group_by_variables, collapse = ", ")))
    print(paste("Debug: assays_to_access =", paste(assays_to_access, collapse = ", ")))
    print(paste("Debug: metadata columns =", paste(colnames(metadata), collapse = ", ")))
  }

  # Generate the required columns for observed_tcr_pairs using the SAME logic as .CreateTcrKeyLookup
  # This ensures consistency between observed_tcr_pairs and lookup tables

  # Use the same logic as .CreateTcrKeyLookup to determine required columns
  get_required_cols <- function(assay_name) {
    is_cdr3_assay <- endsWith(assay_name, "_cdr3")
    chain_type <- gsub("_cdr3", "", assay_name)

    if(is_cdr3_assay) {
      return(.TranslateGroupByVariablesToTcrdist3(chain_type))
    } else {
      return(.TranslateGroupByVariablesToTcrdist3(c(
        paste0(chain_type, "_V"),
        paste0(chain_type, "_J"),
        chain_type
      )))
    }
  }

  required_cols_first <- get_required_cols(assays_to_access[1])
  required_cols_second <- get_required_cols(assays_to_access[2])
  all_required_cols <- unique(c(required_cols_first, required_cols_second))

  if (verbose) {
    print(paste("Debug: required columns for observed_tcr_pairs:", paste(all_required_cols, collapse = ", ")))
  }

  observed_tcr_pairs <- metadata |>
    dplyr::select(dplyr::all_of(all_required_cols)) |>
    dplyr::filter_all(dplyr::all_vars(!is.na(.))) |>
    dplyr::filter_all(dplyr::all_vars(!grepl(",",.))) |>
    dplyr::distinct()

  if (verbose) {
    print(paste("Debug: observed_tcr_pairs has", nrow(observed_tcr_pairs), "rows"))
  }

  #handle Seurat v5 layer system for both assays
  for (assay in assays_to_access) {
    tryCatch({
      if (methods::is(seuratObj_TCR[[assay]], "Assay5")) {
        seuratObj_TCR <- SeuratObject::JoinLayers(seuratObj_TCR, assay = assay)
      }
    }, error = function(e) {
      warning(paste("Could not join layers for assay", assay, ":", e$message))
    })
  }

  #get matrices with fallback for Seurat v5 compatibility
  first_chain_matrix <- tryCatch({
    Seurat::GetAssayData(seuratObj_TCR, assay = assays_to_access[1])
  }, error = function(e) {
    tryCatch({
      Seurat::GetAssayData(seuratObj_TCR, assay = assays_to_access[1], layer = "counts")
    }, error = function(e2) {
      Seurat::GetAssayData(seuratObj_TCR, assay = assays_to_access[1], slot = "data")
    })
  })

  second_chain_matrix <- tryCatch({
    Seurat::GetAssayData(seuratObj_TCR, assay = assays_to_access[2])
  }, error = function(e) {
    tryCatch({
      Seurat::GetAssayData(seuratObj_TCR, assay = assays_to_access[2], layer = "counts")
    }, error = function(e2) {
      Seurat::GetAssayData(seuratObj_TCR, assay = assays_to_access[2], slot = "data")
    })
  })

  #create lookup tables for both chains
  first_chain_lookup <- .CreateTcrKeyLookup(seuratObj_TCR, assays_to_access[1])
  second_chain_lookup <- .CreateTcrKeyLookup(seuratObj_TCR, assays_to_access[2])

  if (verbose) {
    print(paste("first_chain_lookup has", nrow(first_chain_lookup), "rows"))
    print(paste("second_chain_lookup has", nrow(second_chain_lookup), "rows"))
    print(paste("observed_tcr_pairs columns =", paste(colnames(observed_tcr_pairs), collapse = ", ")))
  }

  #create keys for observed_tcr_pairs using the same logic as .CreateTcrKeyLookup
  #apply the same key generation logic as .CreateTcrKeyLookup for both chains

  #for first chain - apply same logic as .CreateTcrKeyLookup
  if(grepl("_cdr3$", assays_to_access[1])) {
    # CDR3-only first chain - extract the first column (CDR3 sequence)
    first_chain_key_col <- required_cols_first[1]
    observed_pairs_with_keys <- observed_tcr_pairs %>%
      dplyr::mutate(first_chain_key = .data[[first_chain_key_col]])
  } else {
    #full first chain - paste V, J, CDR3 with underscore separator
    observed_pairs_with_keys <- observed_tcr_pairs %>%
      dplyr::mutate(first_chain_key = paste(
        .data[[required_cols_first[1]]],  # V gene
        .data[[required_cols_first[2]]],  # J gene
        .data[[required_cols_first[3]]],  # CDR3 sequence
        sep = "_"
      ))
  }

  #for second chain - apply same logic as .CreateTcrKeyLookup
  if(grepl("_cdr3$", assays_to_access[2])) {
    # CDR3-only second chain - extract the first column (CDR3 sequence)
    second_chain_key_col <- required_cols_second[1]
    observed_pairs_with_keys <- observed_pairs_with_keys %>%
      dplyr::mutate(second_chain_key = .data[[second_chain_key_col]])
  } else {
    #full second chain - paste V, J, CDR3 with underscore separator
    observed_pairs_with_keys <- observed_pairs_with_keys %>%
      dplyr::mutate(second_chain_key = paste(
        .data[[required_cols_second[1]]],  # V gene
        .data[[required_cols_second[2]]],  # J gene
        .data[[required_cols_second[3]]],  # CDR3 sequence
        sep = "_"
      ))
  }

  #remove allele annotations (same as in .CreateTcrKeyLookup) to match lookup table format
  observed_pairs_with_keys <- observed_pairs_with_keys %>%
    dplyr::mutate(
      first_chain_key = gsub("\\*01", "", .data$first_chain_key),
      second_chain_key = gsub("\\*01", "", .data$second_chain_key)
    )

  if (verbose) {
    #show sample keys from both sides
    print("Sample keys from observed_pairs_with_keys:")
    print(head(observed_pairs_with_keys[c("first_chain_key", "second_chain_key")], 3))
    print("Sample keys from first_chain_lookup:")
    print(head(first_chain_lookup$key, 3))
    print("Sample keys from second_chain_lookup:")
    print(head(second_chain_lookup$key, 3))
  }

  #map keys to matrix row names using lookups
  valid_pairs <- observed_pairs_with_keys %>%
    dplyr::left_join(first_chain_lookup, by = c("first_chain_key" = "key")) %>%
    dplyr::left_join(second_chain_lookup, by = c("second_chain_key" = "key")) %>%
    dplyr::filter(!is.na(matrix_rowname.x) & !is.na(matrix_rowname.y))

  if (verbose) {
    print(paste("valid_pairs has", nrow(valid_pairs), "rows after join and filter"))
  }

  if (nrow(valid_pairs) == 0) {
    stop("No valid TCR pairs found for multi-chain analysis. Check that metadata columns exist and contain matching data.")
  }

  #create a composite ID for each pair, however this ID needs to map to the "cell barcode" version of the TCR to map with the metadata, rather than the "features" version of the TCR
  #for seurat reasons, the "cellbarcode" supports underscores, and the 'feature' supports hypens.
  #to maintain this when we create the seurat object, we'll use "delimit" as a row/column agnostic delimiter.
  valid_pairs <- valid_pairs %>%
    dplyr::mutate(composite_id = paste0(gsub("-", "_", matrix_rowname.x), "delimit", gsub("-", "_", matrix_rowname.y)))

  #define matrix indices for valid pairs
  first_chain_indices <- match(valid_pairs$matrix_rowname.x, rownames(first_chain_matrix))
  second_chain_indices <- match(valid_pairs$matrix_rowname.y, rownames(second_chain_matrix))

  #pairwise combination of valid pairs
  combos <- expand.grid(i = seq_along(rownames(valid_pairs)), j = seq_along(rownames(valid_pairs)))

  distances <- apply(combos, 1, function(pair) {
    first_dist <- first_chain_matrix[first_chain_indices[pair[1]], first_chain_indices[pair[2]]]
    second_dist <- second_chain_matrix[second_chain_indices[pair[1]], second_chain_indices[pair[2]]]
    first_dist + second_dist
  })
  #create the joint-distance matrix
  combined_matrix <- Matrix::sparseMatrix(
    i = combos$i,
    j = combos$j,
    x = distances,
    dims = c(nrow(valid_pairs), nrow(valid_pairs))
  )
  rownames(combined_matrix) <- colnames(combined_matrix) <- valid_pairs$composite_id

  return(combined_matrix)
}


.TranslateGroupByVariablesToTcrdist3 <- function(group_by_variables){
  group_variables_tcrdist3 <- gsub("TRA_V", "v_a_gene", group_by_variables)
  group_variables_tcrdist3 <- gsub("TRB_V", "v_b_gene", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRG_V", "v_g_gene", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRD_V", "v_d_gene", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRA_J", "j_a_gene", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRB_J", "j_b_gene", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRG_J", "j_g_gene", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRD_J", "j_d_gene", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRA", "cdr3_a_aa", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRB", "cdr3_b_aa", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRG", "cdr3_g_aa", group_variables_tcrdist3)
  group_variables_tcrdist3 <- gsub("TRD", "cdr3_d_aa", group_variables_tcrdist3)
  return(group_variables_tcrdist3)
}


#helper function to create key->rowname mappings for each assay
.CreateTcrKeyLookup <- function(seuratObj_TCR, assay_name) {
  #get metadata columns for this assay type
  is_cdr3_assay <- endsWith(assay_name, "_cdr3")
  chain_type <- gsub("_cdr3", "", assay_name)

  #get translated column names from metadata
  required_cols <- if(is_cdr3_assay) {
    tcrdist3_cols <- .TranslateGroupByVariablesToTcrdist3(chain_type)
    c(tcrdist3_cols)
  } else {
    tcrdist3_cols <- .TranslateGroupByVariablesToTcrdist3(c(
      paste0(chain_type, "_V"),
      paste0(chain_type, "_J"),
      chain_type
    ))
    tcrdist3_cols
  }

  #extract relevant metadata and create composite keys
  metadata <- seuratObj_TCR@meta.data[, required_cols, drop = FALSE]

  keys <- if(is_cdr3_assay) {
    metadata[[1]] #CDR3 sequence
  } else {
    paste(metadata[[1]], metadata[[2]], metadata[[3]], sep = "_")
  }

  #TODO: this (hack) assumes that the allele notation is entirely ceremonial and carries no information.
  #TODO: support alleles? unsure how to do this in the current compute environment though.
  keys <- gsub("\\*01", "", keys)

  matrix_rownames <- tryCatch({
    rownames(Seurat::GetAssayData(seuratObj_TCR, assay = assay_name))
  }, error = function(e) {
    tryCatch({
      rownames(Seurat::GetAssayData(seuratObj_TCR, assay = assay_name, layer = 'counts'))
    }, error = function(e2) {
      rownames(Seurat::GetAssayData(seuratObj_TCR, assay = assay_name, slot = 'data'))
    })
  })

  #create lookup table: key -> matrix row name
  lookup <- data.frame(
    key = keys,
    matrix_rowname = matrix_rownames,
    stringsAsFactors = FALSE
  )
  #distinct in case of duplicates, but should be unnecessary.
  lookup <- dplyr::distinct(lookup, key, .keep_all = TRUE)

  return(lookup)
}

.AddDimensionalityReductions <- function(seuratObj,
                                         pca_result = NULL,
                                         reductionName = NULL,
                                         assayName = NULL,
                                         pcaComponents = 50,
                                         kpcaKernel = 'rbfdot',
                                         proportionOfGraphAsNeighbors = 0.1,
                                         jaccardIndexThreshold = 0.1,
                                         distanceMatrix = NULL) {
  #add PCA/KPCA components and make UMAP
  # Handle both PCA and kernel PCA results
  if (inherits(pca_result, "pca_result")) {
    embeddings <- pca_result$rotated
  } else {
    embeddings <- pca_result@rotated
  }

  rownames(embeddings) <- colnames(seuratObj[[assayName]])
  seuratObj[[reductionName]] <-  Seurat::CreateDimReducObject(embeddings = embeddings,
                                                              assay = assayName,
                                                              key = paste0(reductionName, "_"))
  #take 10% of the graph as nearest neighbors, in the style of conga by default
  k.param <-  round(proportionOfGraphAsNeighbors * ncol(distanceMatrix))

  #TODO: compute a cutoff for the number of components used, but I'm not sure what these distributions look like yet.
  # Handle both PCA and kernel PCA for getting number of components
  if (inherits(pca_result, "pca_result")) {
    n_components = min(c(pcaComponents, nrow(embeddings), ncol(pca_result$rotated)))
  } else {
    n_components = min(c(pcaComponents, nrow(embeddings), ncol(pca_result@rotated)))
  }
  seuratObj <- Seurat::FindNeighbors(seuratObj,
                                     reduction = reductionName,
                                     dims = 1:n_components,
                                     k.param = k.param
  )
  seuratObj <- Seurat::RunUMAP(seuratObj,
                               dims = 1:n_components,
                               reduction = reductionName,
                               umap.method = "uwot",
                               n.neighbors = min(c(nrow(embeddings),30)),
                               reduction.name = paste0(reductionName, "_umap"))
  return(seuratObj)
}
