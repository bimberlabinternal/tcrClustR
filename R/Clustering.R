
utils::globalVariables(
  names = c('matrix_rowname.x', 'matrix_rowname.y', '.data', '.', 'distanceMatrix', 'combined_matrix',
            'key', 'seed', 'seuratObj', 'spikeInDataframe'),
  package = 'tcrClustR',
  add = TRUE
)

#TODO: summarize data lossy-ness vignette.

ClusterTcrs <- function(seuratObj = NULL,
                        seuratObj_TCR = NULL,
                        metadata = NULL,
                        resolutionParameter = 0.1,
                        kpcaComponents = 50,
                        kpcaKernel = "rbfdot",
                        proportionOfGraphAsNeighbors = 0.1,
                        jaccardIndexThreshold = 0.1,
                        seed = 1234,
                        spikeInDataframe =  NULL,
                        computeMultiChain = T) {

  #perform leiden clustering on the single chain distance matrices, create the multichain distance matrices, and cluster them.
  clusteredSeuratObjects <- .DistanceMatrixToClusteredGraphs(seuratObj_TCR = seuratObj_TCR,
                                                             kpcaComponents = kpcaComponents,
                                                             kpcaKernel = kpcaKernel,
                                                             proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                                             jaccardIndexThreshold = jaccardIndexThreshold)
  #Parse the single chain and multi-chain seurat objects, iterate through the assays, and assign cells in the original seuratObj to various clusters
  for (tcr_object in clusteredSeuratObjects) {
    for (assay in SeuratObject::Assays(tcr_object)) {
      print(assay)
      #detect multichain assays
      if (!grepl("_", assay)) {
        #parse out the group_by_variables
        group_by_variables <- c()
        if (endsWith(assay, "CDR3")) {
          group_by_variables <- gsub("CDR3", "", assay)
        } else {
          assayname <- gsub("CDR3", "", assay)
          group_by_variables <- c(paste0(assayname, "_V"), paste0(assayname, "_J"), assayname)
        }
        print(paste0("single chain assay:", assay))
        print(group_by_variables)
        singleChainMetadata <- tcr_object@meta.data
      } else {
        #parse out the group_by_variables
        group_by_variables <- c()
        chains <- strsplit(assay, "_")[[1]]
        for (chain in chains) {
          if (endsWith(chain, "CDR3")) {
            group_by_variables <- c(group_by_variables, gsub("CDR3", "", chain))
          } else {
            chain_name <- gsub("CDR3", "", chain)
            group_by_variables <- c(group_by_variables, paste0(chain_name, "_V"), paste0(chain_name, "_J"), gsub("CDR3", "", chain_name))
          }
        }
        print(paste0("multichain assay:", chains))
        print(group_by_variables)
      }
    }

  }
}



#' @title .DistanceMatrixToClusteredGraphs
#' @description
#' This function takes a Seurat object with TCR distance matrices and computes clustered graphs for each chain.
#' @param seuratObj_TCR Seurat object containing TCR distance matrices.
#' @param kpcaComponents Number of components for kernel PCA. Default is 50.
#' @param kpcaKernel Kernel type for kernel PCA. Default is "rbfdot".
#' @param partitionType Type of partitioning algorithm to use. Default is "CPMVertexPartition".
#' @param proportionOfGraphAsNeighbors Proportion of the graph to consider as neighbors. Default is 0.1.
#' @param jaccardIndexThreshold Jaccard index threshold for pruning edges. Default is 0.1.
#' @param resolutions Vector of resolution parameters for clustering. Default is c(0.1, 0.2, 0.3).
#' @param computeMultiChain Boolean indicating whether to compute multi-chain graphs. Default is TRUE.
#' @return Single Chain and multi-chain Seurat objects

.DistanceMatrixToClusteredGraphs <- function(seuratObj_TCR = NULL,
                                             kpcaComponents = 50,
                                             kpcaKernel = "rbfdot",
                                             partitionType = "CPMVertexPartition",
                                             proportionOfGraphAsNeighbors = 0.1,
                                             jaccardIndexThreshold = 0.1,
                                             resolutions = c(0.1, 0.2, 0.3),
                                             computeMultiChain = T) {

  #check the Assays in the Seurat Object and compute graphs
  assays <- Seurat::Assays(seuratObj_TCR)

  single_chain_graphs <- list()

  #calculate the single chain graphs
  for (assay in assays){
    print(assay)
    #get the distance matrix
    distanceMatrix <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))
    graph_and_kpca_results <- .KpcaAndClustering(distanceMatrix = distanceMatrix,
                                                 kpcaComponents = kpcaComponents,
                                                 kpcaKernel = kpcaKernel,
                                                 proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                                 jaccardIndexThreshold = jaccardIndexThreshold)
    pruned_graph <- graph_and_kpca_results$graph
    print(assay)
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
                                                     verbose = TRUE)
      #add the partition to the seurat object's metadata
      partition_metadata <- data.frame(partition$membership)
      rownames(partition_metadata) <- colnames(seuratObj_TCR[[assay]])

      seuratObj_TCR <- Seurat::AddMetaData(seuratObj_TCR,
                                           partition_metadata,
                                           col.name = paste0("TcrClustR_", assay, "_", resolution))
    }
    #add single-chain KPCA reductions and UMAPs
    reductionName <- paste0("TcrClustR_kpca.", gsub("_",".", assay))
    kpca_result <- graph_and_kpca_results$kpca_result
    seuratObj_TCR <- .AddDimensionalityReductions(seuratObj_TCR,
                                                  kpca_result,
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
    #initialize the vectors to store whether or not the chain is only CDR3
    #or if the V/J segments should be used.
    cdr3_only_chains <- c()
    remaining_chains <- c()
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
    print(chains)
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
      print(paste0("group_variables:", group_by_variables))
      #iterate through the 10X data, merge with the spike-in dataframes and index the metadata by the:
      # 1. observed TRA+TRB combinations in the 10X data
      # 2. provided TRA+TRB combinations in the spikeInData
      # TODO: make this work with a metadata dataframe instead of only with a seurat object
      metadata <- plyr::rbind.fill(seuratObj@meta.data, spikeInDataframe)

      #if there are multiple chains, figure out how to combine them
      if (length(assays_to_access) > 1) {

        combined_matrix <- .ComputeMultiTCRDistanceMatrix(seuratObj_TCR = seuratObj_TCR,
                                                          group_by_variables = group_by_variables,
                                                          assays_to_access = assays_to_access,
                                                          metadata = metadata)

        graph_and_kpca_results <- .KpcaAndClustering(distanceMatrix = as.matrix(combined_matrix),
                                                     kpcaComponents = kpcaComponents,
                                                     kpcaKernel = kpcaKernel,
                                                     proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                                     jaccardIndexThreshold = jaccardIndexThreshold)
        combined_graph <- graph_and_kpca_results$graph

        multi_chain_graphs[[joint_graph]] <- combined_graph

        #create the composite TCR seuratObj distance object/append to the existing one
        if (is.null(seuratObj_TCR_composite)) {
          seuratObj_TCR_composite <- Seurat::CreateSeuratObject(counts = combined_matrix,
                                                                assay = joint_graph)
          seuratObj_TCR_composite$orig.ident <- joint_graph
          seuratObj_TCR_composite <- Seurat::AddMetaData(seuratObj_TCR_composite, rownames(combined_matrix), col.name = "composite_id")
        } else {
          seuratObj_TCR_composite_subsequent_chain_combination <- Seurat::CreateSeuratObject(counts = combined_matrix,
                                                                                             assay = joint_graph)
          seuratObj_TCR_composite_subsequent_chain_combination$orig.ident <- joint_graph
          seuratObj_TCR_composite_subsequent_chain_combination <- Seurat::AddMetaData(seuratObj_TCR_composite_subsequent_chain_combination, rownames(combined_matrix), col.name = "composite_id")
          embeddings <- kpca_result@rotated
          colnames(embeddings) <- paste0("TcrClustR_kpca.", gsub("_",".", joint_graph), "-", seq_len(ncol(embeddings)))
          seuratObj_TCR_composite[[paste0("TcrClustR_kpca.", gsub("_",".", joint_graph))]] <-  Seurat::CreateDimReducObject(embeddings = embeddings,
                                                                                                                            assay = joint_graph,
                                                                                                                            key = "KPCA_")
          seuratObj_TCR_composite <- merge(seuratObj_TCR_composite, seuratObj_TCR_composite_subsequent_chain_combination)
        }

        #add multi-chain KPCA reductions and UMAPs
        reductionName <- paste0("TcrClustR_kpca.", gsub("_",".", joint_graph))
        kpca_result <- graph_and_kpca_results$kpca_result
        assayName <- joint_graph
        seuratObj_TCR_composite <- .AddDimensionalityReductions(seuratObj_TCR_composite,
                                                                kpca_result,
                                                                reductionName,
                                                                assayName = joint_graph,
                                                                distanceMatrix = distanceMatrix
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
                                                         verbose = TRUE)
          #add the partition to the seurat object's metadata
          partition_metadata <- data.frame(partition$membership)
          rownames(partition_metadata) <- colnames(seuratObj_TCR_composite[[joint_graph]])

          if (any(grepl("_cdr3",  joint_graph))) {
            joint_graph_name <- gsub("_cdr3", "CDR3", joint_graph)
          } else {
            joint_graph_name <- joint_graph
          }
          print(joint_graph_name)
          seuratObj_TCR_composite <- Seurat::AddMetaData(seuratObj_TCR_composite, partition_metadata, col.name = paste0("TcrClustR_", joint_graph_name, "_", resolution))
        }
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

.KpcaAndClustering <- function(distanceMatrix = distanceMatrix,
                               kpcaComponents = kpcaComponents,
                               kpcaKernel = kpcaKernel,
                               proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                               jaccardIndexThreshold = jaccardIndexThreshold){
  #perform kernel PCA in the same way that conga does
  kpca_result <- kernlab::kpca(x = distanceMatrix, kernel = kpcaKernel)

  #reduce the data to the first n_components
  n_components <- min( c(kpcaComponents, nrow(distanceMatrix), ncol(kernlab::rotated(kpca_result))))
  reduced_data <- kernlab::rotated(kpca_result)[, 1:n_components]

  #take 10% of the graph as nearest neighbors, in the style of conga by default
  k <-  round(proportionOfGraphAsNeighbors * ncol(distanceMatrix))
  knn_result <- FNN::get.knn(reduced_data, k = round(k))

  #make graph
  edges <- cbind(rep(1:nrow(distanceMatrix), each = k), c(knn_result$nn.index))
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
  return(list(graph = pruned_graph, kpca_result = kpca_result))
}

.ComputeMultiTCRDistanceMatrix <- function(seuratObj_TCR = NULL,
                                           group_by_variables = NULL,
                                           assays_to_access = NULL,
                                           metadata = NULL) {

  observed_tcr_pairs <- metadata |>
    dplyr::select(dplyr::all_of(group_by_variables)) |>
    dplyr::filter_all(dplyr::all_vars(!is.na(.))) |>
    dplyr::filter_all(dplyr::all_vars(!grepl(",",.))) |>
    dplyr::distinct()

  first_chain_matrix <- Seurat::GetAssayData(seuratObj_TCR, assay = assays_to_access[1], layer = "counts")
  second_chain_matrix <- Seurat::GetAssayData(seuratObj_TCR, assay = assays_to_access[2], layer = "counts")

  #populate all possible combinations of metadata features
  first_chain_variables <- seuratObj_TCR@meta.data[,.TranslateGroupByVariablesToTcrdist3(group_by_variables[names(group_by_variables) == assays_to_access[1]]), drop = FALSE]
  second_chain_variables <- seuratObj_TCR@meta.data[,.TranslateGroupByVariablesToTcrdist3(group_by_variables[names(group_by_variables) == assays_to_access[2]]), drop = FALSE]

  #create lookup tables for both chains
  first_chain_lookup <- .CreateTcrKeyLookup(seuratObj_TCR, assays_to_access[1])
  second_chain_lookup <- .CreateTcrKeyLookup(seuratObj_TCR, assays_to_access[2])

  first_chain_type <- gsub("_cdr3$", "", assays_to_access[1])
  second_chain_type <- gsub("_cdr3$", "", assays_to_access[2])
  #generate keys for observed pairs
  observed_pairs_with_keys <- observed_tcr_pairs %>%
    dplyr::mutate(
      first_chain_key = if(grepl("_cdr3$", assays_to_access[1])) {
        .data[[names(group_by_variables[names(group_by_variables) == assays_to_access[1]])]]
      } else {
        paste(
          .data[[paste0(first_chain_type, "_V")]],
          .data[[paste0(first_chain_type, "_J")]],
          .data[[first_chain_type]],
          sep = "_"
        )
      },
      second_chain_key = if(grepl("_cdr3$", assays_to_access[2])) {
        .data[[names(group_by_variables[names(group_by_variables) == assays_to_access[2]])]]
      } else {
        paste(
          .data[[paste0(second_chain_type, "_V")]],
          .data[[paste0(second_chain_type, "_J")]],
          .data[[second_chain_type]],
          sep = "_"
        )
      }
    )

  #map keys to matrix row names using lookups
  valid_pairs <- observed_pairs_with_keys %>%
    dplyr::left_join(first_chain_lookup, by = c("first_chain_key" = "key")) %>%
    dplyr::left_join(second_chain_lookup, by = c("second_chain_key" = "key")) %>%
    dplyr::filter(!is.na(matrix_rowname.x) & !is.na(matrix_rowname.y))

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

  matrix_rownames <- rownames(Seurat::GetAssayData(seuratObj_TCR, assay = assay_name, layer = 'counts'))

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
                                         kpca_result = NULL,
                                         reductionName = NULL,
                                         assayName = NULL,
                                         kpcaComponents = 50,
                                         kpcaKernel = 'rbfdot',
                                         proportionOfGraphAsNeighbors = 0.1,
                                         jaccardIndexThreshold = 0.1,
                                         distanceMatrix = NULL) {
  #add KPCA components and make UMAP
  embeddings <- kpca_result@rotated

  rownames(embeddings) <- paste0(colnames(seuratObj[[assayName]]))
  seuratObj[[reductionName]] <-  Seurat::CreateDimReducObject(embeddings = embeddings,
                                                              assay = assayName,
                                                              key = paste0(reductionName, "_"))
  #take 10% of the graph as nearest neighbors, in the style of conga by default
  k.param <-  round(proportionOfGraphAsNeighbors * ncol(distanceMatrix))

  #TODO: compute a cutoff for the number of components used, but I'm not sure what these distributions look like yet.
  n_components = min(c(kpcaComponents, nrow(embeddings), ncol(kpca_result@rotated)))
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
