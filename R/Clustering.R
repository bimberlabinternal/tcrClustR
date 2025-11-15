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
#' @param resolutionParameters Vector of resolution parameters for clustering. Use to iterate over multiple resolutions.
#' @param pcaComponents Number of components for PCA or kernel PCA. Default is 50.
#' @param kpcaKernel Kernel type for kernel PCA. Default is "rbfdot". Ignored if usePCA is TRUE.
#' @param usePCA Boolean indicating whether to use standard PCA instead of kernel PCA. Default is FALSE.
#' @param partitionType Type of Leiden partitioning algorithm to use. Default is "RBConfigurationVertexPartition".
#' @param proportionOfGraphAsNeighbors Proportion of the graph to consider as neighbors. Default is 0.1.
#' @param jaccardIndexThreshold Jaccard index threshold for pruning edges. Default is 0.1.
#' @param seed Random seed for reproducibility. Default is 1234.
#' @param computeMultiChain Boolean indicating whether to compute multi-chain graphs. Default is TRUE.
#' @param verbose Boolean indicating whether to print verbose debugging output. Default is FALSE.
#' @param neighborsK Fixed number of nearest neighbors for graph construction (Scanpy/Seurat style). Default is 15.
#' @param rbfSigma Width for RBF kernel when converting distances to similarities for KernelPCA and edge weights. Default NULL (heuristic).
#' @param useExactDistanceKNN Boolean indicating whether to build kNN directly from distances (skip KPCA for kNN step). Default is FALSE.
#' @param edgeWeightMode Edge weighting mode: one of "jaccard" (SNN Jaccard), "rbf" (RBF similarity on SNN support), or "jaccard_x_rbf" (product). Default is "jaccard_x_rbf".
#' @param useHdbscanTwoStage Boolean indicating whether to run two-stage HDBSCAN with silhouette-tuned minPts on selected assays. Default is FALSE.
#' @param hdbscanAssays Character vector of assay names to run HDBSCAN on when useHdbscanTwoStage is TRUE. Default is c("TRB_cdr3").
#' @param minPtsRangeNoise Integer vector of minPts values to scan for the noise-removal stage. Default is 2:20.
#' @param minPtsRangeFidelity Integer vector of minPts values to scan for the fidelity stage. Default is 2:20.
#' @return A list containing single-chain and multi-chain Seurat objects with clustering results in their metadata.
#'   Users should perform clonotypic joins to transfer clustering information to their original Seurat object.
#' @export


ClusterTcrs <- function(seuratObj = NULL,
                        seuratObj_TCR = NULL,
                        metadata = NULL,
                        resolutionParameters = c(0.1, 0.2, 0.3),
                        pcaComponents = 50,
                        kpcaKernel = "rbfdot",
                        usePCA = FALSE,
                        partitionType = "RBConfigurationVertexPartition",
                        proportionOfGraphAsNeighbors = 0.1,
                        jaccardIndexThreshold = 0.1,
                        seed = 1234,
                        computeMultiChain = T,
                        verbose = FALSE,
                        neighborsK = 15,
                        rbfSigma = NULL,
                        useExactDistanceKNN = FALSE,
                        edgeWeightMode = "jaccard_x_rbf",
                        useHdbscanTwoStage = FALSE,
                        hdbscanAssays = c("TRB_cdr3"),
                        minPtsRangeNoise = 2:20,
                        minPtsRangeFidelity = 2:20) {

  # normalize resolution parameters (must be numeric vector)
  resolutions_to_use <- as.numeric(resolutionParameters)
  if (any(!is.finite(resolutions_to_use))) {
    stop("resolutionParameters must be numeric.")
  }

  #perform leiden clustering on the single chain distance matrices, create the multichain distance matrices, and cluster them.
  clusteredSeuratObjects <- .DistanceMatrixToClusteredGraphs(seuratObj_TCR = seuratObj_TCR,
                                                             pcaComponents = pcaComponents,
                                                             kpcaKernel = kpcaKernel,
                                                             usePCA = usePCA,
                                                             partitionType = partitionType,
                                                             proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                                             jaccardIndexThreshold = jaccardIndexThreshold,
                                                             seed = seed,
                                                             computeMultiChain = computeMultiChain,
                                                             verbose = verbose,
                                                             neighborsK = neighborsK,
                                                             rbfSigma = rbfSigma,
                                                             useExactDistanceKNN = useExactDistanceKNN,
                                                             edgeWeightMode = edgeWeightMode,
                                                             resolutions = resolutions_to_use,
                                                             useHdbscanTwoStage = useHdbscanTwoStage,
                                                             hdbscanAssays = hdbscanAssays,
                                                             minPtsRangeNoise = minPtsRangeNoise,
                                                             minPtsRangeFidelity = minPtsRangeFidelity)

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
#' @param usePCA Boolean indicating whether to use standard PCA instead of kernel PCA. Default is FALSE.
#' @param partitionType Type of partitioning algorithm to use. Default is "RBConfigurationVertexPartition".
#' @param proportionOfGraphAsNeighbors Proportion of the graph to consider as neighbors. Default is 0.1.
#' @param jaccardIndexThreshold Jaccard index threshold for pruning edges. Default is 0.1.
#' @param resolutions Vector of resolution parameters for clustering. Default is c(0.1, 0.2, 0.3).
#' @param seed Random seed for reproducibility. Default is 1234.
#' @param computeMultiChain Boolean indicating whether to compute multi-chain graphs. Default is TRUE.
#' @param neighborsK Fixed number of nearest neighbors for graph construction (Scanpy/Seurat style). Default is 15.
#' @param rbfSigma Width for RBF kernel when converting distances to similarities for KernelPCA and edge weights. Default NULL (heuristic).
#' @param useExactDistanceKNN Boolean indicating whether to build kNN directly from distances (skip KPCA for kNN step). Default is FALSE.
#' @param edgeWeightMode Edge weighting mode: one of "jaccard", "rbf", or "jaccard_x_rbf". Default is "jaccard_x_rbf".
#' @param useHdbscanTwoStage Boolean indicating whether to run two-stage HDBSCAN with silhouette-tuned minPts on selected assays. Default is FALSE.
#' @param hdbscanAssays Character vector of assay names to run HDBSCAN on when useHdbscanTwoStage is TRUE. Default is c("TRB_cdr3").
#' @param minPtsRangeNoise Integer vector of minPts values to scan for the noise-removal stage. Default is 2:20.
#' @param minPtsRangeFidelity Integer vector of minPts values to scan for the fidelity stage. Default is 2:20.
#' @return Single Chain and multi-chain Seurat objects

.DistanceMatrixToClusteredGraphs <- function(seuratObj_TCR = NULL,
                                             pcaComponents = 50,
                                             kpcaKernel = "rbfdot",
                                             usePCA = FALSE,
                                             partitionType = "RBConfigurationVertexPartition",
                                             proportionOfGraphAsNeighbors = 0.1,
                                             jaccardIndexThreshold = 0.1,
                                             resolutions = c(0.1, 0.2, 0.3),
                                             seed = 1234,
                                             computeMultiChain = T,
                                             verbose = FALSE,
                                             neighborsK = 15,
                                             rbfSigma = NULL,
                                             useExactDistanceKNN = FALSE,
                                             edgeWeightMode = "jaccard_x_rbf",
                                             useHdbscanTwoStage = FALSE,
                                             hdbscanAssays = c("TRB_cdr3"),
                                             minPtsRangeNoise = 2:20,
                                             minPtsRangeFidelity = 2:20) {

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
                                               # proportionOfGraphAsNeighbors kept for b/c; neighborsK takes precedence
                                               proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                               jaccardIndexThreshold = jaccardIndexThreshold,
                                               neighborsK = neighborsK,
                                               rbfSigma = rbfSigma,
                                               useExactDistanceKNN = useExactDistanceKNN,
                                               edgeWeightMode = edgeWeightMode)
    pruned_graph <- graph_and_pca_results$graph
    if (verbose) print(assay)
    single_chain_graphs[[assay]] <- pruned_graph
    for (resolution in resolutions) {
      partition <- leidenbase::leiden_find_partition(pruned_graph,
                                                     partition_type = partitionType,
                                                     initial_membership = NULL,
                                                     edge_weights = igraph::E(pruned_graph)$weight,
                                                     node_sizes = NULL,
                                                     seed = seed,
                                                     resolution_parameter = resolution,
                                                     num_iter = 10,
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
                                                  distanceMatrix = distanceMatrix,
                                                  neighborsK = neighborsK
    )
    #TODO: do a better job on integrating HDBSCAN two-stage into the main clustering flow.
    #consider adding exhausitive clustering methods. Run according to a vector of options?
    if (isTRUE(useHdbscanTwoStage) && assay %in% hdbscanAssays) {
      hdb_res <- .AddHdbscanTwoStageFiltering(
        seuratObj_TCR = seuratObj_TCR,
        assay = assay,
        layer = "counts",
        minPtsRangeNoise = minPtsRangeNoise,
        minPtsRangeFidelity = minPtsRangeFidelity,
        rbfSigma = rbfSigma,
        seed = seed,
        verbose = verbose
      )
      seuratObj_TCR <- hdb_res$seuratObj_TCR
    }
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
                                                   # proportionOfGraphAsNeighbors kept for b/c; neighborsK takes precedence
                                                   proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                                   jaccardIndexThreshold = jaccardIndexThreshold,
                                                   neighborsK = neighborsK,
                                                   rbfSigma = rbfSigma,
                                                   useExactDistanceKNN = useExactDistanceKNN,
                                                   edgeWeightMode = edgeWeightMode)
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
                                                                distanceMatrix = combined_matrix,
                                                                neighborsK = neighborsK
        )

        #TODO: add HDBSCAN two-stage here as well.
        for (resolution in resolutions) {
          partition <- leidenbase::leiden_find_partition(combined_graph,
                                                         partition_type = partitionType,
                                                         initial_membership = NULL,
                                                         edge_weights = igraph::E(combined_graph)$weight,
                                                         node_sizes = NULL,
                                                         seed = seed,
                                                         resolution_parameter = resolution,
                                                         num_iter = 10,
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
                              usePCA = FALSE,
                              proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                              jaccardIndexThreshold = jaccardIndexThreshold,
                              neighborsK = 15,
                              rbfSigma = NULL,
                              useExactDistanceKNN = FALSE,
                              edgeWeightMode = "jaccard_x_rbf"){

  #validate distance matrix dimensions
  if (any(dim(distanceMatrix) == 0)) {
    stop(paste("Distance matrix has zero dimensions:", paste(dim(distanceMatrix), collapse = "x")))
  }

  #validate number of components vs distance matrix dimensions
  if (pcaComponents > ncol(distanceMatrix)) {
   warning("Number of requested components exceeds available dimensions. Setting pcaComponents to the number of columns in the distance matrix.")
   pcaComponents <- ncol(distanceMatrix)
  }

  #determine k using fixed neighborsK (Scanpy/Seurat style). Fallback to proportion if provided.
  n <- nrow(distanceMatrix)
  k <- suppressWarnings(as.integer(neighborsK))
  if (is.na(k) || k < 1) {
    #fallback to proportion
    k <-  round(proportionOfGraphAsNeighbors * ncol(distanceMatrix))
  }
  if (k < 1) k <- 1
  if (k >= n) k <- max(1, n - 1)

  if (usePCA) {
    #use standard PCA
    #NOTE: applying PCA to a distance matrix is not recommended; consider usePCA = FALSE for KernelPCA from distances.
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
    #otherwise, kernel PCA from distances.
    #convert distances to an RBF kernel matrix; if rbfSigma is NULL, use median non-zero distance as heuristic
    d <- distanceMatrix
    #ensure symmetry
    if (!isTRUE(all.equal(d, t(d)))) {
      warning("Distance matrix is not symmetric; forcing symmetry by averaging with transpose.")
      d <- (d + t(d)) / 2
    }
    if (is.null(rbfSigma)) {
      #compute heuristic sigma
      upper <- d[upper.tri(d)]
      upper <- upper[is.finite(upper) & upper > 0]
      if (length(upper) == 0) {
        rbfSigma <- 1
      } else {
        rbfSigma <- stats::median(upper, na.rm = TRUE)
        if (!is.finite(rbfSigma) || rbfSigma <= 0) rbfSigma <- 1
      }
    }
    #compute kernel: K_ij = exp(-(d_ij^2)/(2*sigma^2))
    K <- exp(-(d^2) / (2 * (rbfSigma^2)))
    #construct kernel matrix for kernlab
    K <- kernlab::as.kernelMatrix(K)
    #kpca with precomputed kernel (kernel = "matrix")
    pca_result_obj <- kernlab::kpca(x = K, kernel = "matrix")
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

  #build kNN
  if (useExactDistanceKNN) {
    #build neighbors directly from the distance matrix
    #exclude self by removing the diagonal index
    nn.index <- matrix(NA_integer_, nrow = n, ncol = k)
    for (i in seq_len(n)) {
      #order by increasing distance; remove self
      ord <- order(distanceMatrix[i, ], decreasing = FALSE)
      ord <- ord[ord != i]
      nn.index[i, ] <- utils::head(ord, k)
    }
  } else {
    #use kNN on the KPCA-reduced data
    knn_result <- FNN::get.knn(reduced_data, k = k)
    nn.index <- knn_result$nn.index
  }

  #make initial kNN graph
  edges <- cbind(rep(seq_len(n), each = k), c(nn.index))
  g_knn <- igraph::graph_from_edgelist(edges, directed = FALSE)
  #remove loops
  g_knn <- igraph::simplify(g_knn)

  #compute SNN Jaccard weights explicitly and keep weights
  #NOTE: igraph::similarity returns a dense matrix; for large n this may be memory-intensive. TODO: switch to a sparse SNN builder if needed.
  jaccard_index <- igraph::similarity(g_knn, vids = igraph::V(g_knn), mode = 'all', method = 'jaccard')
  jaccard_index[!is.finite(jaccard_index)] <- 0
  diag(jaccard_index) <- 0
  #prune edges below the jaccard index threshold
  jaccard_index[jaccard_index < jaccardIndexThreshold] <- 0

  #distance-based weighting options
  adj <- jaccard_index
  if (edgeWeightMode != "jaccard") {
    d <- as.matrix(distanceMatrix)
    if (!isTRUE(all.equal(d, t(d)))) {
      d <- (d + t(d)) / 2
    }
    # derive sigma if needed
    sigma_use <- rbfSigma
    if (is.null(sigma_use) || !is.finite(sigma_use) || sigma_use <= 0) {
      upper <- d[upper.tri(d)]
      upper <- upper[is.finite(upper) & upper > 0]
      sigma_use <- if (length(upper) == 0) 1 else stats::median(upper, na.rm = TRUE)
      if (!is.finite(sigma_use) || sigma_use <= 0) sigma_use <- 1
    }
    sim <- exp(-(d^2) / (2 * (sigma_use^2)))
    diag(sim) <- 0
    if (edgeWeightMode == "rbf") {
      #restrict to SNN support after thresholding
      mask <- (jaccard_index > 0)
      adj <- sim * mask
    } else if (edgeWeightMode == "jaccard_x_rbf") {
      adj <- jaccard_index * sim
    }
  }

  #build weighted graph from chosen adjacency
  pruned_graph <- igraph::graph_from_adjacency_matrix(adj, mode = 'undirected', weighted = TRUE, diag = FALSE)

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
                                         distanceMatrix = NULL,
                                         neighborsK = 15) {
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
  #use fixed k for neighbors (Scanpy/Seurat style)
  k.param <-  as.integer(neighborsK)
  if (is.na(k.param) || k.param < 1) {
    #fallback to proportion-based if needed for b/c
    k.param <-  round(proportionOfGraphAsNeighbors * ncol(distanceMatrix))
  }
  k.param <- max(1, min(k.param, nrow(embeddings) - 1))

  #TODO: compute a cutoff for the number of components used, but I'm not sure what these distributions look like yet.
  #handle both PCA and kernel PCA for getting number of components
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
                               n.neighbors = max(15, min(c(nrow(embeddings), k.param))),
                               reduction.name = paste0(reductionName, "_umap"))
  return(seuratObj)
}

.AddHdbscanTwoStageFiltering <- function(seuratObj_TCR,
                                         assay = "TRB_cdr3",
                                         layer = "counts",
                                         minPtsRangeNoise = 2:20,
                                         minPtsRangeFidelity = 2:20,
                                         rbfSigma = NULL,
                                         seed = 1234,
                                         verbose = FALSE) {
  #get distance matrix with Seurat v5 fallbacks
  distance_matrix <- tryCatch({
    as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = layer))
  }, error = function(e1) {
    tryCatch({
      as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay))
    }, error = function(e2) {
      tryCatch({
        as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, slot = "data"))
      }, error = function(e3) {
        stop("Could not extract distance matrix from assay ", assay, ": ", e3$message)
      })
    })
  })

  if (!is.matrix(distance_matrix) || any(dim(distance_matrix) == 0)) {
    stop("Distance matrix is empty or not a matrix for assay: ", assay)
  }

  #ensure symmetry
  symmetric_distance_matrix <- as.matrix(distance_matrix)
  if (!isTRUE(all.equal(symmetric_distance_matrix, t(symmetric_distance_matrix)))) {
    if (verbose) message("Enforcing symmetry on distance matrix by averaging with transpose.")
    symmetric_distance_matrix <- (symmetric_distance_matrix + t(symmetric_distance_matrix)) / 2
  }

  # derive sigma if needed
  if (is.null(rbfSigma) || !is.finite(rbfSigma) || rbfSigma <= 0) {
    upper_triangle_values <- symmetric_distance_matrix[upper.tri(symmetric_distance_matrix)]
    upper_triangle_values <- upper_triangle_values[is.finite(upper_triangle_values) & upper_triangle_values > 0]
    rbfSigma <- if (length(upper_triangle_values) == 0) 1 else stats::median(upper_triangle_values, na.rm = TRUE)
    if (!is.finite(rbfSigma) || rbfSigma <= 0) rbfSigma <- 1
  }

  # RBF similarity kernel from distances
  kernel_matrix <- exp(-(symmetric_distance_matrix^2) / (2 * (rbfSigma^2)))
  kernel_matrix <- kernlab::as.kernelMatrix(kernel_matrix)

  #compute Kernel PCA features from the RBF kernel
  kpca_basis <- kernlab::kpca(x = kernel_matrix, kernel = "matrix")
  kpca_features <- kernlab::rotated(kpca_basis)
  if (is.null(rownames(kpca_features))) {
    if (!is.null(colnames(symmetric_distance_matrix))) {
      rownames(kpca_features) <- colnames(symmetric_distance_matrix)
    } else if (!is.null(rownames(symmetric_distance_matrix))) {
      rownames(kpca_features) <- rownames(symmetric_distance_matrix)
    } else {
      rownames(kpca_features) <- paste0("idx_", seq_len(nrow(kpca_features)))
    }
  }
  #map IDs once for both passes
  distance_ids <- if (!is.null(colnames(symmetric_distance_matrix))) colnames(symmetric_distance_matrix) else rownames(symmetric_distance_matrix)
  if (is.null(distance_ids)) distance_ids <- rownames(kpca_features)
  cells_all_distance_ids <- distance_ids

  #TODO: remove and depend on this.
  if (!requireNamespace("clusterCrit", quietly = TRUE)) {
    stop("Package 'clusterCrit' is required for silhouette optimization. Please install it.")
  }

  # helper function to compute silhouette on non-noise subset; returns NA when invalid
  .silhouette_for_labels <- function(labels, distance_mat) {
    non_noise_idx <- which(labels != 0)
    if (length(non_noise_idx) < 3) return(NA_real_)
    labels_non_noise <- labels[non_noise_idx]
    if (length(unique(labels_non_noise)) < 2) return(NA_real_)
    labels_reindexed <- as.integer(factor(labels_non_noise))
    crit <- try(clusterCrit::intCriteria(as.matrix(distance_mat[non_noise_idx, non_noise_idx, drop = FALSE]),
                                         labels_reindexed,
                                         c("Silhouette")), silent = TRUE)
    if (inherits(crit, "try-error") || is.null(crit$silhouette)) return(NA_real_)
    as.numeric(crit$silhouette)
  }
  #TODO: set seed higher up the stack
  set.seed(seed)

  #first pass: tune minPts for noise removal on full KPCA features
  if (verbose) message("Tuning HDBSCAN minPts for Stage 1 (noise removal) over ", paste(range(minPtsRangeNoise), collapse = ":"))
  silhouette_scores_noise <- sapply(minPtsRangeNoise, function(mp) {
    res <- dbscan::hdbscan(kpca_features, minPts = mp, cluster_selection_epsilon = 0)
    .silhouette_for_labels(res$cluster, symmetric_distance_matrix)
  })
  silhouette_scores_noise_clean <- ifelse(is.na(silhouette_scores_noise), -Inf, silhouette_scores_noise)
  selected_minPts_noise <- minPtsRangeNoise[which.max(silhouette_scores_noise_clean)]
  if (verbose) message("Selected minPts (Stage 1) = ", selected_minPts_noise)

  hdbscan_stage1 <- dbscan::hdbscan(kpca_features, minPts = selected_minPts_noise, cluster_selection_epsilon = 0)

  #optional heatmap for stage 1
  try({
    stage1_cluster_labels <- hdbscan_stage1$cluster
    unique_stage1_clusters <- sort(unique(stage1_cluster_labels[stage1_cluster_labels != 0]))
    cluster_palette_stage1 <- if (length(unique_stage1_clusters) > 0) tryCatch({ NatParksPalettes::natparks.pals("Charmonix", length(unique_stage1_clusters)) }, error = function(e) { grDevices::rainbow(length(unique_stage1_clusters)) }) else character(0)
    names(cluster_palette_stage1) <- as.character(unique_stage1_clusters)
    # include noise (0) color if present
    if (any(stage1_cluster_labels == 0)) {
      cluster_palette_stage1 <- c(`0` = "#BDBDBD", cluster_palette_stage1)
    }
    heatmap_stage1 <- tcrClustR:::.TCRDistanceHeatmap(
      seuratObj_TCR = seuratObj_TCR,
      assay = assay,
      cluster_info = stage1_cluster_labels,
      cluster_colors = cluster_palette_stage1
    )
    ComplexHeatmap::draw(heatmap_stage1)
  }, silent = TRUE)

  #write first-pass clusters to metadata for current assay cells only
  assay_cells <- colnames(seuratObj_TCR[[assay]])
  stage1_mapped_cells <- gsub("-", "_", cells_all_distance_ids)
  cells_to_write_stage1 <- intersect(stage1_mapped_cells, assay_cells)
  if (length(cells_to_write_stage1) == 0) {
    warning("No matching cell names found for Stage 1 between distance IDs and current assay cells; no first-pass HDBSCAN clusters were added.")
  } else {
    idx_in_all <- match(cells_to_write_stage1, stage1_mapped_cells)
    labels_to_write_stage1 <- hdbscan_stage1$cluster[idx_in_all]
    metadata_col_name_stage1 <- paste0("hdbscan_", assay, "_primary")
    metadata_to_add_stage1 <- data.frame(value = labels_to_write_stage1, row.names = cells_to_write_stage1)
    seuratObj_TCR <- Seurat::AddMetaData(seuratObj_TCR, metadata = metadata_to_add_stage1, col.name = metadata_col_name_stage1)
  }

  #filter non-noise points (for refinement) is replaced by rescue logic on noise-only points
  noise_indices <- which(hdbscan_stage1$cluster == 0)

  #if no noise, write rescue column identical to primary and exit
  if (length(noise_indices) == 0) {
    assay_cells <- colnames(seuratObj_TCR[[assay]])
    all_mapped_cells <- stage1_mapped_cells
    cells_to_write_rescue <- intersect(all_mapped_cells, assay_cells)
    if (length(cells_to_write_rescue) == 0) {
      warning("No matching cell names found between distance IDs and current assay cells; rescue column mirrors primary but nothing matched.")
    } else {
      idx_all <- match(cells_to_write_rescue, all_mapped_cells)
      rescue_labels_to_write <- hdbscan_stage1$cluster[idx_all]
      metadata_col_name_rescue <- paste0("hdbscan_", assay, "_rescue")
      metadata_to_add_rescue <- data.frame(value = rescue_labels_to_write, row.names = cells_to_write_rescue)
      seuratObj_TCR <- Seurat::AddMetaData(seuratObj_TCR, metadata = metadata_to_add_rescue, col.name = metadata_col_name_rescue)
    }

    return(list(
      seuratObj_TCR = seuratObj_TCR,
      kept_cells = character(0),
      stage1_minPts = selected_minPts_noise,
      stage2_minPts = NA_integer_,
      stage1_silhouette_scores = silhouette_scores_noise,
      stage2_silhouette_scores = NA
    ))
  }

  #rescue pass: run HDBSCAN only on Stage-1 noise
  kpca_features_noise <- kpca_features[noise_indices, , drop = FALSE]
  cells_noise_distance_ids <- distance_ids[noise_indices]

  if (verbose) message("Tuning HDBSCAN minPts for Rescue over ", paste(range(minPtsRangeFidelity), collapse = ":"))
  silhouette_scores_rescue <- sapply(minPtsRangeFidelity, function(mp) {
    res <- dbscan::hdbscan(kpca_features_noise, minPts = mp, cluster_selection_epsilon = 0)
    .silhouette_for_labels(res$cluster, symmetric_distance_matrix[noise_indices, noise_indices, drop = FALSE])
  })
  silhouette_scores_rescue_clean <- ifelse(is.na(silhouette_scores_rescue), -Inf, silhouette_scores_rescue)
  selected_minPts_rescue <- minPtsRangeFidelity[which.max(silhouette_scores_rescue_clean)]
  if (verbose) message("Selected minPts (Rescue) = ", selected_minPts_rescue)

  hdbscan_rescue <- dbscan::hdbscan(kpca_features_noise, minPts = selected_minPts_rescue, cluster_selection_epsilon = 0)
  rescue_labels_noise <- hdbscan_rescue$cluster

  #merge rescue labels back onto all cells (appending to the end of the 1st pass labels)
  merged_labels_all <- hdbscan_stage1$cluster
  max_primary <- ifelse(any(merged_labels_all > 0), max(merged_labels_all), 0)
  rescue_pos <- which(rescue_labels_noise > 0)
  if (length(rescue_pos) > 0) {
    merged_labels_all[noise_indices[rescue_pos]] <- rescue_labels_noise[rescue_pos] + max_primary
  }

  #map merged rescue labels back to metadata for current assay cells only
  assay_cells <- colnames(seuratObj_TCR[[assay]])
  all_mapped_cells <- stage1_mapped_cells
  cells_to_write_rescue <- intersect(all_mapped_cells, assay_cells)
  if (length(cells_to_write_rescue) == 0) {
    warning("No matching cell names found between distance IDs and current assay cells; no HDBSCAN rescue clusters were added.")
  } else {
    idx_all <- match(cells_to_write_rescue, all_mapped_cells)
    labels_to_write_rescue <- merged_labels_all[idx_all]
    metadata_col_name_rescue <- paste0("hdbscan_", assay, "_rescue")
    metadata_to_add_rescue <- data.frame(value = labels_to_write_rescue, row.names = cells_to_write_rescue)
    seuratObj_TCR <- Seurat::AddMetaData(seuratObj_TCR, metadata = metadata_to_add_rescue, col.name = metadata_col_name_rescue)
  }

  #optional heatmap for rescue subset only (1st pass's 'noise' cells)
  try({
    noise_mapped_cells <- gsub("-", "_", cells_noise_distance_ids)
    subset_cells_metadata <- intersect(noise_mapped_cells, assay_cells)
    subset_features <- gsub("_", "-", subset_cells_metadata)

    unique_rescue_clusters <- sort(unique(rescue_labels_noise[rescue_labels_noise != 0]))
    cluster_palette_rescue <- if (length(unique_rescue_clusters) > 0) tryCatch({ NatParksPalettes::natparks.pals("Charmonix", length(unique_rescue_clusters)) }, error = function(e) { grDevices::rainbow(length(unique_rescue_clusters)) }) else character(0)
    names(cluster_palette_rescue) <- as.character(unique_rescue_clusters)
    if (any(rescue_labels_noise == 0)) {
      cluster_palette_rescue <- c(`0` = "#BDBDBD", cluster_palette_rescue)
    }

    match_idx <- match(subset_cells_metadata, noise_mapped_cells)
    cluster_info_subset <- rescue_labels_noise[match_idx]

    heatmap_rescue <- tcrClustR:::.TCRDistanceHeatmap(
      seuratObj_TCR = subset(seuratObj_TCR, cells = subset_cells_metadata, features = subset_features),
      assay = assay,
      cluster_info = cluster_info_subset,
      cluster_colors = cluster_palette_rescue
    )
    ComplexHeatmap::draw(heatmap_rescue)
  }, silent = TRUE)

  return(list(
    seuratObj_TCR = seuratObj_TCR,
    kept_cells = gsub("-", "_", cells_noise_distance_ids),
    stage1_minPts = selected_minPts_noise,
    stage2_minPts = selected_minPts_rescue,
    stage1_silhouette_scores = silhouette_scores_noise,
    stage2_silhouette_scores = silhouette_scores_rescue
  ))
}


#' @title .DianaClustering
#' @description Internal function to perform DIANA (Divisive Analysis) hierarchical clustering
#' on a distance matrix and cut the dendrogram at a specified height.
#' @param distanceMatrix Distance matrix (as matrix or dist object)
#' @param cutHeight Height at which to cut the dendrogram. Default is 20.
#' @param verbose Boolean indicating whether to print verbose output. Default is FALSE.
#' @return A list containing the clustering assignments and dendrogram object
#' @keywords internal
.DianaClustering <- function(distanceMatrix, cutHeight = 20, verbose = FALSE) {
  if (verbose) {
    message("Running DIANA clustering with cutHeight = ", cutHeight)
  }
  
  #convert to dist object if necessary
  if (is.matrix(distanceMatrix)) {
    distObj <- stats::as.dist(distanceMatrix)
  } else {
    distObj <- distanceMatrix
  }
  
  #run DIANA clustering
  diana_result <- cluster::diana(distObj, diss = TRUE)
  
  #cut dendrogram at specified height
  cluster_assignments <- stats::cutree(stats::as.hclust(diana_result), h = cutHeight)
  
  if (verbose) {
    message("DIANA clustering produced ", length(unique(cluster_assignments)), " clusters")
  }
  
  return(list(
    clustering = cluster_assignments,
    diana_object = diana_result,
    cutHeight = cutHeight
  ))
}


#' @title .ThresholdClustersBySize
#' @description Internal function to filter clusters by minimum number of unique clones
#' @param clusterAssignments Named vector of cluster assignments (names are cell/clone IDs)
#' @param minClusterSize Minimum number of unique clones required to keep a cluster. Default is 2.
#' @param verbose Boolean indicating whether to print verbose output. Default is FALSE.
#' @return Named vector of cluster assignments with singletons/small clusters removed (set to 0)
#' @keywords internal
.ThresholdClustersBySize <- function(clusterAssignments, minClusterSize = 2, verbose = FALSE) {
  if (verbose) {
    message("Thresholding clusters with minimum size = ", minClusterSize)
  }
  
  #count cluster sizes
  cluster_sizes <- table(clusterAssignments)
  
  #identify clusters below threshold
  small_clusters <- names(cluster_sizes[cluster_sizes < minClusterSize])
  
  #set small clusters to 0 (noise/unassigned)
  filtered_assignments <- clusterAssignments
  filtered_assignments[clusterAssignments %in% small_clusters] <- 0
  
  if (verbose) {
    n_removed <- sum(clusterAssignments %in% small_clusters)
    message("Removed ", length(small_clusters), " clusters with < ", minClusterSize, " clones")
    message("Total clones removed: ", n_removed)
    message("Remaining clusters: ", length(unique(filtered_assignments[filtered_assignments != 0])))
  }
  
  return(filtered_assignments)
}


#' @title .WriteClusteringResultsToParquet
#' @description Internal function to write clustering results in columnar format to parquet file
#' @param metadata Data frame containing TCR metadata with columns for chains and gene segments
#' @param clusterAssignments Named vector of cluster assignments
#' @param assayName Name of the assay used for clustering
#' @param groupingVariable Name of the metadata column used for grouping (e.g., "tissue")
#' @param variableValue Value of the grouping variable for this batch (e.g., "lung")
#' @param clusterSizeThreshold Minimum cluster size threshold used
#' @param clusteringMethod Name of the clustering method used (e.g., "DIANA", "Leiden")
#' @param clusteringParams List of additional clustering parameters
#' @param outputPath Path to write the parquet file
#' @param chains Vector of chain names to include (e.g., c("TRA", "TRB"))
#' @param verbose Boolean indicating whether to print verbose output. Default is FALSE.
#' @return Path to the written parquet file
#' @keywords internal
.WriteClusteringResultsToParquet <- function(metadata,
                                              clusterAssignments,
                                              assayName,
                                              groupingVariable = NULL,
                                              variableValue = NULL,
                                              clusterSizeThreshold,
                                              clusteringMethod,
                                              clusteringParams = list(),
                                              outputPath,
                                              chains = c("TRA", "TRB"),
                                              verbose = FALSE) {
  if (verbose) {
    message("Writing clustering results to parquet: ", outputPath)
  }
  
  n_rows <- length(clusterAssignments)
  
  #determine which chains are actually used in this assay (for multi-chain analysis)
  #parse assay name to detect multi-chain combinations
  assay_normalized <- gsub("CDR3", "_cdr3", assayName)
  assay_parts <- strsplit(assay_normalized, "_")[[1]]
  
  #identify which chains are present in the assay name
  all_possible_chains <- c("TRA", "TRB", "TRG", "TRD")
  chains_in_assay <- all_possible_chains[sapply(all_possible_chains, function(chain) {
    any(grepl(chain, assay_parts, ignore.case = FALSE))
  })]
  
  #determine if this is a multi-chain assay (more than one chain type)
  is_multichain <- length(chains_in_assay) > 1
  
  if (is_multichain) {
    chains_used <- paste(sort(chains_in_assay), collapse = "+")
  } else if (length(chains_in_assay) == 1) {
    chains_used <- chains_in_assay[1]
  } else {
    #fallback if we can't parse the assay name
    chains_used <- "unknown"
  }
  
  if (verbose) {
    message("  Assay: ", assayName)
    message("  Detected chains: ", paste(chains_in_assay, collapse = ", "))
    message("  Multi-chain: ", is_multichain)
    message("  chains_used: ", chains_used)
  }
  
  #create results data frame with correct number of rows
  results_df <- data.frame(
    grouping_variable = rep(groupingVariable, n_rows),
    variable_value = rep(variableValue, n_rows),
    assay = rep(assayName, n_rows),
    chains_used = rep(chains_used, n_rows),
    stringsAsFactors = FALSE
  )
  
  #add ALL possible chain columns (TRA, TRB, TRG, TRD) with V, J, and CDR3
  #This ensures consistent schema across all outputs
  for (chain in all_possible_chains) {
    v_col <- paste0(chain, "_V")
    j_col <- paste0(chain, "_J")
    cdr3_col <- chain
    
    #check if columns exist in metadata AND if this chain is used in the assay
    #if chain is not used in assay, force to NA regardless of metadata
    chain_is_used <- chain %in% chains_in_assay
    
    if (chain_is_used && v_col %in% colnames(metadata)) {
      results_df[[v_col]] <- metadata[[v_col]][match(names(clusterAssignments), rownames(metadata))]
    } else {
      results_df[[v_col]] <- rep(NA_character_, n_rows)
    }
    
    if (chain_is_used && j_col %in% colnames(metadata)) {
      results_df[[j_col]] <- metadata[[j_col]][match(names(clusterAssignments), rownames(metadata))]
    } else {
      results_df[[j_col]] <- rep(NA_character_, n_rows)
    }
    
    if (chain_is_used && cdr3_col %in% colnames(metadata)) {
      results_df[[cdr3_col]] <- metadata[[cdr3_col]][match(names(clusterAssignments), rownames(metadata))]
    } else {
      results_df[[cdr3_col]] <- rep(NA_character_, n_rows)
    }
  }
  
  #add cluster assignments
  results_df$Cluster <- as.character(clusterAssignments)
  
  #add clustering parameters
  results_df$Cluster_Size_Threshold <- rep(clusterSizeThreshold, n_rows)
  results_df$Clustering_Method <- rep(clusteringMethod, n_rows)
  
  #add additional clustering parameters as JSON string or separate columns
  if (length(clusteringParams) > 0) {
    for (param_name in names(clusteringParams)) {
      results_df[[paste0("param_", param_name)]] <- rep(as.character(clusteringParams[[param_name]]), n_rows)
    }
  }
  
  #write to parquet
  arrow::write_parquet(results_df, outputPath)
  
  if (verbose) {
    message("Wrote ", nrow(results_df), " rows to ", outputPath)
  }
  
  return(outputPath)
}
