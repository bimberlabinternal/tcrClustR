


ClusterTcrs <- function(seuratObj = NULL,
                        seuratObj_TCR = NULL,
                        metadata = NULL,
                        resolutionParameter = 0.1,
                        kpcaComponents = 50,
                        kpcaKernel = "rbfdot",
                        proportionOfGraphAsNeighbors = 0.1,
                        jaccardIndexThreshold = 0.1,
                        seed = 1234,
                        spikeInDataframe =  NULL) {


  pruned_graphs <- .DistanceMatrixToGraphs(seuratObj_TCR = seuratObj_TCR,
                                           kpcaComponents = kpcaComponents,
                                           kpcaKernel = kpcaKernel,
                                           proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                           jaccardIndexThreshold = jaccardIndexThreshold)

  #iterate through resolutions and add them to the metadata


  #iterate through resolutions and graphs and add them to a dataframe
  resolutions_computed <- c()
  i = 1
  for (graph in names(pruned_graphs)) {
    #initialize the vectors to store whether or not the chain is only CDR3
    #or if the V/J segments should be used.
    cdr3_only_chains <- c()
    remaining_chains <- c()
    group_by_variables <- c()
    assays_to_access <- c()

    #the graphs are named with the chains (TRA if V+J+CDR3, or TRACDR3 if only CDR3) separated by underscores
    #(e.g. TRACDR3_TRB for only TRA's CDR3 and V+J+CDR3 for the beta chain)
    #parse these strings and assign them to how they should inform the dplyr grouping
    #to determine unique combinations of TRA+TRB observed in the data.
    #if (length(strsplit(graph, "_")) > 1 | any(grepl("CDR3", graph))) {
    #  cdr3_only_chains <- gsub("CDR3", "", strsplit(graph, "_")[grepl("CDR3", graph)][[1]])
    #} else if (!any(grepl("CDR3", graph))){
    #  remaining_chains <- strsplit(graph, "_")[!grepl("CDR3", graph)][[1]]
    #} else {
    #  remaining_chains <- graph
    #}


    if (length(strsplit(graph, "_")[[1]]) > 1 ) {
      chains <- strsplit(graph, "_")[[1]]
    } else {
      chains <- graph
    }
    print(chains)
    i = i + 1
    #chains <- c(cdr3_only_chains, remaining_chains)

    for (chain in chains) {
      print(chain)
      #populate the variables that will be used to group the metadata in the scRNASeq data
      #if (chain %in% cdr3_only_chains) {
        if (chain == "TRACDR3") {
          group_by_variables <- c(group_by_variables, "TRA")
          names(group_by_variables) <- c(names(group_by_variables)[names(group_by_variables) != ""], "TRA_cdr3")
          assays_to_access <- c(assays_to_access, "TRA_cdr3")
        } else if (chain == "TRBCDR3") {
          group_by_variables <- c(group_by_variables, "TRB")
          names(group_by_variables) <- c(names(group_by_variables)[names(group_by_variables) != ""], "TRB_cdr3")
          assays_to_access <- c(assays_to_access, "TRB_cdr3")
        } else if (chain == "TRGCDR3") {
          group_by_variables <- c(group_by_variables, "TRG")
          names(group_by_variables) <- c(names(group_by_variables)[names(group_by_variables) != ""], "TRG_cdr3")
          assays_to_access <- c(assays_to_access, "TRG_cdr3")
        } else if (chain == "TRDCDR3") {
          group_by_variables <- c(group_by_variables, "TRD")
          names(group_by_variables) <- c(names(group_by_variables)[names(group_by_variables) != ""], "TRD_cdr3")
          assays_to_access <- c(assays_to_access, "TRD_cdr3")
        } else if (chain == "TRA") {
          group_by_variables <- c(group_by_variables, "TRA_V", "TRA_J", "TRA")
          names(group_by_variables) <- c(names(group_by_variables)[names(group_by_variables) != ""], rep("TRA", 3))
          assays_to_access <- c(assays_to_access, "TRA")
        } else if (chain == "TRB") {
          group_by_variables <- c(group_by_variables, "TRB_V", "TRB_J", "TRB")
          names(group_by_variables) <- c(names(group_by_variables)[names(group_by_variables) != ""], rep("TRB", 3))
          assays_to_access <- c(assays_to_access, "TRB")
        } else if (chain == "TRG") {
          group_by_variables <- c(group_by_variables, "TRG_V", "TRG_J", "TRG")
          names(group_by_variables) <- c(names(group_by_variables)[names(group_by_variables) != ""], rep("TRG"), 3)
          assays_to_access <- c(assays_to_access, "TRG")
        } else if (chain == "TRD") {
          group_by_variables <- c(group_by_variables, "TRD_V", "TRD_J", "TRD")
          names(group_by_variables) <- c(names(group_by_variables)[names(group_by_variables) != ""], rep("TRD"), 3)
          assays_to_access <- c(assays_to_access, "TRD")
        }
      }
    print(paste0("group_variables:", group_by_variables))
    #iterate through the 10X data, merge with the spike-in dataframes and index the metadata by the:
    # 1. observed TRA+TRB combinations in the 10X data
    # 2. provided TRA+TRB combinations in the spikeInData
    # TODO: make this work with a metadata dataframe instead of only with a seurat object
    metadata <- plyr::rbind.fill(seuratObj@meta.data, spikeInDataframe)


    #if there are multiple chains, figure out how to combine them
    if (length(assays_to_access > 1)) {

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


    }
    #translate group_by_variables for compatibility with tcrdist3
    group_variables_tcrdist3 <- .TranslateGroupByVariablesToTcrdist3(group_by_variables)


  }
}

# pruned_graph <- graphs[[graph]]
#
# for (resolution in resolutionParameter) {
#   partition <- leidenbase::leiden_find_partition(pruned_graph,
#                                                  partition_type = c("CPMVertexPartition", "ModularityVertexPartition", "RBConfigurationVertexPartition", "RBERVertexPartition", "SignificanceVertexPartition", "SurpriseVertexPartition"),
#                                                  initial_membership = NULL,
#                                                  edge_weights = NULL,
#                                                  node_sizes = NULL,
#                                                  seed = seed,
#                                                  resolution_parameter = resolution,
#                                                  num_iter = 2,
#                                                  verbose = TRUE)
#
#   seuratObj_TCR <- Seurat::AddMetaData(seuratObj_TCR, partition$membership, col.name = paste0("TcrCluster_", graph, "_", resolution))
# }

#joiningScheme parsing

#we'll depend on using CellBarcode to resolve the TCR clustering using the v/j and/or CDR3s
#ensure these CellBarcodes are present.
# if (!("CellBarcode" %in% colnames(seuratObj@meta.data))) {
#   seuratMetadata <- seuratObj@meta.data |>
#     tibble::rownames_to_column(var = "CellBarcode")
# } else {
#   seuratMetadata <- seuratObj@meta.data
# }
# #if both TRA and TRB are present, use both chains to assign cells to clusters
# if ("cdr3_a_aa" %in% colnames(metadata) & "cdr3_b_aa" %in% colnames(metadata)) {
#   joinedMetadata <- metadata |>
#     dplyr::select(dplyr::all_of(c("cdr3_a_aa", "cdr3_b_aa", paste0("TcrClustR_",resolutions_computed)))) |>
#     dplyr::distinct()
#
#   if (!("CellBarcode" %in% colnames(seuratObj@meta.data))) {
#     seuratMetadata <- seuratObj@meta.data |>
#       tibble::rownames_to_column(var = "CellBarcode")
#   } else {
#     seuratMetadata <- seuratObj@meta.data
#   }
#
#   joinedMetadata <- dplyr::left_join(seuratMetadata, joinedMetadata, by = c("TRA" = "cdr3_a_aa", "TRB" = "cdr3_b_aa")) |>
#     dplyr::select(dplyr::all_of(c("CellBarcode", "TRA", "TRB", paste0("TcrCluster_",resolutions_computed)))) |>
#     dplyr::filter(!is.na(TRA) & !is.na(TRB)) |>
#     unique.data.frame()
#   rownames(joinedMetadata) <- joinedMetadata$CellBarcode
#
#   #TODO: figure out how to handle CDR3s that are assigned to multiple clusters
#   #This can happen when you compute distances with tcrdist3 using full TCR data (e.g. v and j segments)
#   #then use the cdr3s to assign clusters (of which there can be duplicates).
#   seuratObj <- Seurat::AddMetaData(seuratObj, joinedMetadata)
#
#   #if only one chain is present, use that chain to assign cells to clusters
# } else if ("cdr3_a_aa" %in% colnames(metadata) | "cdr3_b_aa" %in% colnames(metadata)){
#   #determine the chain present
#   if ("cdr3_a_aa" %in% colnames(metadata)) {
#     cdr3_variable <- "cdr3_a_aa"
#     chain <- "TRA"
#   } else if ("cdr3_b_aa" %in% colnames(metadata)) {
#     cdr3_variable <- "cdr3_b_aa"
#     chain <- "TRB"
#   } else {
#     stop("Error: unable to determine which chain is present in the metadata. Please ensure that either cdr3_a_aa or cdr3_b_aa is present in the metadata.")
#   }
#   #if only one chain is present, use that chain to assign cells to clusters
#   joinedMetadata <- metadata |>
#     dplyr::select(dplyr::all_of(c(chain, paste0("TcrCluster_",resolutions_computed)))) |>
#     dplyr::distinct()
#
#
#
# }




#' @title DistanceMatrixToGraph

.DistanceMatrixToGraphs <- function(seuratObj_TCR = NULL,
                       kpcaComponents = 50,
                       kpcaKernel = "rbfdot",
                       proportionOfGraphAsNeighbors = 0.1,
                       jaccardIndexThreshold = 0.1) {

  #check the Assays in the Seurat Object and compute graphs
  assays <- Seurat::Assays(seuratObj_TCR)

  combinations <- combn(assays, 2, simplify = FALSE)
  graphs <- list()

  #calculate the single chain graphs
  for (assay in assays){
    #get the distance matrix
    distanceMatrix <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = assay, layer = "counts"))
    pruned_graph <- .KpcaAndClustering(distanceMatrix = distanceMatrix,
                                kpcaComponents = kpcaComponents,
                                kpcaKernel = kpcaKernel,
                                proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                jaccardIndexThreshold = jaccardIndexThreshold)
    if (any(grepl("_cdr3",  assay))) {
      assay <- gsub("_cdr3", "CDR3", assay)
    } else {
      assay <- assay
    }
    graphs[[assay]] <- pruned_graph
  }

  #calculate all (non-self) joint chain graphs
  for (combo in combinations) {
    #ensure the first three letters of the assays are not identical
    #tcrdist3 computes the distance matrix for the full TCR and the CDR3 only, so we don't want to duplicate the CDR3's effect.
    if (substr(combo[1], 1, 3) == substr(combo[2], 1, 3)) {
      next
    }
    #get the distance matrix
    distanceMatrix_1 <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = combo[1], layer = "counts"))
    distanceMatrix_2 <- as.matrix(Seurat::GetAssayData(seuratObj_TCR, assay = combo[2], layer = "counts"))
    distanceMatrix <- as.matrix(distanceMatrix_1 + distanceMatrix_2)
    pruned_graph <- .KpcaAndClustering(distanceMatrix = distanceMatrix,
                                kpcaComponents = kpcaComponents,
                                kpcaKernel = kpcaKernel,
                                proportionOfGraphAsNeighbors = proportionOfGraphAsNeighbors,
                                jaccardIndexThreshold = jaccardIndexThreshold)
    #if the assay contains "_cdr3" in the name, this avoids names like "TRA_cdr3_TRB" in favor of "TRACDR3_TRB"
    #to easy splitting on underscores if one were to need to do so.
    if (any(grepl("_cdr3",  combo))) {
      combination_name <- gsub("_cdr3", "CDR3", combo)
    } else {
      combination_name <- combo
    }
    graphs[[paste(combination_name, collapse = "_")]] <- pruned_graph
  }
  return(graphs)
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
