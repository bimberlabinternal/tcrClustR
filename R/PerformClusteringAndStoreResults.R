#' @title PerformClusteringAndStoreResults
#' @description A wrapper method designed to provide a start-to-finish pipeline to cluster/store TCR families for a given seurat object
#'
#' @param seuratObj The Seurat object
#' @param outputDir Directory to write parquet files. Default = "./tcrclustering_output/".
#' @param outputPrefix Optional prefix for output filenames. Default = "clustering".
#' @param verbose Logical. Print progress. Default = TRUE.
#' @return A seurat object with clustering information stored as assays
#' @export
PerformClusteringAndStoreResults <- function(seuratObj, outputDir, outputPrefix = 'TCR.clustering', chains, verbose = TRUE, minimumCloneSize) {
  rdsOutputPath <- file.path(outputDir, ".tcrdist3DistanceMatrices")

  distDistOutput <- RunTcrdist3(
    inputData = seuratObj,
    calculateChainPairs = TRUE,
    chains = chains,
    minimumCloneSize = 2,
    rdsOutputPath = rdsOutputPath,
    verbose = verbose
  )
   
  print("TCR distance matrices computed successfully!")
  print(paste0("Available assays: ", paste(SeuratObject::Assays(seuratObj_TCR), collapse = ", ")))
  print(paste0("Number of cells in TCR object: ", ncol(seuratObj_TCR)))

  clustering_result <- RunTcrClustering(
    seuratObj_TCR = seuratObj_TCR,
    assaysToCluster = "TRB",
    resolutionParameters = resolutionParameters,
    usePCA = usePCA,
    dianaHeight = 20,
    clusterSizeThreshold = 1,
    clusteringMethod = "DIANA",
    pcaComponents = pcaComponents,
    outputDir = file.path(outputPath, "tcrclustering_output"),
    pythonExecutable = pythonExecutable,
    verbose = verbose,
    stripAlleles = TRUE
  )

  print("TCR clustering completed successfully!")

  #By default, RunTcrClustering writes parquet files to "./tcrclustering_output/".
  #If you set a custom outputDir in RunTcrClustering(outputDir = "./my_output/"), use that here.
  seurat_with_families <- JoinClusteringResults(
    seuratObj = seuratObj,
    parquetDir = file.path(outputPath, "tcrclustering_output"),
    filePattern = "\\.parquet$",
    metadataColumnPrefix = "TcrFamily",
    stripAlleles = TRUE,
    overwriteExisting = TRUE,
    verbose = TRUE
  )

  #check that you have a UMAP-able number of classes:
  unique_classes <- length(levels(seurat_with_families@meta.data[[new_cols[1]]]))
  print(paste0("Number of unique TCR families (including LowFrequency and No_TCR_Data): ", unique_classes))

  DimPlot(
    seurat_with_families,
    reduction = "umap",
    group.by = new_cols,
    label = TRUE,
    pt.size = 1
  )

  return(seuratObj)
}