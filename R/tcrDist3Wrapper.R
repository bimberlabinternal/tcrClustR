
utils::globalVariables(
  names = c('SubjectId', 'TRA_V', 'TRA_J', 'TRB_V', 'TRB_J', 'CloneNames', 'count'),
  package = 'tcrClustR',
  add = TRUE
)

#' @title RunTcrdist3
#' @description This function runs the tcrdist3 pipeline on a set of CDR3 sequences in a Seurat Object.
#'
#' @param inputData Either a seurat Object containing TCR information or a dataframe containing metadata
#' @param organism Organism to use for tcrdist3. Default is 'human'.
#' @param chains Vector of TCR chains to include in the analysis. Default is c("TRA", "TRB").
#' @param spikeInDataframe Data frame containing known CDR3s and gene segments to be included in the clustering. Default is NULL.
#' @param summarizeClones Pass-through boolean controlling whether to summarize clones into a frequency (by SubjectId, TRA, TRB, TRA_V, and TRB_J). Default is TRUE.
#' @param imputeCloneNames Pass-through boolean controlling whether to impute clone names if they are missing. Existing clone names will be inherited. Default is TRUE.
#' @param minimumClonesPerSubject Minimum number of clones per subject to include in the analysis. Default is 2.
#' @param multichain Boolean controlling whether to compute joint multi-chain distance matrices for observed chain combinations. Default is FALSE.
#' @param rdsOutputPath Path to the output directory for the RDS files containing the distance matrices. Default is "./tcrdist3DistanceMatrices/".
#' @param pythonExecutable Path to the python executable. Default is reticulate::py_exe().
#' @param debugTcrdist3 Boolean controlling whether to run tcrdist3 in debug mode. Default is TRUE.
#' @param verbose Boolean controlling whether to display processing steps. Default is FALSE.
#' @import Matrix
#' @importFrom methods as
#' @importFrom SeuratObject CreateSeuratObject Assays RenameAssays GetAssayData JoinLayers
#' @importFrom Seurat AddMetaData
#' @importFrom methods is
#' @examples
#' \dontrun{
#'   seuratObj_TCR<- RunTcrdist3(inputData = seuratObj@meta.data,
#'               chains = c("TRA", "TRB"),
#'               summarizeClones = T,
#'               imputeCloneNames = T,
#'               minimumClonesPerSubject = 2,
#'               multichain = FALSE,
#'               rdsOutputPath = "./tcrdist3DistanceMatrices/",
#'               pythonExecutable = reticulate::py_exe(),
#'               verbose = FALSE)
#'
#'   # Example with multichain analysis
#'   seuratObj_TCR<- RunTcrdist3(inputData = seuratObj,
#'               chains = c("TRA", "TRB", "TRG", "TRD"),
#'               multichain = TRUE,
#'               verbose = TRUE)
#'
#'     spikeInDataframe <- data.frame(CloneNames = rep(1:3),
#'                                  TRA_V = c("TRAV1-2", "TRAV1-2", "TRAV1-2"),
#'                                  TRA_J = c("TRAJ33", "TRAJ20", "TRAJ33"),
#'                                  TRA = c("CAVRDSNYQLIW", "CAVSLQDYKLSF", "CAVRDSNYQLIW"),
#'                                  TRB_V = c("TRBV6-4", "TRBV6-4", "TRBV6-4"),
#'                                  TRB_J = c("TRBJ1-1", "TRBJ2-1", "TRBJ2-3"),
#'                                  TRB = c("CASSAAAAAAAAFF", "CASSVVVVVVVVQF", "CASSWWWWWWWWQY"))
#' }
#'
#' @export

#TODO: flesh out examples demonstrating formatting requirements for spikeInDataframe


RunTcrdist3 <- function(inputData = NULL,
                        organism = 'human',
                        chains = c("TRA", "TRB"),
                        spikeInDataframe = NULL,
                        summarizeClones = T,
                        imputeCloneNames = T,
                        minimumClonesPerSubject = 2,
                        multichain = FALSE,
                        rdsOutputPath = "./tcrdist3DistanceMatrices/",
                        pythonExecutable = NULL,
                        debugTcrdist3 = TRUE,
                        verbose = FALSE) {

  if (is.null(inputData)) {
    stop("Please provide either a Seurat Object or the Seurat Object's metadata as input.")
  }

  if (typeof(inputData) == 'S4' && class(inputData)[1] == 'Seurat' ){
    metadata <- inputData@meta.data
  }

  if (!is.data.frame(metadata)) {
    stop('Expected inputData to be either a Seurat object or a data.frame')
  }

  if (is.null(pythonExecutable)) {
    pythonExecutable <- reticulate::py_exe()

    #if reticulate doesn't provide a valid executable, try system default
    if (is.null(pythonExecutable) || pythonExecutable == "" || !file.exists(pythonExecutable)) {
      pythonExecutable <- Sys.which("python3")
      if (pythonExecutable == "") {
        pythonExecutable <- Sys.which("python")
      }
    }

    if (pythonExecutable == "" || !file.exists(pythonExecutable)) {
      stop("No valid Python executable found. Please specify pythonExecutable parameter or ensure Python is in PATH.")
    }
  }

  #validate Python environment and required packages
  if (verbose) {
    print(paste("Using Python executable:", pythonExecutable))
  }

  # NOTE: The package is installed as 'tcrdist3' but imported as 'tcrdist'
  #test if required packages are available
  testCmd <- paste0(pythonExecutable, " -c \"import tcrdist; import rpy2; print('Python environment OK')\"")
  testResult <- tryCatch({
    system(testCmd, intern = TRUE)
  }, error = function(e) {
    warning(paste("Python environment test failed:", e$message))
  })

  if (all(is.null(testResult))) {
    stop("Python environment may not have required packages (tcrdist3, rpy2). This may cause failures.")
  }

  postFormattingMetadataCsvOutput <- tempfile(fileext = ".csv")
  postFormattingMetadataCsvOutput <- gsub(postFormattingMetadataCsvOutput, pattern = '\\\\', replacement = '/')

  rdsOutputPath <- R.utils::getAbsolutePath(rdsOutputPath)
  rdsOutputPath <- gsub(rdsOutputPath, pattern = '\\\\', replacement = '/')

  pythonExecutable <- R.utils::getAbsolutePath(pythonExecutable)

  formatted_metadata <- FormatMetadataForTcrDist3(metadata = metadata,
                                        outputPath = dirname(postFormattingMetadataCsvOutput),
                                        outputCsvName = basename(postFormattingMetadataCsvOutput),
                                        chains = chains,
                                        summarizeClones = summarizeClones,
                                        imputeCloneNames = imputeCloneNames,
                                        minimumClonesPerSubject = minimumClonesPerSubject,
                                        spikeInDataframe = spikeInDataframe,
                                        verbose = verbose
  )

  #convert chains into a string for parsing in python
  chainsString <- ""
  for (chain in chains) {
    if (chain == "TRA") {
      chainsString <- paste0(chainsString, "alpha")
      } else if (chain == "TRB") {
        chainsString <- paste0(chainsString, "beta")
      } else if (chain == "TRG") {
        chainsString <- paste0(chainsString, "gamma")
      } else if (chain == "TRD") {
        chainsString <- paste0(chainsString, "delta")
      } else {
        warning(paste0("Chain Type ", chain, " not recognized. Skipping."))
      }
    }

  #run tcrdist3 in python and return individual RDS files with the results (they may get very large)

  #read the python script template from tcrClustR and write them to a tempfile
  template <- readr::read_file(system.file("scripts/tcrdist3TcrDistances.py", package = "tcrClustR"))
  script <- tempfile(fileext = ".py")
  readr::write_file(template, script)

  #create distance matrix output directory if it doesn't exist
  if (!dir.exists(rdsOutputPath)) {
    dir.create(rdsOutputPath)
  }
  if (verbose) message("Creating tcrdist3 distance matrices in the following directory: ", rdsOutputPath)
  #format and write the python function to the end of the script
  command <- paste0("writeTcrDistances(csv_path = '", postFormattingMetadataCsvOutput,
                   "', organism = '", organism,
                   "', chainsString = '", chainsString,
                    #TODO: review what DB this is using, and probably conditionalize based on organisms
                   "', db_file = 'alphabeta_gammadelta_db.tsv",
                   "', rds_output_path = '", rdsOutputPath,
                   "', debug ='", ifelse(!is.null(debugTcrdist3) && debugTcrdist3, yes = 'True', no = 'False'),
                   "')")
  readr::write_file(command, script, append = TRUE)
  system2(pythonExecutable, script)

  #return a seurat object, with the distance matrices implemented as assays
  seuratObj_TCR <- NULL

  for (chain in chains) {
    #read the RDS file
    if (chain == "TRA") {
      chain_tcrdist3 <- "alpha"
    } else if (chain == "TRB") {
      chain_tcrdist3 <- "beta"
    } else if (chain == "TRG") {
      chain_tcrdist3 <- "gamma"
    } else if (chain == "TRD") {
      chain_tcrdist3 <- "delta"
    } else {
      warning(paste0("Chain Type ", chain, " not recognized. Skipping."))
      next
    }
    #process the full length TCR distance matrix
    rdsFile <- paste0(rdsOutputPath, "/pw_", chain_tcrdist3, ".rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise 'full' tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_full_length <- readRDS(rdsFile)
    colnames(distanceMatrix_full_length) <- paste0(chain, "_", seq_len(ncol(distanceMatrix_full_length)))
    rownames(distanceMatrix_full_length) <- paste0(chain, "_", seq_len(nrow(distanceMatrix_full_length)))

    #process the CDR3 only TCR distance matrix
    #grab the first letter of the chain (distance matrices for the cdr3 are stored as "pw_cdr3_b_aa.rds" for a beta chain)
    chain_cdr3_id <- tolower(strsplit(chain_tcrdist3, split = "")[[1]][1])
    rdsFile <- paste0(rdsOutputPath, "/pw_cdr3_", chain_cdr3_id, "_aa.rds")
    if (!file.exists(rdsFile)) {
      stop(paste0("Pairwise CDR3 tcrdist3 distance matrix RDS file not found: ", rdsFile))
    }
    distanceMatrix_CDR3 <- readRDS(rdsFile)
    colnames(distanceMatrix_CDR3) <- paste0(chain, "_", seq_len(ncol(distanceMatrix_CDR3)), "_cdr3")
    rownames(distanceMatrix_CDR3) <- paste0(chain, "_", seq_len(nrow(distanceMatrix_CDR3)), "_cdr3")

    #add the distance matrices to the Seurat object
    if (is.null(seuratObj_TCR)){
      formatted_metadata_for_chain <- FormatMetadataForTcrDist3(metadata = metadata,
                                                      chains = chain,
                                                      organism = organism,
                                                      summarizeClones = summarizeClones,
                                                      imputeCloneNames = imputeCloneNames,
                                                      minimumClonesPerSubject = minimumClonesPerSubject,
                                                      spikeInDataframe = spikeInDataframe,
                                                      verbose = FALSE
      )

      seuratObj_TCR <- SeuratObject::CreateSeuratObject(counts = as(distanceMatrix_full_length, "dgCMatrix"), assay =  chain, metadata = formatted_metadata_for_chain)
      seuratObj_TCR_CDR3 <- SeuratObject::CreateSeuratObject(counts = as(distanceMatrix_CDR3, "dgCMatrix"), assay = paste0(chain, "_cdr3"))
      seuratObj_TCR_CDR3 <- Seurat::AddMetaData(seuratObj_TCR_CDR3, metadata = formatted_metadata_for_chain)
      seuratObj_TCR <- merge(seuratObj_TCR, seuratObj_TCR_CDR3)
    } else {
      print('Adding!')
      seuratObj_TCR_subsequentChain <- SeuratObject::CreateSeuratObject(counts = as(distanceMatrix_full_length, "dgCMatrix"), assay =  chain)
      seuratObj_TCR_subsequentChain <- Seurat::AddMetaData(seuratObj_TCR_subsequentChain, metadata = formatted_metadata_for_chain)
      seuratObj_TCR_CDR3_subsequentChain <- SeuratObject::CreateSeuratObject(counts = as(distanceMatrix_CDR3, "dgCMatrix"), assay = paste0(chain, "_cdr3"))
      seuratObj_TCR_CDR3_subsequentChain <- Seurat::AddMetaData(seuratObj_TCR_CDR3_subsequentChain, metadata = formatted_metadata_for_chain)
      seuratObj_TCR_subsequentChain <- merge(seuratObj_TCR_subsequentChain, seuratObj_TCR_CDR3_subsequentChain)
      seuratObj_TCR <- merge(seuratObj_TCR, seuratObj_TCR_subsequentChain)
    }
  }

  #add multichain tcrdist matrices if requested
  if (multichain && length(chains) > 1) {
    if (verbose) message("Computing multichain distance matrices...")

    #create joint distance matrices for observed chain combinations
    #Note: matrices represent distances between observations (rows in formatted metadata)
    #that have information for both chains in the combination
    chain_combinations <- utils::combn(chains, 2, simplify = FALSE)

    for (combo in chain_combinations) {
      chain1 <- combo[1]
      chain2 <- combo[2]

      #skip identical chains (trivial case)
      if (chain1 == chain2) {
        if (verbose) message(paste("Skipping", chain1, "+", chain2, "combination - identical chains"))
        next
      }

      #sanity check to make sure the assays exist before pulling
      if (!chain1 %in% SeuratObject::Assays(seuratObj_TCR) || !chain2 %in% SeuratObject::Assays(seuratObj_TCR)) {
        if (verbose) message(paste("Skipping", chain1, "+", chain2, "combination - matrices not found"))
        next
      }

      #get distance matrices for both chains
      dist_matrix1 <- SeuratObject::GetAssayData(seuratObj_TCR, assay = chain1, layer = "counts")
      dist_matrix2 <- SeuratObject::GetAssayData(seuratObj_TCR, assay = chain2, layer = "counts")

      if (verbose) {
        message(paste("Matrix dimensions for", chain1, ":", nrow(dist_matrix1), "x", ncol(dist_matrix1)))
        message(paste("Matrix dimensions for", chain2, ":", nrow(dist_matrix2), "x", ncol(dist_matrix2)))
      }

      # TODO: GW, can we just drop this entirely?
      #for different-sized matrices (almost certainly the case), identify observations with both chains
      #read the formatted metadata to map sequence indices to observations
      if (!exists("formatted_metadata")) {
        if (file.exists(postFormattingMetadataCsvOutput)) {
          formatted_metadata <- readr::read_csv(postFormattingMetadataCsvOutput, show_col_types = FALSE)
        } else {
          if (verbose) message("Warning: Cannot create multichain matrices - formatted metadata file not found")
          next
        }
      }

      #whole-dataset chain pairing vectors
      # Map chain names to tcrdist3 column names
      chain1_col <- NULL
      chain2_col <- NULL

      if (chain1 == "TRA") {
        chain1_col <- "cdr3_a_aa"
      } else if (chain1 == "TRB") {
        chain1_col <- "cdr3_b_aa"
      } else if (chain1 == "TRG") {
        chain1_col <- "cdr3_g_aa"
      } else if (chain1 == "TRD") {
        chain1_col <- "cdr3_d_aa"
      }

      if (chain2 == "TRA") {
        chain2_col <- "cdr3_a_aa"
      } else if (chain2 == "TRB") {
        chain2_col <- "cdr3_b_aa"
      } else if (chain2 == "TRG") {
        chain2_col <- "cdr3_g_aa"
      } else if (chain2 == "TRD") {
        chain2_col <- "cdr3_d_aa"
      }

      if (is.null(chain1_col) || is.null(chain2_col)) {
        if (verbose) message(paste("Skipping", chain1, "+", chain2, "combination - unsupported chain types"))
        next
      }

      if (!chain1_col %in% colnames(formatted_metadata) ||
          !chain2_col %in% colnames(formatted_metadata)) {
        if (verbose) message(paste("Skipping", chain1, "+", chain2, "combination - chain columns not found in metadata"))
        if (verbose) message(paste0("Looking for columns: ", chain1_col, ", ", chain2_col))
        if (verbose) message(paste0("Available columns: ", paste(colnames(formatted_metadata), collapse = ", ")))
        next
      }

      #find rows (cells) that have non-NA values for both chains
      both_chains_present <- !is.na(formatted_metadata[[chain1_col]]) &
                            !is.na(formatted_metadata[[chain2_col]]) &
                            formatted_metadata[[chain1_col]] != "" &
                            formatted_metadata[[chain2_col]] != ""

      if (sum(both_chains_present) == 0) {
        if (verbose) message(paste("Skipping", chain1, "+", chain2, "combination - no observations with both chains"))
        next
      }

      #get the subset of metadata with both chains and set the dimensions of the joint matrix
      dual_chain_metadata <- formatted_metadata[both_chains_present, ]
      n_observations <- nrow(dual_chain_metadata)

      if (verbose) message(paste("Found", n_observations, "observations with both", chain1, "and", chain2))

      #create unique sequence mappings for each chain
      #each chain's distance matrix rows/cols correspond to unique sequences in that chain
      chain1_sequences <- unique(dual_chain_metadata[[chain1_col]])
      chain2_sequences <- unique(dual_chain_metadata[[chain2_col]])

      #create mapping from sequences to matrix indices
      chain1_seq_to_idx <- setNames(seq_along(chain1_sequences), chain1_sequences)
      chain2_seq_to_idx <- setNames(seq_along(chain2_sequences), chain2_sequences)

      #create all possible matrix combinations (full + CDR3) for this chain pair
      #this allows for mixed comparisons like full TRA + CDR3 TRB

      #init available matrix types for each chain
      chain1_matrices <- list()
      chain2_matrices <- list()

      #add full-length matrices
      chain1_matrices[[chain1]] <- dist_matrix1
      chain2_matrices[[chain2]] <- dist_matrix2

      #add CDR3 matrices (if they exist)
      cdr3_assay1 <- paste0(chain1, "_cdr3")
      cdr3_assay2 <- paste0(chain2, "_cdr3")

      if (cdr3_assay1 %in% SeuratObject::Assays(seuratObj_TCR)) {
        chain1_matrices[[paste0(chain1, "_cdr3")]] <- SeuratObject::GetAssayData(seuratObj_TCR, assay = cdr3_assay1, layer = "counts")
      }

      if (cdr3_assay2 %in% SeuratObject::Assays(seuratObj_TCR)) {
        chain2_matrices[[paste0(chain2, "_cdr3")]] <- SeuratObject::GetAssayData(seuratObj_TCR, assay = cdr3_assay2, layer = "counts")
      }

      #create joint matrices for all combinations
      for (matrix1_name in names(chain1_matrices)) {
        for (matrix2_name in names(chain2_matrices)) {

          matrix1 <- chain1_matrices[[matrix1_name]]
          matrix2 <- chain2_matrices[[matrix2_name]]

          #initialize joint distance matrix (observations x observations)
          joint_matrix <- matrix(0, nrow = n_observations, ncol = n_observations)

          # for each pair of observations, compute joint distance
          for (i in 1:n_observations) {
            for (j in 1:n_observations) {
              #get sequences for observation i
              seq1_chain1 <- dual_chain_metadata[[chain1_col]][i]
              seq1_chain2 <- dual_chain_metadata[[chain2_col]][i]

              #get sequences for observation j
              seq2_chain1 <- dual_chain_metadata[[chain1_col]][j]
              seq2_chain2 <- dual_chain_metadata[[chain2_col]][j]

              #get matrix indices for these sequences
              idx1_chain1 <- chain1_seq_to_idx[seq1_chain1]
              idx2_chain1 <- chain1_seq_to_idx[seq2_chain1]
              idx1_chain2 <- chain2_seq_to_idx[seq1_chain2]
              idx2_chain2 <- chain2_seq_to_idx[seq2_chain2]

              #extract distances from matrices
              dist_chain1 <- 0
              dist_chain2 <- 0

              if (!is.na(idx1_chain1) && !is.na(idx2_chain1) &&
                  idx1_chain1 <= nrow(matrix1) && idx2_chain1 <= ncol(matrix1)) {
                dist_chain1 <- matrix1[idx1_chain1, idx2_chain1]
              }

              if (!is.na(idx1_chain2) && !is.na(idx2_chain2) &&
                  idx1_chain2 <= nrow(matrix2) && idx2_chain2 <= ncol(matrix2)) {
                dist_chain2 <- matrix2[idx1_chain2, idx2_chain2]
              }

              #fill in joint matrix
              joint_matrix[i, j] <- dist_chain1 + dist_chain2
            }
          }

          #create assay name based on matrix types
          joint_assay_name <- paste0(matrix1_name, "_", matrix2_name)
          colnames(joint_matrix) <- paste0(joint_assay_name, "_", seq_len(ncol(joint_matrix)))
          rownames(joint_matrix) <- paste0(joint_assay_name, "_", seq_len(nrow(joint_matrix)))

          #create separate seurat object for joint matrix
          seuratObj_joint <- SeuratObject::CreateSeuratObject(counts = as(joint_matrix, "dgCMatrix"),
                                                              assay = joint_assay_name)
          seuratObj_joint <- Seurat::AddMetaData(seuratObj_joint, metadata = dual_chain_metadata)

          #merge with main object
          seuratObj_TCR <- merge(seuratObj_TCR, seuratObj_joint)

          if (verbose) message(paste("Created joint distance matrix:", joint_assay_name))
        }
      }
    }
  }

  unlink(postFormattingMetadataCsvOutput)

  return(seuratObj_TCR)
}
