# tcrClustR

An R package for clustering and analyzing T-cell receptor (TCR) sequences to identify 'TCR families' via sequence similarity. This package uses [tcrdist3](https://github.com/kmayerb/tcrdist3) for TCR distance calculations and provides flexible clustering algorithms to identify groups of functionally related TCRs. Similar in concept to GLIPH and CoNGA, but with more direct control over clustering parameters.

Full documentation: https://bimberlabinternal.github.io/tcrClustR/

## Table of Contents
* [Quick Start](#quick-start)
* [Overview](#overview)
* [Installation](#installation)
* [Usage Examples](#usage-examples)
* [Workflows](#workflows)
* [Output Schemas](#output-file-formats)
* [Known Issues](#issues)
* [Development Guidelines](#developers)

## Quick Start

The fastest way to cluster TCR data:

```r
library(tcrClustR)

#Step 1: compute TCR distance matrices (stored as Seurat assays)
seuratObj_TCR <- RunTcrdist3(
  inputData = seuratObj,
  chains = c("TRA", "TRB"),
  minimumCloneSize = 2
)

#Step 2: cluster and export results to parquet files
results <- RunTcrClustering(
  seuratObj_TCR = seuratObj_TCR,
  clusteringMethod = "DIANA",
  outputDir = "./tcrclustering_output/"
)

#Step 3: join parquet clustering results back to the original Seurat object
seurat_with_families <- JoinClusteringResults(
  seuratObj = seuratObj,
  parquetFiles = results$parquet_files,
  metadataColumnPrefix = "TcrFamily"
)

#Parquet outputs contain cluster assignments for single and paired assays.
#Single-chain: chain, v_gene, j_gene, CDR3, Cluster, ...
#Paired-chain: chain_1/2, v_gene_1/2, j_gene_1/2, CDR3_1/2, Cluster, ...
```

## Overview

### What is TCR Clustering?

T-cell receptor (TCR) clustering groups TCR sequences based on similarity metrics, enabling identification of functionally related TCR families from single-cell sequencing data. TCRs with similar sequences (CDR3 regions and V/J gene segments) may recognize the same or related antigens.

### Workflow

1. **Format & Validate**: Clean TCR metadata, filter low-quality clones
2. **Compute Distances**: Calculate pairwise TCR distances using tcrdist3 (BLOSUM62 matrix)
3. **Cluster**: Apply network-based (Leiden) or hierarchical (DIANA) clustering
4. **Export**: Generate parquet tables with cluster assignments for downstream analysis
5. **Visualize**: Create heatmaps and plots to explore clustering results

### Key Features

- **Flexible clustering**: DIANA (hierarchical) or Leiden (network-based) algorithms
- **Single & paired chain analysis**: TRA, TRB, or combined TRA+TRB distances
- **Automatic data filtering**: Remove NA and concatenated values (optional)
- **Allele normalization**: Strip allele suffixes like `*01` for cleaner analysis
- **Compact outputs**: Parquet files optimized for large-scale analysis

## Installation

### From GitHub (Recommended)

```r
# Install devtools if needed
if (!require("devtools")) install.packages("devtools")

# Install tcrClustR
devtools::install_github("bimberlabinternal/tcrClustR")
```

### Docker (Includes Python Dependencies)

```bash
docker pull ghcr.io/bimberlabinternal/tcrclustr:latest
```

### Python Dependencies

tcrClustR requires Python 3.8+ with tcrdist3 and related packages. The package includes tools to simplify setup.

#### Automated Setup (Recommended)

Use the built-in helper function to validate and install Python dependencies:

```r
library(tcrClustR)

#check and install Python dependencies automatically
SetupPythonEnvironment()

#or just validate without installing
SetupPythonEnvironment(installMissing = FALSE)

#use specific Python executable
SetupPythonEnvironment(pythonExecutable = "/path/to/python3")
```

This function:
- Validates Python installation (requires 3.8+)
- Checks for required modules (tcrdist3, pandas, numpy, rpy2)
- Installs missing packages from `requirements.txt`


#### Manual Setup

If you prefer manual installation:

```bash
# install individual packages
pip install pandas numpy scikit-learn rpy2
pip install git+https://github.com/kmayerb/tcrdist3.git@0.2.2

#optional: install from requirements.txt in this repo
pip install -r requirements.txt
```

Set the Python path in R if needed:

```r
Sys.setenv(RETICULATE_PYTHON = "/path/to/python3")
```

#### Troubleshooting Python Issues

If you encounter Python-related errors:

1. **Run the setup helper**: `SetupPythonEnvironment(verbose = TRUE)`
2. **Check Python version**: Must be 3.8 or higher
3. **Verify tcrdist3**: `python3 -c 'import tcrdist; print(tcrdist.__version__)'`
4. **Check reticulate config**: `reticulate::py_config()`
5. **Review error logs**: The package now captures and displays Python stderr/stdout

Common error messages and solutions:

```r
# Error: "Missing required Python modules: tcrdist"
# Solution: Run SetupPythonEnvironment() to install

# Error: "No valid Python executable found"
# Solution: Install Python 3.8+ or specify path:
SetupPythonEnvironment(pythonExecutable = "/usr/bin/python3")

# Error: "tcrdist3_gene_segments.txt generation failed"
# Solution: Check Python dependencies with verbose mode:
FormatMetadataForTcrDist3(..., verbose = TRUE)
```

## Usage Examples

### Example 1: Primary Workflow (RunTcrClustering -> JoinClusteringResults)

**`RunTcrClustering()`** is the recommended function for most users who have already subset their data:

```r
library(tcrClustR)

#Compute distances
seuratObj_TCR <- RunTcrdist3(
  inputData = seuratObj,
  chains = c("TRA", "TRB"),
  minimumCloneSize = 2
)

#Cluster with automatic filtering and allele stripping
results <- RunTcrClustering(
  seuratObj_TCR = seuratObj_TCR,
  clusteringMethod = "DIANA",       # or "Leiden"
  dianaHeight = 20,                  # cut height for DIANA
  clusterSizeThreshold = 2,          # min cluster size
  filterInvalidClones = TRUE,        # remove NA/concatenated (default)
  stripAlleles = TRUE,               # remove *01 suffixes (default)
  outputDir = "./clustering_output/"
)

#Join results back to the original Seurat object
seurat_with_families <- JoinClusteringResults(
  seuratObj = seuratObj,
  parquetFiles = results$parquet_files,
  metadataColumnPrefix = "TcrFamily",
  stripAlleles = TRUE
)

#Read parquet directly if desired
library(arrow)
clusters <- read_parquet(results$parquet_files[1])
head(clusters)
```

### Example 2: Advanced Clustering with ClusterTcrs (Seurat-integrated)

For interactive analysis with full Seurat integration:

```r
#Compute distances
seuratObj_TCR <- RunTcrdist3(
  inputData = seuratObj,
  chains = c("TRA", "TRB"),
  minimumCloneSize = 2
)

#Cluster with multiple resolutions
clustered_results <- ClusterTcrs(
  seuratObj_TCR = seuratObj_TCR,
  resolutionParameters = c(0.1, 0.2, 0.5),
  usePCA = TRUE,
  pcaComponents = 50
)

# Visualize
TCRDistanceHeatmaps(
  seuratObj_TCR = clustered_results$singleChainSeuratObject,
  resolution = 0.1
)
```

### Example 3: Format and Validate TCR Metadata

```r
# Clean and validate TCR data before clustering
formatted_metadata <- FormatMetadataForTcrDist3(
  metadata = seuratObj@meta.data,
  chains = c("TRA", "TRB"),
  minimumCloneSize = 2
)
```

## Workflows

### Interactive Analysis

For exploratory analysis with RMarkdown:

```r
# Export example workflow template
GetExampleMarkdown(dest = 'tcrClustR_workflow.Rmd')

# Or view built-in vignettes
browseVignettes("tcrClustR")
```

### Programmatic Pipeline

Complete workflow for automated analysis:

```r
# 1. Format metadata
formatted <- FormatMetadataForTcrDist3(
  metadata = seuratObj@meta.data,
  chains = c("TRA", "TRB"),
  minimumCloneSize = 2
)

# 2. Compute distances
tcr_obj <- RunTcrdist3(
  inputData = seuratObj,
  chains = c("TRA", "TRB"),
  minimumCloneSize = 2
)

# 3. Cluster
results <- RunTcrClustering(
  seuratObj_TCR = tcr_obj,
  outputDir = "./results/"
)
```

## Output File Formats

### Parquet Tables (from RunTcrClustering)

**Single-chain format** (TRA, TRB, TRA_cdr3, TRB_cdr3):
```
chain | v_gene  | j_gene  | CDR3          | Cluster | Cluster_Size_Threshold | Clustering_Method
TRB   | TRBV7-9 | TRBJ2-1 | CASSLGQAYEQYF | 1       | 2                     | DIANA
```

**Paired-chain format** (TRA_TRB, TRACDR3_TRB):
```
chain_1 | v_gene_1 | j_gene_1 | CDR3_1 | chain_2 | v_gene_2 | j_gene_2 | CDR3_2       | Cluster
TRA     | TRAV1-1  | TRAJ1-1  | CAAA   | TRB     | TRBV7-9  | TRBJ2-1  | CASSLGQAYEQYF| 1
```

### Seurat Object Metadata (from ClusterTcrs)

Clustering results stored as columns in metadata:
- `TcrClustR_{assay}_{resolution}` - cluster assignments
- Dimensionality reductions stored in `@reductions` slot (e.g., `TcrClustR_pca.TRA`)

## <a name="issues">Known Issues</a>

- Memory: tcrdist3 switches to sparse matrices for n > 10,000 clones
- Python path: set `RETICULATE_PYTHON` environment variable if tcrdist3 fails
- Seurat v5: always call `JoinLayers()` before accessing assay data in v5 objects
- Gene alleles: current implementation can optionally strip allele notation (e.g., `TRBV7-9*01` → `TRBV7-9`)

## <a name="developers">Development Guidelines</a>

### Running Tests

```r
# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-RunTcrClustering.R")
```

### Building Documentation

```r
# Update documentation
devtools::document()

# Build vignettes
devtools::build_vignettes()

# Check package
devtools::check()
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Make changes with tests
4. Run `devtools::check()` (0 errors/warnings/notes)
5. Submit a pull request

See [copilot-instructions.md](.github/copilot-instructions.md) for detailed development patterns.

## Citation

If you use tcrClustR in your research, please cite:

```
[Citation information to be added]
```

## License

This project is licensed under the MIT License - see [LICENSE.md](LICENSE.md) for details.

## Acknowledgments

- **tcrdist3**: Mayer-Blackwell et al. (https://github.com/kmayerb/tcrdist3)
- **Seurat**: Hao et al., Cell 2021
- Bimber Laboratory, Oregon Health & Science University
