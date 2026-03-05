# tcrClustR

An R package for clustering and analyzing T-cell receptor (TCR) sequences to identify 'TCR families' via sequence similarity. This package uses [tcrdist3](https://github.com/kmayerb/tcrdist3) for TCR distance calculations and provides flexible clustering algorithms to identify groups of functionally related TCRs. Similar in concept to GLIPH and CoNGA, but with more direct control over clustering parameters.

Full documentation: https://bimberlabinternal.github.io/tcrClustR/

## Table of Contents
* [Quick Start](#quick-start)
* [Overview](#overview)
* [Installation](#installation)
* [Dirichlet Process Analysis](#dirichlet-process-analysis)
* [Known Issues](#issues)
* [Development Guidelines](#developers)

## Quick Start

The fastest way to cluster TCR data:

```r
library(tcrClustR)

# Step 1: compute TCR distance matrices (stored in seuratObj@misc$TCR_Distances)
seuratObj <- CalculateTcrDistances(
  inputData = seuratObj,
  chains = c("TRA", "TRB"),
  minimumCloneSize = 2,
  calculateChainPairs = TRUE
)

# Step 2: cluster TCRs via DIANA and store results in metadata
seuratObj <- RunTcrClustering(
  seuratObj_TCR = seuratObj,
  dianaHeight = 20,
  clusterSizeThreshold = 1
)

# Cluster assignments live in metadata columns like TRB_fl_ClusterIdx
DimPlot(seuratObj, reduction = "umap", group.by = "TRB_fl_ClusterIdx", label = TRUE)

# Retrieve a raw distance matrix
distance_mat <- GetDistanceMatrix(seuratObj, chains = "TRA")
```

## Overview

### What is TCR Clustering?

T-cell receptor (TCR) clustering groups TCR sequences based on similarity metrics, enabling identification of functionally related TCR families from single-cell sequencing data. TCRs with similar sequences (CDR3 regions and V/J gene segments) may recognize the same or related antigens.

### Workflow

1. **Format & Validate**: Clean TCR metadata, filter low-quality clones
2. **Compute Distances**: Calculate pairwise TCR distances via tcrdist3 (BLOSUM62 matrix)
3. **Cluster**: Apply DIANA hierarchical clustering
4. **Analyse & Tune**: Use Dirichlet process mixture models to discover natural modes in the distance distribution and guide `dianaHeight` selection
5. **Visualize**: Heatmaps, histograms, cluster-mean error bars, and mixing-proportion charts

### Key Features

- **DIANA clustering**: Divisive hierarchical clustering with a tuneable height cutoff
- **Dirichlet process analysis**: Non-parametric Gaussian mixture modelling of pairwise distances to identify natural cluster modes and inform `dianaHeight`
- **Single & paired chain analysis**: TRA, TRB, TRG, TRD, or combined TRA+TRB distances
- **Quantile-stratified downsampling**: Preserves rare modes in the tails of distance distributions
- **Automatic data filtering**: Remove NA and concatenated values (optional)

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
pip install git+https://github.com/bimberlabinternal/tcrdist3.git@0.3

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

```

For exploratory analysis with RMarkdown:

```r
# Export example workflow template
GetExampleMarkdown(dest = 'tcrClustR_workflow.Rmd')

# Or view built-in vignettes
browseVignettes("tcrClustR")
```

## Dirichlet Process Analysis

`DirichletClusterAnalysis()` fits a non-parametric Gaussian Dirichlet process mixture to within-group pairwise TCR distances. The discovered cluster means (mu) and spreads (sigma) reveal the natural modes of the distance distribution, which you can use to select an informed `dianaHeight` cutoff for `RunTcrClustering()`.

### Basic Usage

```r
# Fit DP mixture models per group
dp <- DirichletClusterAnalysis(
  seuratObj   = seuratObj,
  assayName   = "TRA_fl",
  splitField  = "Population",
  maxSamples  = 1000,
  nIterations = 500
)

# Two diagnostic plots (combine with patchwork)
library(patchwork)
(PlotClusterMeans(dp) + Seurat::NoLegend()) +
  PlotMixingProportions(dp) +
  plot_layout(guides = "collect")

# Inspect the cluster parameter table
glimpse(dp$cluster_summary)
#Rows: 11
#Columns: 6
#$ Cluster          <int> 1, 2, 3, 4, 5, 6, 1, 2,…
#$ Mu               <dbl> 135.98787492, 31.812791…
#$ Sigma            <dbl> 4.9843801, 3.4540618, 1…
#$ MixingProportion <dbl> 0.79400000, 0.17100000,…
#$ PointsPerCluster <int> 794, 171, 31, 1, 2, 1, …
#$ Group            <chr> "MR1-5-OP-RU-Tet", "MR1…
```

### Quantile-Stratified Downsampling

TCR distance distributions often have small, rare modes in the tails that would be missed by uniform random sampling. `DirichletClusterAnalysis()` uses quantile-stratified (n-tile) downsampling by default: distances are divided into `nBins` equal-frequency bins and up to `samplesPerBin` values are drawn from each, preserving the full distributional shape. Because the Dirichlet Proces fitting can be computationally taxing, `maxSamples` downsamples the quantile-stratified population if `nBins` * `samplesPerBin` > `maxSamples`.

Increase `nBins` and `samplesPerBin` (and `maxSamples`) to improve resolution of rare modes:

```r
dp <- DirichletClusterAnalysis(
  seuratObj     = seuratObj,
  assayName     = "TRA_fl",
  splitField    = "Population",
  maxSamples    = 1000,
  nIterations   = 500,
  nBins         = 20,
  samplesPerBin = 150
)
```

### Output

`DirichletClusterAnalysis()` returns a `tcrDirichletResult` list containing:

| Field | Description |
|---|---|
| `cluster_summary` | Tidy data.frame: one row per group × cluster with `Mu`, `Sigma`, `MixingProportion`, `PointsPerCluster`, `Group` |
| `models` | Named list of raw `dirichletprocess` model objects for downstream inspection |
| `assayName` | The distance assay that was analysed |
| `splitField` | The metadata column used for grouping |

### Companion Plots

- **`PlotClusterMeans(dp)`** — Error-bar plot of mu ± sigma per cluster, dodged by group. The y-axis corresponds directly to TCR distance and can be compared against `dianaHeight`.
- **`PlotMixingProportions(dp)`** — Dodged bar chart of cluster mixing weights. Clusters with high weight and low mu identify well-supported clonotype families.

### Helper: ExtractGroupDistanceVectors

If you need the raw per-group distance vectors without fitting a DP model (e.g., for your own downstream analysis):

```r
vecs <- ExtractGroupDistanceVectors(
  seuratObj = seuratObj,
  assayName = "TRB_cdr3",
  splitField = "Tissue"
)
# Returns a named list of numeric vectors, one per group
hist(vecs[["Spleen"]], breaks = 50)
```

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
testthat::test_file("tests/testthat/test-runTcrClustering.R")
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


## Citation

If you use tcrClustR in your research, please cite:

```
[Citation information to be added]
```

## License

This project is licensed under the GPL (>= 3) License — see [LICENSE.md](LICENSE.md) for details.

## Acknowledgments

- **tcrdist3**: Mayer-Blackwell et al. (https://github.com/kmayerb/tcrdist3)
- **Seurat**: Hao et al., Cell 2021
- Bimber Laboratory, Oregon Health & Science University
