![R Build and Checks](https://github.com/bimberlabinternal/tcrclustr/workflows/R%20Build%20and%20Checks/badge.svg)

# tcrClustR
An R package designed to cluster and analyze TCR sequences, specifically to the point of defining 'TCR families' via clustering. This package primarily uses tcrdist3 for TCR distance calculations and is similar in concept to GLIPH (Grouping of Lymphocyte Interactions by Paratope Hotspots) and CoNGA (Clonotype Neighbor-graph Analysis), but exposes the clustering parameters more directly. [Please see our documentation for more detail](https://bimberlabinternal.github.io/tcrClustR/).

## Table of Contents
* [Overview](#overview)
* [Example Usage](#example)
* [Installation](#installation)
* [Known Issues](#issues)
* [Development Guidelines](#developers)

### <a name="overview">Overview</a>

T-cell receptor (TCR) clustering is a computational method that groups TCR sequences based on measures of their sequence similarity, allowing for the identification of functionally related TCR families within single-cell TCR-seq data. The general idea is to compute distance matrices between TCR sequences using tcrdist3, then apply network-based clustering methods to identify groups of similar TCRs that may recognize the same or related antigens.

TCR clustering analysis begins with TCR sequence data, typically extracted from single-cell RNA-seq experiments that include TCR information. The sequences are processed to extract key components (CDR3 regions, V/J gene segments), and distance metrics are calculated between all pairs of TCRs using tcrdist3, which accounts for sequence similarity, based on the BLOSUM substitution matrix. 

Once distance matrices are computed, clustering algorithms must be applied to identify TCR families or groups. This package provides several functions:
- Data formatting and quality control for TCR sequence data, including metadata cleaning and validation for tcrdist3 compatibility
- A streamlined interface to run tcrdist3 for TCR distance matrix computation, supporting both single-chain (TRA or TRB) and multi-chain analyses
- Network-based clustering using methods like Leiden clustering to identify TCR families from distance matrices
- Visualization tools including heatmaps and histograms to explore TCR similarity patterns and clustering results

Each step of the workflow can either be run interactively in R (through the terminal or RStudio), or it can be executed as a pipeline that processes TCR data and generates clustering results with accompanying visualizations.

### <a name="example">Example Usage</a>

Below are the primary functions of tcrClustR needed to process and cluster TCR sequence data:
```r
# Example 1: Format TCR metadata for tcrdist3 analysis
formatted_metadata <- FormatMetadataForTcrDist3(
  metadata = seuratObj@meta.data,
  chains = c("TRA", "TRB"),
  cleanMetadata = TRUE,
  summarizeClones = TRUE
)

# Example 2: Run tcrdist3 to compute TCR distance matrices
# This casts the data into a "per-clonotype" grouping, rather than "per-cell". 
seuratObj_TCR <- RunTcrdist3(
  seuratObj = seuratObj,
  chains = c("TRA", "TRB"),
  cleanMetadata = TRUE,
  minimumClonesPerSubject = 2,
  rdsOutputPath = "./tcrdist3DistanceMatrices/"
)

# Example 3: Cluster TCRs using network-based methods
clustered_results <- ClusterTcrs(
  seuratObj_TCR = seuratObj_TCR,
  resolutionParameter = 0.1,
  usePCA = TRUE,
  pcaComponents = 50
)

# Example 4: Generate TCR distance heatmaps
TCRDistanceHeatmaps(
  seuratObj_TCR = seuratObj_TCR,
  resolution = 0.1,
  annotate_clusters = TRUE
)

# Example 5: Plot TCR distance histograms by cluster
TCRDistanceHistograms(
  seuratObj_TCR = seuratObj_TCR,
  resolution = 0.1
)
```

Or export/save a template RMarkdown file outlining the default workflow, which can be run interactively or as part of a pipeline:
 
```r
# Get example workflow (when implemented)
GetExampleMarkdown(dest = 'tcrClustR_template.rmd')
```

Finally, the workflow can be executed using wrapper functions that process TCR data and generate clustering results with QC reports:
 
```r
# Complete TCR clustering workflow using tcrdist3
seuratObj_TCR <- RunTcrdist3(
  seuratObj = seuratObj, 
  chains = c("TRA", "TRB"),
  cleanMetadata = TRUE,
  minimumClonesPerSubject = 2,
  rdsOutputPath = "./tcr_analysis/"
)

clustered_results <- ClusterTcrs(
  seuratObj_TCR = seuratObj_TCR,
  resolutionParameter = 0.1
)
```
### <a name="installation">Installation</a>

```{r}
# Make sure to update your Rprofile to include Bioconductor repos, such as adding this line to ~/.Rprofile:
local({options(repos = BiocManager::repositories())})

#Latest version:
devtools::install_github(repo = 'bimberlabinternal/tcrclustr', ref = 'master', dependencies = TRUE, upgrade = 'always')
```

**Python Dependencies:** This package requires Python with tcrdist3 and rpy2 packages for TCR distance calculations. You can install these using:

```bash
pip install tcrdist3 rpy2
```

Pre-packaged Docker images with all needed dependencies installed can be found on our [GitHub Packages page](https://github.com/orgs/BimberLabInternal/packages/container/package/tcrclustr). We recommend using a specific release, which you can do using tags:    

```
docker pull ghcr.io/bimberlabinternal/tcrclustr:latest
```

### <a name="issues">Known Issues</a>

If you receive an error related to Python environment or tcrdist3 package:
```
"Python environment may not have required packages (tcrdist3, rpy2)"
```
Please ensure that tcrdist3 and rpy2 are installed in your Python environment:
```bash
pip install tcrdist3 rpy2
```

**Memory considerations:** TCR distance matrix computation can be memory-intensive for large datasets. Consider:
- Increasing your R memory limit using `memory.limit()` on Windows
- Using the `minimumClonesPerSubject` parameter to filter rare clones
- Processing subsets of your data if memory limitations persist

### <a name="developers">Development Guidelines</a>

* New development should occur on a branch, and go through a Pull Request before merging into the master branch.  [See here for information on the pull request workflow](https://guides.github.com/introduction/flow/).  Ideally PRs would be reviewed by another person.  For the PR, please review the set of changed files carefully to make sure you are only merging the changes you intend.   

* New functions should have [Roxygen2 documentation](https://kbroman.org/pkg_primer/pages/docs.html).

* As part of each PR, you should run 'devtools::document()' to update documentation and include these changes with your commits.

* It is a good idea to run 'R CMD check' locally to make sure your changes will pass.  [See here for more information](http://r-pkgs.had.co.nz/check.html)

* Code should only be merged after the build and tests pass.  The master branch should always be stable.

* New features should ideally have at least a basic test (see [R testthat](http://r-pkgs.had.co.nz/tests.html)).  There is existing test data in ./tests/testdata.  This can be expanded, but please be conscious about file size and try to reuse data across tests if appropriate.

* When adding new TCR distance algorithms or clustering methods, ensure they follow the established patterns for data input/output and integrate properly with the Seurat object framework used throughout the package. The package is primarily designed around tcrdist3 for distance calculations.

