# Multi-stage Dockerfile for tcrClustR
# Stage 1: Base - System dependencies and base tools
# Stage 2: Deps - R packages and Python packages (cached layer)
# Stage 3: Runtime - Final application build

# ============================================================================
# Stage 1: Base - System dependencies
# ============================================================================
ARG BASE_IMAGE=rocker/r-base:4.4.2
FROM ${BASE_IMAGE} AS base

ARG DEBIAN_FRONTEND=noninteractive
ARG SKIP_BASE_DEPS=false

# Install system dependencies (skipped if using pre-built base image)
# To use a pre-built base: --build-arg BASE_IMAGE=ghcr.io/.../base-deps:tag --build-arg SKIP_BASE_DEPS=true
RUN if [ "$SKIP_BASE_DEPS" = "false" ]; then \
    apt-get update && apt-get install -y \
        build-essential \
        libcurl4-openssl-dev \
        libssl-dev \
        uuid-dev \
        libxml2-dev \
        libgpgme11-dev \
        squashfs-tools \
        libseccomp-dev \
        r-cran-devtools \
        libsqlite3-dev \
        libgit2-dev \
        pkg-config \
        git-all \
        wget \
        libbz2-dev \
        zlib1g-dev \
        python3-dev \
        libffi-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libmbedtls-dev \
        cargo \
        libmagick++-dev \
        libudunits2-dev \
        libgsl-dev \
        libtbb-dev \
        cmake \
        libcairo2-dev \
        libgpg-error-dev \
        libgmp-dev \
        ca-certificates \
        && rm -rf /var/lib/apt/lists/*; \
    else \
        echo "Skipping base dependency installation (using pre-built base image)"; \
    fi

# ============================================================================
# Stage 2: Deps - Install R and Python dependencies
# ============================================================================
FROM base AS deps

ARG DEBIAN_FRONTEND=noninteractive
ARG GH_PAT='NOT_SET'

# Install Python packages in a virtual environment
RUN apt-get update && \
    apt-get install -y python3-pip python3-venv && \
    python3 -m venv /opt/venv && \
    . /opt/venv/bin/activate && \
    pip install --upgrade pip && \
    pip install --no-cache-dir \
        numpy scipy scikit-learn scikit-misc matplotlib tqdm sympy \
        setuptools pandas pyyaml scanpy rpy2 && \
    pip install --no-cache-dir git+https://github.com/kmayerb/tcrdist3.git@0.2.2 && \
    # Install conga - remove if exists (for base image compatibility)
    rm -rf /conga && \
    mkdir -p /conga && \
    cd /conga && \
    git clone https://github.com/phbradley/conga.git && \
    cd conga/tcrdist_cpp && \
    make && \
    cd ../ && \
    pip install -e . && \
    cd / && \
    rm -rf /var/lib/apt/lists/*

# Add virtual environment to PATH so Python scripts can find packages
ENV PATH="/opt/venv/bin:$PATH"

# Install R dependencies
RUN apt-get update && apt-get install -y r-base r-base-dev && \
    if [ "${GH_PAT}" != 'NOT_SET' ]; then \
        echo 'Setting GH_PAT'; \
        export GITHUB_PAT="${GH_PAT}"; \
    fi && \
    Rscript -e "install.packages(c('remotes', 'devtools', 'BiocManager', 'pryr', 'rmdformats', 'knitr', 'logger', 'Matrix', 'kernlab', 'tidyverse', 'Seurat', 'leidenbase', 'igraph', 'FNN'), lib='/usr/local/lib/R/site-library', dependencies=TRUE, ask = FALSE, upgrade = 'always')" && \
    echo "local({options(repos = BiocManager::repositories())})" >> ~/.Rprofile && \
    Rscript -e "BiocManager::install('ComplexHeatmap', ask = FALSE, update = TRUE)" && \
    Rscript -e "install.packages(c('clusterCrit', 'dbscan', 'cluster', 'arrow', 'RColorBrewer', 'patchwork'), lib='/usr/local/lib/R/site-library', dependencies=TRUE, ask = FALSE)" && \
    Rscript -e "remotes::install_github('kevinsblake/NatParksPalettes', lib='/usr/local/lib/R/site-library')" && \
    rm -rf /var/lib/apt/lists/* /tmp/downloaded_packages/ /tmp/*.rds

# ============================================================================
# Stage 3: Runtime - Build and install tcrClustR package
# ============================================================================
# Use a pre-built deps image when provided to dramatically speed up runtime builds
# Fallback to the local `deps` stage if DEPS_IMAGE is not set (e.g., local builds)
# Note: Dockerfile variable expansion doesn't support ${VAR:-default} reliably; set default at ARG.
ARG DEPS_IMAGE=deps
FROM ${DEPS_IMAGE} AS runtime

ARG DEBIAN_FRONTEND=noninteractive
ARG GH_PAT='NOT_SET'

# Copy application code
ADD . /tcrClustR

# Build and install tcrClustR
# upgrade = 'never' is specified because the dependencies install/upgrading is already handled in the deps stage
RUN cd /tcrClustR && \
    R CMD build . && \
    Rscript -e "BiocManager::install(ask = F, upgrade = 'never');" && \
    Rscript -e "install.packages('devtools', dependencies=TRUE, lib='/usr/local/lib/R/site-library'); \
    devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'never');" && \
    R CMD INSTALL --build *.tar.gz && \
    rm -Rf /tmp/downloaded_packages/ /tmp/*.rds

ENV NUMBA_CACHE_DIR=/work/numba_cache
ENV MPLCONFIGDIR=/work/mpl_cache

ENTRYPOINT ["/bin/bash"]
