FROM rocker/r-base:4.4.2

ARG DEBIAN_FRONTEND=noninteractive
ARG GH_PAT='NOT_SET'

ADD . /tcrClustR

RUN apt-get update && apt-get install -y \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    uuid-dev \
    libxml2-dev \
    libgpgme11-dev \
    squashfs-tools \
    libseccomp-dev \
    libsqlite3-dev \
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
    libjpeg-dev && \
    mkdir /TCR_Python && \
    cd /TCR_Python && \
    wget https://github.com/python/cpython/archive/refs/tags/v3.8.10.tar.gz && \
    tar -zxvf v3.8.10.tar.gz && \
    cd cpython-3.8.10 && \
    ./configure --prefix=/TCR_Python && \
    cd /TCR_Python/cpython-3.8.10 && \
    make && \
    make install && \
    /TCR_Python/bin/pip3 --no-cache-dir install numpy scipy scikit-learn scikit-misc matplotlib tqdm sympy setuptools pandas pyyaml scanpy && \
    /TCR_Python/bin/pip3 --no-cache-dir install git+https://github.com/kmayerb/tcrdist3.git@0.2.2 && \
    # Install conga:
    mkdir /conga  && \
    cd /conga && \
    git clone https://github.com/phbradley/conga.git && \
    cd conga/tcrdist_cpp && \
    make && \
    cd ../ && \
    /TCR_Python/bin/pip3 install -e . && \
    cd / && \
    chmod -R 777 /TCR_Python


#install R and dependencies
RUN apt-get update && apt-get install -y r-base r-base-dev && \
    if [ "${GH_PAT}" != 'NOT_SET' ]; then \
        echo 'Setting GH_PAT'; \
        export GITHUB_PAT="${GH_PAT}"; \
    fi && \
    Rscript -e "install.packages(c('remotes', 'devtools', 'BiocManager', 'pryr', 'rmdformats', 'knitr', 'logger', 'Matrix', 'kernlab', 'tidyverse', 'Seurat', 'leidenbase', 'igraph', 'FNN'), dependencies=TRUE, ask = FALSE, upgrade = 'always')" && \
    echo "local({options(repos = BiocManager::repositories())})" >> ~/.Rprofile


#build tcrClustR
RUN cd /tcrClustR && \
    R CMD build . && \
    Rscript -e "BiocManager::install(ask = F, upgrade = 'always');" && \
    Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" && \
    R CMD INSTALL --build *.tar.gz && \
    rm -Rf /tmp/downloaded_packages/ /tmp/*.rds

ENV NUMBA_CACHE_DIR=/work/numba_cache
ENV MPLCONFIGDIR=/work/mpl_cache

ENTRYPOINT ["/bin/bash"]
