# Base image with R and essential system tools
FROM rocker/r-ver:4.3.2

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libtiff5-dev \
    libjpeg-dev \
    libpng-dev \
    libgit2-dev \
    libglpk-dev \
    libxt-dev \
    libx11-dev \
    r-cran-xml \
    r-cran-rgl \
    r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Bioconductor and CRAN packages
RUN Rscript -e "install.packages(c('BiocManager', 'tidyverse', 'data.table', 'glue', 'stringr', 'readxl', 'writexl', 'RColorBrewer'))" \
    && Rscript -e "install.packages(c('vegan', 'remotes'))" \
    && Rscript -e "BiocManager::install(c('DESeq2', 'edgeR', 'limma', 'WGCNA'), ask=FALSE, update=FALSE)"

# Copy your R script into the container
COPY WGCNA_CHEMM_v1.R /home/ruser/WGCNA_CHEMM_v1.R

# Set working directory and default command
WORKDIR /home/ruser
CMD ["Rscript", "WGCNA_CHEMM_v1.R"]
