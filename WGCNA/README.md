# WGCNA Docker Container

Dockerfile for running Weighted Gene Co-expression Network Analysis (WGCNA) with supporting R packages in a reproducible, containerized environment. Built for use on HPC clusters via Apptainer/Singularity.

## Included packages

**Bioconductor:** WGCNA, DESeq2, edgeR, limma

**CRAN:** tidyverse, vegan, data.table, glue, stringr, readxl, writexl, RColorBrewer, remotes

**Base image:** [rocker/r-ver:4.3.2](https://hub.docker.com/r/rocker/r-ver) (R 4.3.2)

## Build and deploy

### 1. Build Docker image locally

```bash
docker build -t wgcna .
```

### 2. Save image as tar archive

```bash
docker save wgcna -o wgcna.tar
```

### 3. Transfer to HPC

```bash
scp wgcna.tar username@hpc-cluster:/path/to/containers/
```

### 4. Convert to Apptainer SIF on HPC

```bash
apptainer build wgcna.sif docker-archive://wgcna.tar
```

### 5. Run on HPC

```bash
apptainer exec wgcna.sif Rscript your_script.R
```

Or within a SLURM job script:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32GB
#SBATCH --time=2:00:00

apptainer exec /path/to/wgcna.sif Rscript your_script.R
```

## Usage notes

The Dockerfile copies a default R script (`WGCNA_CHEMM_v1.R`) into the container. To run your own scripts, use `apptainer exec` with your script path rather than the container's default `CMD`, or modify the Dockerfile before building.

## License

MIT
