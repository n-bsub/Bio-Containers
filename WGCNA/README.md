# WGCNA Docker Container

Dockerfile and analysis script for running Weighted Gene Co-expression Network Analysis (WGCNA) with supporting R packages in a reproducible, containerized environment. Built for use on HPC clusters via Apptainer/Singularity.

## Overview

The included R script (`WGCNA_CHEMM_v1.R`) performs a complete WGCNA pipeline on metatranscriptomic data:

- TMM normalization and voom transformation (edgeR/limma) of gene count matrices
- Optional filtering to fungal-annotated genes using SwissProt BLAST annotations
- Sample clustering and soft threshold power selection for scale-free topology
- Signed hybrid network construction with blockwise module detection
- Module eigengene computation and module–trait correlation heatmaps
- Gene-module mapping with functional annotation export

## Included packages

**Bioconductor:** WGCNA, DESeq2, edgeR, limma

**CRAN:** dplyr, tibble, tidyr, stringr, ggplot2, data.table, glue, readxl, writexl, vegan, RColorBrewer, matrixStats, remotes

**Base image:** [rocker/r-ver:4.3.2](https://hub.docker.com/r/rocker/r-ver) (R 4.3.2)

## Repository contents

```
wgcna/
├── Dockerfile              Container recipe with all R dependencies
├── WGCNA_CHEMM_v1.R        Complete WGCNA analysis pipeline
└── README.md
```

## Input files

The script reads input from a directory specified via the `INPUT_DIR` environment variable. The following files are expected:

| File | Description |
|---|---|
| `RSEM.gene.counts.matrix` | Trinity/RSEM raw gene-level count matrix (tab-delimited) |
| `RSEM.gene.TMM.EXPR.matrix` | TMM-normalized expression matrix (tab-delimited) |
| `Blast_GO_wide_CHEMM.RDS` | SwissProt BLAST annotations with GO terms (R object) |
| `Metadata_wgcna.txt` | Sample metadata with trait variables for module–trait correlation (tab-delimited, must contain a `Samples` or `Sample` column) |

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

Pass the input directory as an environment variable:

```bash
apptainer exec --env INPUT_DIR=/path/to/input/data wgcna.sif Rscript WGCNA_CHEMM_v1.R
```

Or within a SLURM job script:

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=32GB
#SBATCH --time=2:00:00

export INPUT_DIR=/path/to/input/data

apptainer exec /path/to/wgcna.sif Rscript /home/ruser/WGCNA_CHEMM_v1.R
```

## Outputs

The script generates the following output files in the working directory:

| File | Description |
|---|---|
| `gene_module_mapping.txt` | Gene-to-module assignments with functional annotations |
| `modules/<color>_module.txt` | Individual module gene lists |
| `Module_Trait_Correlation_Heatmap_Dust.pdf` | Module–trait relationship heatmap |
| `Module_Trait_Correlations.csv` | Correlation values (modules × traits) |
| `Module_Trait_Pvalues.csv` | P-values for module–trait correlations |
| `MEs.RDS` | Module eigengenes (R object) |
| `norm_data.RDS` | Normalized expression matrix (R object) |
| `module_colors.RDS` | Module color assignments (R object) |
| `voom_halla_sprot.txt` | Voom-normalized data with SwissProt IDs for downstream integration |
| `sft_plot.pdf` | Scale-free topology fit and mean connectivity plots |
| `samples_cluster.pdf` | Sample dendrogram |

## Usage notes

The Dockerfile copies `WGCNA_CHEMM_v1.R` into the container at `/home/ruser/`. To run your own WGCNA script instead, use `apptainer exec` with your script path, the container provides all necessary R packages regardless of which script is run.

## License

MIT
