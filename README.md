# Containers

Dockerfiles and analysis scripts for reproducible bioinformatics workflows on HPC clusters via Apptainer/Singularity.

## Custom containers

### [WGCNA-GOseq Containerized Pipeline](wgcna/)

Complete pipeline for Weighted Gene Co-expression Network Analysis with per-module GO enrichment. Chains a custom WGCNA container (voom normalization, signed hybrid network construction, module detection, module–trait correlation) with a Trinity container (GOseq enrichment), orchestrated by a SLURM job script. See the [wgcna/ README](wgcna/README.md) for full documentation.

## Pre-built containers

### Trinity (v2.12)

Pre-built image pulled from Docker Hub and converted to Apptainer for HPC use. Used in the WGCNA–GOseq pipeline for per-module GOseq enrichment via Trinity's `run_GOseq.pl`.

**Source:**
- [Dockerfile](https://github.com/trinityrnaseq/trinityrnaseq/tree/master/Docker)
- [Build script](https://github.com/trinityrnaseq/trinityrnaseq/blob/master/Docker/build_docker.sh)
- [Docker Hub](https://hub.docker.com/r/trinityrnaseq/trinityrnaseq)

**Pull and convert:**

```bash
docker pull trinityrnaseq/trinityrnaseq:2.12.0
docker save trinityrnaseq/trinityrnaseq:2.12.0 -o trinity.tar
scp trinity.tar username@hpc-cluster:/path/to/containers/
apptainer build trinityrnaseq_2.12.sif docker-archive://trinity.tar
```

**Citation:**

> Grabherr MG, Haas BJ, Yassour M, et al. Full-length transcriptome assembly from RNA-Seq data without a reference genome. *Nature Biotechnology* 29, 644–652 (2011). https://doi.org/10.1038/nbt.1883

## General workflow: Docker to Apptainer on HPC

```bash
# Build or pull image
docker build -t <image_name> .
# or
docker pull <repository>:<tag>

# Save as tar
docker save <image_name> -o <image_name>.tar

# Transfer to HPC
scp <image_name>.tar username@hpc-cluster:/path/to/containers/

# Convert to Apptainer SIF
apptainer build <image_name>.sif docker-archive://<image_name>.tar

# Run
apptainer exec <image_name>.sif <command>
```

## License

MIT
