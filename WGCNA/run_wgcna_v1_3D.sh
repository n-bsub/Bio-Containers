#!/bin/bash 
#SBATCH --nodes=1  
#SBATCH --ntasks-per-node=48 
#SBATCH --time=3:00:00 
#SBATCH --job-name=run_wgcna
#SBATCH --account=
#SBATCH --mail-type=ALL 


set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"

# Point to the same script you call inside Apptainer:
RSCRIPT="$SUBMIT_DIR/WGCNA_CHEMM_v1.3D.R"          # adjust R script name

SCRIPT_NAME="$(basename "$RSCRIPT")"         

# RUN_ID = "<script name>_<jobid>"
RUN_ID="${SCRIPT_NAME}_${SLURM_JOB_ID:-local}"

# Where to write outputs (change to "$SUBMIT_DIR/$RUN_ID" if you don't want a 'runs' parent)
RUN_DIR="$SUBMIT_DIR/runs/$RUN_ID"
mkdir -p "$RUN_DIR"
mkdir -p "$RUN_DIR/modules"

INPUT_DIR="$SUBMIT_DIR/runs/wgcna_input"

# Log everything into the run folder
exec > >(tee -a "$RUN_DIR/run.log") 2>&1

# Launch container with CWD inside the run folder
module load apptainer 2>/dev/null || true
apptainer exec --cleanenv \
  --bind "$SUBMIT_DIR:$SUBMIT_DIR" \
  --bind "$RUN_DIR:$RUN_DIR" \
  --pwd  "$RUN_DIR" \
  --env  INPUT_DIR="$INPUT_DIR" \
  /fs/ess/PAS1182/Neeraja_CHEMM/WGCNA/wgcna_chemm.sif Rscript "$RSCRIPT"

echo "[info] WGCNA step finished."

# 2) Run GOseq per module
SIF_TRINITY="trinityrnaseq_2.12.sif"   

# Checks
[[ -f "$RUN_DIR/wgcna_input_gene_list_background.txt" ]] || { echo "[error] Missing $RUN_DIR/wgcna_input_gene_list_background.txt"; exit 1; }
[[ -f "$INPUT_DIR/go_annotations_CHEMM_G.txt" ]] || { echo "[error] Missing $INPUT_DIR/go_annotations_CHEMM_G.txt"; exit 1; }
[[ -f "$INPUT_DIR/CD_HIT_EST_CHEMM_all.gene_lengths.txt" ]] || { echo "[error] Missing $INPUT_DIR/CD_HIT_EST_CHEMM_all.gene_lengths.txt"; exit 1; }

shopt -s nullglob

for file in "$RUN_DIR"/modules/*_module.txt; do
  base="$(basename "$file")"
  echo "[info] GOseq for $base"

  apptainer exec --no-home --cleanenv \
    --bind "$RUN_DIR:$RUN_DIR" \
    --bind "$INPUT_DIR:$INPUT_DIR" \
    --pwd  "$RUN_DIR" \
    "$SIF_TRINITY" \
    /usr/local/bin/Analysis/DifferentialExpression/run_GOseq.pl \
      --genes_single_factor "$file" \
      --GO_assignments      "$INPUT_DIR/go_annotations_CHEMM_G.txt" \
      --lengths             "$INPUT_DIR/CD_HIT_EST_CHEMM_all.gene_lengths.txt" \
      --background          "$RUN_DIR/wgcna_input_gene_list_background.txt"
done

mkdir -p "$RUN_DIR/Enrichment_results" 

# Collect outputs (regardless of where GOseq wrote them)
mv -f "$RUN_DIR"/modules/*.GOseq.enriched "$RUN_DIR/Enrichment_results/" 2>/dev/null || true
mv -f "$RUN_DIR"/modules/*.GOseq.depleted "$RUN_DIR/Enrichment_results/" 2>/dev/null || true
mv -f "$RUN_DIR"/*.GOseq.enriched         "$RUN_DIR/Enrichment_results/" 2>/dev/null || true
mv -f "$RUN_DIR"/*.GOseq.depleted         "$RUN_DIR/Enrichment_results/" 2>/dev/null || true

echo "[info] GOseq complete → $RUN_DIR/Enrichment_results"

RSCRIPT_COMBINE="$SUBMIT_DIR/combine_enrichment.R"
if [[ -s "$RSCRIPT_COMBINE" ]]; then
  echo "[info] Combining GOseq outputs..."
  apptainer exec --cleanenv \
    --bind "$SUBMIT_DIR:$SUBMIT_DIR" \
    --bind "$RUN_DIR:$RUN_DIR" \
    --pwd  "$RUN_DIR" \
    wgcna_chemm.sif Rscript "$RSCRIPT_COMBINE"
  echo "[info] Combined files are in: $RUN_DIR/GO_enrichment"
else
  echo "[warn] $RSCRIPT_COMBINE not found; skipping combination."
fi

