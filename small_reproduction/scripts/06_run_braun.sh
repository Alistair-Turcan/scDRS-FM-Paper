#!/usr/bin/env bash
# 06_run_braun.sh  — Braun human developing brain (Linnarsson) (human).
# Original: slurms/419_linnarsson.slurm. That slurm's gs=() for GWAS is shadowed
# by the microglia block (an upstream bug); the notebooks actually consume BOTH
# the 8 brain GWAS traits AND the 5 microglia signatures, so we run both here,
# writing to a single output dir (results/real/braun_micro) as nb10 expects.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"
source "$(dirname "${BASH_SOURCE[0]}")/_trait_lists.sh"

H5AD="${SUBSET_DIR}/Braun/human_dev_layers_100k.h5ad"
COV="${SUBSET_DIR}/Braun/human_dev_layers_100k.cov"
OUT="${REAL_OUT}/braun_micro"; mkdir -p "${OUT}"

echo "[06a] Braun brain GWAS (${#BRAIN_TRAITS[@]} traits, scDRS-FM / MAGIC + ctrl)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT}" "${GS_SPLIT}" "${BRAIN_TRAITS[@]}" \
  --h5ad_species human --flag_filter --flag_raw_count --imputation magic --include_ctrl_score \
  2>&1 | tee "${LOG_DIR}/real_braun_gwas.log"

echo "[06b] Braun microglia signatures (${#MICROGLIA_SIGS[@]}, scDRS-FM / MAGIC + ctrl)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT}" "${MICROGLIA}" "${MICROGLIA_SIGS[@]}" \
  --h5ad_species human --flag_filter --flag_raw_count --imputation magic --include_ctrl_score \
  2>&1 | tee "${LOG_DIR}/real_braun_micro.log"

echo "[06] Done -> ${OUT}"
