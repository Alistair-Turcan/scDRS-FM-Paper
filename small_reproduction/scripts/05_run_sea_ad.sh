#!/usr/bin/env bash
# 05_run_sea_ad.sh  — SEA-AD middle temporal gyrus (human).
# Two runs (original: slurms/21_sea_ad.slurm + slurms/27_micro.slurm):
#   (a) 8 brain GWAS traits, scDRS-FM (MAGIC)            -> results/real/sea_ad_brain
#   (b) 5 microglia signatures, scDRS-FM (MAGIC)         -> results/real/sea_ad_micro
# Feeds nb10 (SEA-AD analysis). Note: SEA_AD X is already normalized, so no
# --flag_raw_count (matches the original slurm).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"
source "$(dirname "${BASH_SOURCE[0]}")/_trait_lists.sh"

H5AD="${SUBSET_DIR}/SEA_AD/combined_healthy_filtered.h5ad"
COV="${SUBSET_DIR}/SEA_AD/combined_healthy_filtered.cov"

# (a) brain GWAS
OUT_BR="${REAL_OUT}/sea_ad_brain"; mkdir -p "${OUT_BR}"
echo "[05a] SEA-AD brain GWAS (${#BRAIN_TRAITS[@]} traits, scDRS-FM / MAGIC)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT_BR}" "${GS_SPLIT}" "${BRAIN_TRAITS[@]}" \
  --h5ad_species human --flag_filter --imputation magic \
  2>&1 | tee "${LOG_DIR}/real_sea_ad_brain.log"

# (b) microglia signatures
OUT_MG="${REAL_OUT}/sea_ad_micro"; mkdir -p "${OUT_MG}"
echo "[05b] SEA-AD microglia signatures (${#MICROGLIA_SIGS[@]}, scDRS-FM / MAGIC)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT_MG}" "${MICROGLIA}" "${MICROGLIA_SIGS[@]}" \
  --h5ad_species human --flag_filter --imputation magic \
  2>&1 | tee "${LOG_DIR}/real_sea_ad_micro.log"

echo "[05] Done -> ${OUT_BR} , ${OUT_MG}"
