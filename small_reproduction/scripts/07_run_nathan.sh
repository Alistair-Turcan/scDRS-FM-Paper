#!/usr/bin/env bash
# 07_run_nathan.sh  — Nathan et al. T cells (human).
# Two runs (original: slurms/223_nathan.slurm):
#   (a) 15 immune GWAS traits, scDRS-FM (MAGIC) + ctrl   -> results/real/nathan_immune
#   (b) 52 T-cell phenotype signatures, scDRS-FM (MAGIC) -> results/real/nathan_tcell
# Feeds nb09 (as an independent replication dataset) and nb11.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"
source "$(dirname "${BASH_SOURCE[0]}")/_trait_lists.sh"

H5AD="${SUBSET_DIR}/Nathan/raw.h5ad"
COV="${SUBSET_DIR}/Nathan/raw.cov"

OUT_IMM="${REAL_OUT}/nathan_immune"; mkdir -p "${OUT_IMM}"
echo "[07a] Nathan immune GWAS (${#IMMUNE_TRAITS[@]} traits, scDRS-FM / MAGIC + ctrl)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT_IMM}" "${GS_SPLIT}" "${IMMUNE_TRAITS[@]}" \
  --h5ad_species human --flag_filter --flag_raw_count --imputation magic --include_ctrl_score \
  2>&1 | tee "${LOG_DIR}/real_nathan_immune.log"

OUT_TC="${REAL_OUT}/nathan_tcell"; mkdir -p "${OUT_TC}"
echo "[07b] Nathan T-cell phenotypes (${#TCELL_SIGS[@]} signatures, scDRS-FM / MAGIC)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT_TC}" "${TCELL_PHENO}" "${TCELL_SIGS[@]}" \
  --h5ad_species human --flag_filter --flag_raw_count --imputation magic \
  2>&1 | tee "${LOG_DIR}/real_nathan_tcell.log"

echo "[07] Done -> ${OUT_IMM} , ${OUT_TC}"
