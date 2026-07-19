#!/usr/bin/env bash
# 04_run_soskic.sh  — Soskic T-cell activation (human).
# Two runs (original: slurms/26_soskic.slurm + slurms/27_tcell.slurm):
#   (a) 15 immune GWAS traits, scDRS-FM (MAGIC)              -> results/real/soskic_immune
#   (b) 52 T-cell phenotype signatures, scDRS-FM (MAGIC)     -> results/real/soskic_tcell
# Feeds nb09 (Soskic t-cell analysis) and nb11 (phenotype decomposition).
#
# The .scdrs_ct cell-type tables consumed by nb09 are produced downstream by
# scripts/13_make_ct_tables.sh (or the notebook's own rescoring helpers).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"
source "$(dirname "${BASH_SOURCE[0]}")/_trait_lists.sh"

H5AD="${SUBSET_DIR}/Soskic/soskic_100k.h5ad"
COV="${SUBSET_DIR}/Soskic/soskic.cov"

# (a) immune GWAS
OUT_IMM="${REAL_OUT}/soskic_immune"; mkdir -p "${OUT_IMM}"
echo "[04a] Soskic immune GWAS (${#IMMUNE_TRAITS[@]} traits, scDRS-FM / MAGIC)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT_IMM}" "${GS_SPLIT}" "${IMMUNE_TRAITS[@]}" \
  --h5ad_species human --flag_filter --flag_raw_count --imputation magic \
  2>&1 | tee "${LOG_DIR}/real_soskic_immune.log"

# (b) T-cell phenotype signatures
OUT_TC="${REAL_OUT}/soskic_tcell"; mkdir -p "${OUT_TC}"
echo "[04b] Soskic T-cell phenotypes (${#TCELL_SIGS[@]} signatures, scDRS-FM / MAGIC)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT_TC}" "${TCELL_PHENO}" "${TCELL_SIGS[@]}" \
  --h5ad_species human --flag_filter --flag_raw_count --imputation magic \
  2>&1 | tee "${LOG_DIR}/real_soskic_tcell.log"

echo "[04] Done -> ${OUT_IMM} , ${OUT_TC}"
