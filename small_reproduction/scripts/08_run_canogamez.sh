#!/usr/bin/env bash
# 08_run_canogamez.sh  — Cano-Gamez et al. naive/memory CD4 T cells (human).
# Two runs (original: slurms/419_canogamez.slurm):
#   (a) 15 immune GWAS traits, scDRS-FM (MAGIC) + ctrl    -> results/real/canogamez_immune
#   (b) 52 T-cell phenotype signatures, scDRS-FM (MAGIC)  -> results/real/canogamez_tcell
# Feeds nb09 (independent replication dataset) and nb11.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"
source "$(dirname "${BASH_SOURCE[0]}")/_trait_lists.sh"

H5AD="${SUBSET_DIR}/Cano_Gamez/obj_raw.h5ad"
COV="${SUBSET_DIR}/Cano_Gamez/canogamez.cov"

OUT_IMM="${REAL_OUT}/canogamez_immune"; mkdir -p "${OUT_IMM}"
echo "[08a] Cano-Gamez immune GWAS (${#IMMUNE_TRAITS[@]} traits, scDRS-FM / MAGIC + ctrl)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT_IMM}" "${GS_SPLIT}" "${IMMUNE_TRAITS[@]}" \
  --h5ad_species human --flag_filter --flag_raw_count --imputation magic --include_ctrl_score \
  2>&1 | tee "${LOG_DIR}/real_canogamez_immune.log"

OUT_TC="${REAL_OUT}/canogamez_tcell"; mkdir -p "${OUT_TC}"
echo "[08b] Cano-Gamez T-cell phenotypes (${#TCELL_SIGS[@]} signatures, scDRS-FM / MAGIC)"
run_scdrs_fm "${H5AD}" "${COV}" "${OUT_TC}" "${TCELL_PHENO}" "${TCELL_SIGS[@]}" \
  --h5ad_species human --flag_filter --flag_raw_count --imputation magic \
  2>&1 | tee "${LOG_DIR}/real_canogamez_tcell.log"

echo "[08] Done -> ${OUT_IMM} , ${OUT_TC}"
