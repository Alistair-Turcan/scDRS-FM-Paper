#!/usr/bin/env bash
# 02_run_ts_facs.sh  — Tabula Sapiens FACS (human), all 75 GWAS traits.
# scDRS baseline with imputation=none + ctrl scores.  Original: slurms/222_ts_facs.slurm
# Feeds nb07 (TS_FACS heatmaps) and nb08 (genetic correlation).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

H5AD="${SUBSET_DIR}/TS_FACS/ts_facs.h5ad"
COV="${SUBSET_DIR}/TS_FACS/ts.cov"
OUT="${REAL_OUT}/ts_facs"
mkdir -p "${OUT}"

mapfile -t traits < <(find "${GS_SPLIT}" -maxdepth 1 -type f -printf "%f\n" | sort)
echo "[02] Scoring ${#traits[@]} traits on TS_FACS (scDRS / none + ctrl)"

run_scdrs_fm \
  "${H5AD}" "${COV}" "${OUT}" "${GS_SPLIT}" \
  "${traits[@]}" \
  --h5ad_species human \
  --flag_filter \
  --flag_raw_count \
  --imputation none \
  --include_ctrl_score \
  2>&1 | tee "${LOG_DIR}/real_ts_facs.log"

echo "[02] Done -> ${OUT}"
