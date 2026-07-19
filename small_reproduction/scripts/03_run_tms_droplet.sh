#!/usr/bin/env bash
# 03_run_tms_droplet.sh  — TMS Droplet (mouse), all 75 GWAS traits.
# scDRS baseline with imputation=none + ctrl scores.  Original: slurms/222_tmsd.slurm
# Feeds nb07 (TMS_Droplet heatmaps) and nb08 (genetic correlation).
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

H5AD="${SUBSET_DIR}/TMS_Droplet/TMS_Droplet.h5ad"
COV="${SUBSET_DIR}/TMS_Droplet/tms_droplet.cov"
OUT="${REAL_OUT}/tms_droplet"
mkdir -p "${OUT}"

mapfile -t traits < <(find "${GS_SPLIT}" -maxdepth 1 -type f -printf "%f\n" | sort)
echo "[03] Scoring ${#traits[@]} traits on TMS_Droplet (scDRS / none + ctrl)"

run_scdrs_fm \
  "${H5AD}" "${COV}" "${OUT}" "${GS_SPLIT}" \
  "${traits[@]}" \
  --h5ad_species mouse \
  --flag_filter \
  --flag_raw_count \
  --imputation none \
  --include_ctrl_score \
  2>&1 | tee "${LOG_DIR}/real_tms_droplet.log"

echo "[03] Done -> ${OUT}"
