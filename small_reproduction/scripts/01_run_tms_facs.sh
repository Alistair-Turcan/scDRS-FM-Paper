#!/usr/bin/env bash
# 01_run_tms_facs.sh  — TMS FACS (mouse), all 75 GWAS traits.
# scDRS-FM (MAGIC imputation).  Original: slurms/26_tms_facs.slurm
#
# Produces the primary marginal + conditional scores that drive nb01, nb05, nb06,
# nb08, and the new runtime/memory figure. This is the headline scDRS-FM run.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

H5AD="${SUBSET_DIR}/TMS_FACS/TMS_FACS.h5ad"
COV="${SUBSET_DIR}/TMS_FACS/TMS_FACS.cov"
OUT="${REAL_OUT}/tms_facs"
mkdir -p "${OUT}"

# All 75 traits = every filename in gs_split/ (one per trait)
mapfile -t traits < <(find "${GS_SPLIT}" -maxdepth 1 -type f -printf "%f\n" | sort)
echo "[01] Scoring ${#traits[@]} traits on TMS_FACS (scDRS-FM / MAGIC)"

run_scdrs_fm \
  "${H5AD}" "${COV}" "${OUT}" "${GS_SPLIT}" \
  "${traits[@]}" \
  --h5ad_species mouse \
  --flag_filter \
  --flag_raw_count \
  --imputation magic \
  2>&1 | tee "${LOG_DIR}/real_tms_facs.log"

echo "[01] Done -> ${OUT}"
