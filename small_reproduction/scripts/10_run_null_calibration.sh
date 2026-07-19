#!/usr/bin/env bash
# 10_run_null_calibration.sh  — FDR/type-I calibration on null gene sets (TMS_FACS).
# Original: slurms/419_null.slurm (SLURM array 0-11 = 4 gene-set families x 3 sizes).
#
#   families = { all, highmean, highvar, overdispersed }
#   sizes    = { 100, 500, 1000 } genes
#   reps     = 0..99  (100 random null gene sets per family x size)
#
# Null gene sets live in data/gene_sets/gs_files_null/ named
#   {family}_ngene{size}_rep{rep}. scDRS-FM (MAGIC). Used to confirm the FDP
# stays at/below nominal (calibration panel context for nb03).
#
# WARNING: this is 12 * 100 = 1200 null traits. It is OPTIONAL for the headline
# results; run only if you want the null-calibration panel. Set NULL_GS_DIR to the
# directory of null .gs files.
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

NULL_GS_DIR="${NULL_GS_DIR:-${DATA_DIR}/gene_sets/gs_files_null}"
H5AD="${SUBSET_DIR}/TMS_FACS/TMS_FACS.h5ad"
COV="${SUBSET_DIR}/TMS_FACS/TMS_FACS.cov"
OUT="${RESULTS}/null/scdrs+_results"; mkdir -p "${OUT}"

if [[ ! -d "${NULL_GS_DIR}" ]]; then
  echo "[10] NULL_GS_DIR not found: ${NULL_GS_DIR}"
  echo "     Generate null gene sets first (see notebooks/02) or set NULL_GS_DIR. Skipping."
  exit 0
fi

families=(all highmean highvar overdispersed)
sizes=(100 500 1000)
for prefix in "${families[@]}"; do
  for ngene in "${sizes[@]}"; do
    traits=()
    for rep in $(seq 0 99); do traits+=("${prefix}_ngene${ngene}_rep${rep}"); done
    echo "[10] ${prefix}_ngene${ngene}: ${#traits[@]} null traits"
    run_scdrs_fm "${H5AD}" "${COV}" "${OUT}" "${NULL_GS_DIR}" "${traits[@]}" \
      --h5ad_species human --flag_filter --flag_raw_count --imputation magic
  done
done
echo "[10] Done -> ${OUT}"
