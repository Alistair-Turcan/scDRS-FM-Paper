#!/usr/bin/env bash
# 00_make_10k_subsets.sh
# Subsample every real dataset to EXACTLY 10,000 cells (fixed seed=0), preserving
# raw counts and the full gene set, and align each covariate file to the subset.
# Mirrors the 10k-scale reproduction design (see docs/DATA.md).
#
# Requires: download/download_figshare_data.sh already run (populates data/extracted/).
# Output: data/subsets_10k/<Dataset>/{<h5ad>,<cov>} + subset_summary.json
#
# Run all datasets:      bash scripts/00_make_10k_subsets.sh
# Run one dataset:       bash scripts/00_make_10k_subsets.sh TMS_FACS
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"

export DATA_DIR SUBSET_DIR
export DATA_EXTRACTED="${DATA_EXTRACTED:-${DATA_DIR}/extracted}"
export LOCAL_TMP="${LOCAL_TMP:-/tmp/subset_tmp}"

"${PYTHON}" "$(dirname "${BASH_SOURCE[0]}")/make_10k_subsets.py" "${1:-}"
echo "[00] Subsets written to ${SUBSET_DIR}"
