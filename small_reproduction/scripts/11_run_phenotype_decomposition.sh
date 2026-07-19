#!/usr/bin/env bash
# 11_run_phenotype_decomposition.sh  — gradient decomposition for nb11.
# Runs the UNMODIFIED scDRS-FM package function
# scdrs_fm.decompose_gradients.run_decomposition (CLI: run_decompose_gradients.py)
# to produce the {trait}.marg-marg.decomposition.tsv.gz files that
# notebooks/11_phenotype_analysis.ipynb reads.
#
# Config (user-approved reproduction settings):
#   * phenotype predictors = 24 functional signatures (microglia dataset -> 5 microglia sigs)
#   * outcome scores       = MAGIC-imputed
#   * cells                = ALL cells in each 10k subset
#
# Prereq: the T-cell/microglia phenotype scoring runs (scripts 04b/05b/07b/08b)
# and the immune ct scoring (13_make_ct_tables.sh) must have produced the
# *_full.gz control-augmented score tables this step consumes.
#
# Run all datasets:   bash scripts/11_run_phenotype_decomposition.sh
# One dataset:        bash scripts/11_run_phenotype_decomposition.sh soskic
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/_config.sh"
export REPO_ROOT SUBSET_DIR RESULTS SCDRS_FM_HOME
export LOCAL_TMP="${LOCAL_TMP:-/tmp}"

"${PYTHON}" "$(dirname "${BASH_SOURCE[0]}")/run_decomposition.py" "$@"
echo "[11] Decomposition done."
